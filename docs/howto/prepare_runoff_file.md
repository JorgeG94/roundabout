# Preparing a runoff file for `[ocean.dataovr] lrunoff / frunoff`

The model reads runoff from a file that is ALREADY on its grid, with every
non-zero value on a wet cell (a value on a land cell aborts at the read: that
mass would vanish silently). `tools/om_runoff_to_model.f90` builds such a file
from a regular lon-lat product in kg m-2 s-1 (JRA55-do `friver` for liquid,
`licalvf` for frozen; Tsujino et al., 2018).

```bash
nvfortran -O2 $(nf-config --fflags) -o om_runoff_to_model tools/om_runoff_to_model.f90 $(nf-config --flibs)
./om_runoff_to_model ocean_hgrid.nc topog.nc friver_year.nc friver runoff_liquid.nc --msub 4
./om_runoff_to_model ocean_hgrid.nc topog.nc licalvf_year.nc licalvf runoff_frozen.nc
./om_runoff_to_model --selftest        # synthetic mosaic, no data needed
```

What it does: every source cell is turned into a mass rate with its exact
spherical area and cut into `msub x msub` equal-area sub-cells; each sub-cell
hands its share to the model T cell that contains it (nearest-centre
quadrilateral test; the nearest centre when the point lies within two cell
radii of it). A sub-cell farther than that from every model cell is OUTSIDE the
model domain (a regional grid, a polar hole): its share is DISCARDED, and the
discarded total (kg/s) and sub-cell count are printed per record. Use
`--global` on a global grid: any discard is then fatal. Shares that land on
model land cells (`depth <= --min-depth`, default 0, in `topog.nc`) move to the
nearest wet coastal cell. It prints the totals for every record and aborts if
in-domain + discarded differs from the source by more than 1e-12 relative or if
any land cell is non-zero. It also reports, per record, the mass moved from land to the coast, and at the end
how many land cells carried mass and the farthest centre-to-centre move in km: a large
value is a river put on the wrong coast, so look at it. Missing values (fill, NaN, above 1e10) count as zero and are
reported. The output keeps the source time axis and its `units`; point
`lrunoff_file` / `frunoff_file` at it with `lrunoff_var = "runoff"` and the
`time_mode` you want (`linear`, `cyclic`, `static`). It needs
`[ocean.forcing] enable_components = true` and `freshwater = "virtual"` or
`"mass"`. One file per tool call: concatenate years offline (a multi-file IAF
list is not supported yet).

## Spreading the rivers

A river is a point source in the product (one or two 0.25-degree cells) and, put on one 1-degree
coastal cell, it is a large flux into a thin column: the Amazon of the synthetic test file falls on
one 10 m minimum-depth cell at up to 5.4 m/day of fresh water. Real coupled systems therefore
distribute river discharge over the coastal ocean instead of releasing it into one cell. The tool
does it offline, after the remap and the land-to-coast move, and the model reader is untouched.

```bash
./om_runoff_to_model ... --spread-radius-km 500 --spread-efold-km 150      # the defaults
./om_runoff_to_model ... --spread-radius-km 0                              # no spreading: the unspread file, bit for bit
./om_runoff_to_model ... --spread-max-rate 0.5 --spread-hard-km 1500       # widen sources that stay above 0.5 m/day
```

- Every cell's runoff mass rate (kg/s) is shared over the WET cells reachable from it through
  wet cells only (4-neighbour flood fill, wrapping in i on a periodic mosaic and mirrored across a
  tripolar fold) within `--spread-radius-km` great-circle distance, with weights
  `exp(-d/--spread-efold-km)` normalised to sum to 1. Nothing is lost or created: the global kg/s
  before and after is asserted equal per record (1e-12 relative), and the extra report lines give
  the number of cells with runoff and the peak rate in m/day before and after.
- Because the fill only walks wet cells, runoff cannot cross a land bridge (Panama, a strait one
  cell wide) into the other basin, and a closed one-cell bay keeps its own runoff. This closes
  the "nearest coast can cross a land bridge" limit for the spreading step; the earlier
  land-to-coast move (which only places mass that fell on land) still picks the nearest coastal
  cell by great-circle distance.
- Defaults 500 km and 150 km: on the 1-degree tripolar mosaic they take the Amazon from 5.4 to
  0.49 m/day over 67 cells, and the Congo from 0.54 to 0.05 (Ganges 0.45 to 0.07, Orinoco 0.87 to 0.08;
  synthetic test file). 300/100 gave 1.0 m/day for the Amazon, 500/250 is nearly flat over
  the disc (0.25 m/day) and puts mass 500 km offshore at the same weight as at the mouth.
  Pick the e-folding near the model's shelf width; the radius should be 3 to 4 e-foldings.
- `--spread-max-rate M_PER_DAY` (default off, fresh-water m/day at 1000 kg/m3): after spreading, a
  source whose targets still exceed the rate is widened (radius x1.5 per level, at most 8 levels and
  `--spread-hard-km`, default 1500) and the loop repeats; cells still above it are reported.
- Cost: the operator is linear and static, built once per source cell (a short list of target
  cell and weight), then applied per record. On 360 x 320 the whole run is dominated by the
  point-map step (about one minute for a 0.25-degree source); the operator adds a second or
  so. At 1440 x 1080 (0.25 degree) a 500 km disc is about 15 000 cells, so a few thousand
  source cells need a few 1e7 list entries (hundreds of MB); reduce the radius or e-folding
  there. Not measured at that size.

## Limits

- Runoff files are refused together with an ice-shelf cavity.
- The coastal test wraps in i on a periodic mosaic and mirrors across a tripolar north fold.
  "Nearest coast" (the move of mass that fell on land) is by great-circle distance between cell
  centres, so it can still pick the wrong basin across a land bridge: check the redirect report.
  The spreading step does not cross land.
- The tool's wet test is `depth > --min-depth`; the model's wet mask may differ (a minimum
  depth, filled basins). Pass the same minimum depth, or a file with a value on a cell the model
  treats as land aborts the run at the read.
- Source cell areas use `R = 6.378e6 m`. A product built with `6.371e6 m` differs by about 0.2 %
  in the global total; the check is internal to the tool and cannot see it.
- A point poleward of the mosaic (a lon-lat mosaic ending at, say, 80N) is outside and is
  discarded.
