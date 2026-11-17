# Preparing a runoff file for `&ocean_dataovr_nml lrunoff / frunoff`

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
`&ocean_forcing_nml enable_components = .true.` and `freshwater = "virtual"` or
`"mass"`. One file per tool call: concatenate years offline (a multi-file IAF
list is not supported yet).

## Limits

- Runoff files are refused together with an ice-shelf cavity.
- The coastal test wraps in i on a periodic mosaic and mirrors across a tripolar north fold.
  "Nearest coast" is by great-circle distance between cell centres, so it can cross a land
  bridge into the wrong basin: check the coastal cells that receive large rivers.
- The tool's wet test is `depth > --min-depth`; the model's wet mask may differ (a minimum
  depth, filled basins). Pass the same minimum depth, or a file with a value on a cell the model
  treats as land aborts the run at the read.
- Source cell areas use `R = 6.378e6 m`. A product built with `6.371e6 m` differs by about 0.2 %
  in the global total; the check is internal to the tool and cannot see it.
- A point poleward of the mosaic (a lon-lat mosaic ending at, say, 80N) is outside and is
  discarded.
