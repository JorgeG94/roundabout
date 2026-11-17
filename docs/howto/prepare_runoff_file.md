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
quadrilateral test; nearest centre otherwise), so the global total is preserved
exactly. Shares that land on model land cells (`depth <= 0` in `topog.nc`) move
to the nearest wet coastal cell. It prints the total before and after for every
record and aborts if they differ by more than 1e-12 relative or if any land
cell is non-zero. Missing values (fill, NaN, above 1e10) count as zero and are
reported. The output keeps the source time axis and its `units`; point
`lrunoff_file` / `frunoff_file` at it with `lrunoff_var = "runoff"` and the
`time_mode` you want (`linear`, `cyclic`, `static`). It needs
`&ocean_forcing_nml enable_components = .true.` and `freshwater = "virtual"` or
`"mass"`. One file per tool call: concatenate years offline (a multi-file IAF
list is not supported yet).
