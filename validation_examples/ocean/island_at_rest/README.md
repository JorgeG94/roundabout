# Island-at-rest — static land-mask quiescent trap

Headline correctness check for the ocean **static land-masking** foundation
(CHUNK A). Mirrors the seamount quiescent trap that caught three latent
dyn-core bugs: a domain that *must* stay at rest, where any motion is a bug.

## Setup (`island_at_rest.nml`)

- 32 × 32 closed basin @ 5 km (160 km square), 4 σ-layers.
- `topo_config = "island"`: flat 1000 m basin with a **central LAND square**
  (`b = 0 < LAND_DEPTH_THRESHOLD`). `slope_scale = 0.2` is reused as the land
  half-fraction → the middle 40 % of each axis is land.
- Uniform `T = 15 °C`, `S = 35 PSU` ⇒ density horizontally and vertically
  uniform ⇒ **zero APE**; the BPG over uniform ρ cancels exactly.
- f-plane `f = 0`, no wind, no surface flux, no sponges, closed walls.

## Expected outcome (GPU run — the orchestrator runs this, not Stage-6)

The static land mask makes every interior land face a hard, free-slip,
no-normal-flow wall (the 6 face metrics `dy_cu, idxCu, dxCu, dx_cv, idyCv,
dyCv` are zeroed at land faces; land T-cells are held at finite reference
state — `h_layer = H_VANISHED`, T/S held, velocities 0). With that:

- `max|u|`, `max|v|`, `max|η|` stay at **machine round-off** for the entire
  run — no spurious coastline pressure-gradient / flux / vorticity.
- **Mass conserved to ~1e-13**: the wet domain is closed (zero transport
  through every land face).

## Failure interpretation

- Energy / `max|u|` growing → a land-mask bug: an un-zeroed face metric, a
  missing velocity reset, or PPM `h_face` contamination by a land neighbour
  next to the coast.
- Note the velocity-reset and mirror-h pieces (C1/C2/C4) are **CHUNK B**;
  CHUNK A ships the metric-zeroing + masks + finite-land seeding. At rest
  with zero APE and zero forcing, CHUNK A alone should already hold the basin
  at round-off (nothing drives the un-reset velocities). A non-trivial
  growth here would indicate a CHUNK-A gap rather than a CHUNK-B refinement.
