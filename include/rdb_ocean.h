/* Ocean C ABI — lifecycle (P1) + accessors and the error ring (P2).
 *
 * Hand-written; no C header ever existed for the pre-carve-out FFI (ctypes
 * was the only client there — see tmp_local_artifacts/python_ffi_scope/
 * 03_removed_ffi_precedent.md, section 3). This one exists so a non-Python
 * C/C++ caller has a real prototype set to compile against, and so the
 * Fortran `bind(c, name=...)` signatures in src/api/rdb_ocean_api.F90 have
 * a single documented source of truth to be checked against.
 *
 * P1 scope: create / step / destroy, plus enough introspection (time, step
 * count, grid shape, one scalar mass diagnostic, working precision) to prove
 * a solver actually advanced.
 *
 * P2 scope (this revision): the error ring (rdb_ocean_last_error /
 * rdb_flush_logs), rdb_ocean_refresh_host, ~17 raw-pointer state
 * getters, tracer-by-name access, 8 narrow setters, and a kinetic-energy
 * scalar. See docs/ocean_python_api_plan.md S3 and
 * tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md D3/D4
 * for the D<->H contract these follow: getters return (ptr, extents,
 * generation) WITHOUT refreshing the host copy themselves — a caller must
 * compare `generation` against what it last saw and call
 * rdb_ocean_refresh_host() before dereferencing a stale pointer.
 * `generation` is the handle's outer-step count at the moment of the call.
 * Raw pointers point at FULL (ghost-inclusive) arrays; ghost width comes
 * from rdb_ocean_get_grid_info(). Tracer arrays are the raw `h*Tr` store
 * (NOT concentration) — dividing by `h_layer` (with an H_VANISHED guard) is
 * the caller's job, not this library's.
 *
 * P2.5 scope (this revision): pre-create GEOMETRY INJECTION
 * (Oceananigans-style — docs/ocean_python_api_plan.md S5b) via a
 * two-phase create: rdb_ocean_create_pending() builds + validates a
 * config and allocates a handle WITHOUT running setup or mapping the
 * device, then any of rdb_ocean_stage_bathymetry() /
 * _stage_metrics() / _stage_topology() may be called (in any order, any
 * subset) to inject an in-memory bathymetry array, MOM6-style supergrid
 * metrics, and/or grid topology (per-dimension periodicity), and finally
 * rdb_ocean_create_finalize() runs setup (consuming whatever was
 * staged) through device mapping — exactly like the tail of
 * rdb_ocean_create_from_string(), which remains the normal one-call
 * path for a caller with nothing to inject and is entirely unaffected by
 * this phase. rdb_ocean_required_halo() is a stateless query (no
 * handle) folding the library's scattered minimum-nghost rules so a
 * caller sizing a grid's halo agrees with what create_finalize() will
 * itself enforce.
 *
 * Link against librdb_core.so (RDB_BUILD_SHARED=ON) or the static
 * core_rdb archive; both export these symbols (the API sources compile
 * into core_objs, which backs both).
 */
#ifndef RDB_OCEAN_H
#define RDB_OCEAN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle. Always initialise to NULL before create() and check for
 * NULL after every call — every entry point below leaves the handle
 * untouched (or nulls it, for destroy) on failure. */
typedef void *rdb_ocean_handle;

/* Status codes returned by every entry point below.
 *
 * 0           success.
 * 1-5         config/setup path (rdb_ocean_status::OCEAN_STATUS_ERR_*,
 *             shared with the Fortran-only namelist-file create path —
 *             P0 of the API plan).
 * 10-19       handle-lifecycle + P2-accessor path (this API only); 15/16
 *             are P2.5 geometry injection (bathymetry sign, pending-state).
 * 20-22       restart-resume mismatch (schema/decomp/grid) — see
 *             rdb_ocean_status.F90.
 *
 * The numeric values are load-bearing: they are the literal
 * `integer, parameter` values in src/core/ocean/state/rdb_ocean_status.F90.
 * Keep this block in sync with that file if it ever grows.
 *
 * A non-zero status here is only the STAGE. For the SPECIFIC reason, call
 * rdb_ocean_last_error(0, buf, cap) immediately after — see below. */
#define RDB_OCEAN_OK 0
#define RDB_OCEAN_ERR_CONFIG_PARSE 1
#define RDB_OCEAN_ERR_CONFIG_VALIDATE 2
#define RDB_OCEAN_ERR_SETUP 3
#define RDB_OCEAN_ERR_IC_SEED 4
#define RDB_OCEAN_ERR_IO 5
#define RDB_OCEAN_ERR_BAD_HANDLE 10
#define RDB_OCEAN_ERR_ALREADY_EXISTS 11
#define RDB_OCEAN_ERR_NOT_INITIALISED 12
#define RDB_OCEAN_ERR_NOT_FOUND 13
#define RDB_OCEAN_ERR_BAD_SHAPE 14
#define RDB_OCEAN_ERR_BATHYMETRY_SIGN 15
#define RDB_OCEAN_ERR_NOT_PENDING 16
#define RDB_OCEAN_ERR_RESTART_SCHEMA 20
#define RDB_OCEAN_ERR_RESTART_DECOMP 21
#define RDB_OCEAN_ERR_RESTART_GRID 22

/* ---- Lifecycle ---------------------------------------------------- */

/* Build a config from an in-memory namelist string (newline-separated
 * records, no filesystem touch), set up the ocean dyn-core, map it onto the
 * device, and hand back a handle in *handle_out.
 *
 * Fails with RDB_OCEAN_ERR_ALREADY_EXISTS if an ocean handle is already
 * live in this process — multi-instance is not supported (see the "Single
 * live ocean handle" note in src/api/rdb_ocean_api.F90). On any failure,
 * *handle_out is set to NULL.
 *
 * nml_text need not be NUL-terminated; nml_len is authoritative. */
int rdb_ocean_create_from_string(const char *nml_text, int nml_len,
                                    rdb_ocean_handle *handle_out);

/* Advance n_steps fixed-dt outer steps. n_steps <= 0 is a successful no-op. */
int rdb_ocean_step(rdb_ocean_handle handle, int n_steps);

/* Idempotent: NULL, an already-destroyed handle, or garbage is a no-op
 * success — safe for a Python __del__ (or any caller) to invoke blind.
 * On a live handle: unwinds device residency, releases host allocations,
 * frees the handle, and clears the single-live-handle guard. *handle is set
 * to NULL on return either way. */
int rdb_ocean_destroy(rdb_ocean_handle *handle);

/* ---- P2.5: pre-create geometry injection ----------------------------
 * Two-phase create: create_pending() -> any subset of stage_*() (any
 * order) -> create_finalize(). Geometry MUST land before create_finalize()
 * (which runs engine setup through device mapping) — there is no
 * post-create equivalent to stage_metrics()/stage_topology(), and
 * rdb_ocean_set_bathymetry() (below, P2) is a NARROWER mid-run
 * perturbation, not a re-seed (it does not redo land masking / the ALE
 * z_ref table / n_inner derivation the way create_finalize()'s pipeline
 * does).
 */

/* Sign conventions for rdb_ocean_stage_bathymetry()'s `convention`
 * argument. No default — required and checked: an unrecognised value is
 * RDB_OCEAN_ERR_SETUP, a value that normalises to zero wet cells is
 * RDB_OCEAN_ERR_BATHYMETRY_SIGN. Roundabout's `%barotropic%b` is a
 * POSITIVE-DOWN depth (a 4000 m-deep cell is +4000); GEBCO/ETOPO/
 * Oceananigans GridFittedBottom ship a NEGATIVE-DOWN height (the same
 * cell is -4000). */
#define RDB_BATHY_CONVENTION_DEPTH_POSITIVE_DOWN 1
#define RDB_BATHY_CONVENTION_HEIGHT_POSITIVE_UP 2

/* Phase 1 of 2: build + validate a config and allocate a handle exactly
 * like rdb_ocean_create_from_string(), but stop there — does NOT run
 * setup or map the device. Claims the single-live-handle guard
 * immediately. Fails with RDB_OCEAN_ERR_ALREADY_EXISTS under the same
 * condition as rdb_ocean_create_from_string(). On any failure,
 * *handle_out is set to NULL (the same "destroy the whole handle" F9
 * contract). */
int rdb_ocean_create_pending(const char *nml_text, int nml_len,
                                rdb_ocean_handle *handle_out);

/* Stage an interior-sized (nx_p, ny_p) bathymetry array on a PENDING
 * handle (from create_pending(), before create_finalize()). Consumed by
 * setup inside create_finalize(): sign-normalised to positive-down depth
 * per `convention`, wet-fraction validated (RDB_OCEAN_ERR_BATHYMETRY_SIGN
 * on zero wet cells — see the convention macros above), written into the
 * interior, and ghost-filled by constant extrapolation — this OVERRIDES
 * whatever `topo_config` the namelist specifies. Shape is validated
 * against the grid inside create_finalize() (the grid does not exist yet
 * at this call), so a mismatch surfaces there as
 * RDB_OCEAN_ERR_BAD_SHAPE, not here.
 * Requires a PENDING handle: RDB_OCEAN_ERR_NOT_PENDING otherwise
 * (including on an already-finalised handle). */
int rdb_ocean_stage_bathymetry(rdb_ocean_handle handle,
                                  const double *b_data, int nx_p, int ny_p,
                                  int convention);

/* Stage in-memory MOM6-style supergrid metric arrays on a PENDING handle
 * — the same layout rdb_ocean_metrics::metrics_assemble_from_supergrid_arrays
 * consumes: x/y (degrees, shape (nxp,nyp)), dx (m, shape (nx,nyp)), dy (m,
 * shape (nxp,ny)), area (m^2, shape (nx,ny)), where nxp=2*nx_phys+1,
 * nyp=2*ny_phys+1, nx=2*nx_phys, ny=2*ny_phys. Consumed inside
 * create_finalize(), bypassing `&ocean_grid_nml grid_config` entirely —
 * this is the SAME assembler the NetCDF supergrid reader and the analytic
 * tripolar generator use, so a mosaic file and these arrays (built from
 * the same node/segment data) produce identical metrics. Shape is
 * validated against the grid inside create_finalize().
 * Requires a PENDING handle: RDB_OCEAN_ERR_NOT_PENDING otherwise. */
int rdb_ocean_stage_metrics(rdb_ocean_handle handle, const double *x,
                               const double *y, const double *dx,
                               const double *dy, const double *area, int nxp,
                               int nyp, int nx, int ny);

/* Stage Oceananigans-style grid topology (per-dimension periodicity —
 * "the grid owns periodicity", not the per-edge &ocean_bc_nml tags) on a
 * PENDING handle. Nonzero = periodic on that axis. Consumed inside
 * create_finalize() AFTER the namelist edge tags are parsed, OVERRIDING
 * whatever they derived for periodic_x/periodic_y and back-filling the
 * edge tags on every axis marked periodic (both edges together — a
 * west/east or south/north mismatch is not representable through this
 * call). An axis NOT marked periodic here keeps whatever physical BC the
 * namelist gave it (WALL/OPEN/...).
 * Requires a PENDING handle: RDB_OCEAN_ERR_NOT_PENDING otherwise. */
int rdb_ocean_stage_topology(rdb_ocean_handle handle, int periodic_x,
                                int periodic_y);

/* Phase 2 of 2: complete a handle started by create_pending() (optionally
 * staged with any of the stage_*() calls above in between) — runs setup
 * (consuming whatever was staged) through device mapping, exactly like
 * the tail of rdb_ocean_create_from_string(). On success, *handle is
 * UNCHANGED (same pointer value, now fully initialised) — use it with
 * every ordinary lifecycle/query/getter/setter call below. On failure,
 * the WHOLE handle is destroyed (F9) and *handle is set to NULL — same
 * out-null convention as rdb_ocean_destroy(), since (unlike
 * create_from_string(), whose caller only ever sees a handle on success)
 * a caller here already holds a handle value going in.
 * Requires a PENDING handle: RDB_OCEAN_ERR_NOT_PENDING on a handle
 * that was never created via create_pending(), or was already finalised. */
int rdb_ocean_create_finalize(rdb_ocean_handle *handle);

/* Stateless (no handle — like rdb_working_precision()): the minimum
 * nghost for a given scheme/topology selection, folding the library's six
 * scattered minimum-nghost rules (PV-advection WENO order, tracer
 * reconstruction WENO order, periodic topology, tripolar fold,
 * MPI-decomposed run, kappa-shear vertex form) behind ONE query, so a
 * caller sizing a grid's halo before create() agrees with what
 * create_finalize() will itself enforce. A zero-length scheme string
 * means "unset" (baseline 2, the centered/PPM/non-periodic/single-rank
 * floor). Nonzero = true for the four topology/scheme flags. Always
 * succeeds (returns a halo width, never a status code). */
int rdb_ocean_required_halo(const char *pv_adv_scheme,
                               int pv_adv_scheme_len,
                               const char *tracer_recon,
                               int tracer_recon_len, int periodic,
                               int tripolar_fold, int decomposed,
                               int kappa_shear_at_vertex);

/* ---- Queries -------------------------------------------------------
 * Every query below returns RDB_OCEAN_ERR_BAD_HANDLE on a bad handle and
 * RDB_OCEAN_ERR_NOT_INITIALISED on a handle that never finished create().
 */

/* Current simulation time, in seconds. */
int rdb_ocean_get_time(rdb_ocean_handle handle, double *t_out);

/* Outer-step counter. */
int rdb_ocean_get_step_count(rdb_ocean_handle handle, int *step_out);

/* Physical (ghost-excluded) grid shape + ghost width. */
int rdb_ocean_get_grid_info(rdb_ocean_handle handle, int *nx, int *ny,
                               int *nz, int *nghost);

/* Total water mass over the physical domain (sum(h_layer) * dx * dy *
 * reference density) — the one scalar diagnostic exposed at this phase, so
 * a caller can prove the solver advanced (and, in a closed basin, that it is
 * conserving mass) without any state-array accessor (those are P2). */
int rdb_ocean_get_total_mass(rdb_ocean_handle handle, double *m_out);

/* Bytes per working-precision float (4 or 8) — no handle required; resolve
 * float32 vs float64 once at library load instead of guessing. */
int rdb_working_precision(void);

/* ---- P2: the error ring --------------------------------------------
 * A ring, not a single slot: several create()-path failures thread a
 * specific inner reason up through an outer wrapper, and a single slot
 * would keep only the wrapper's generic message. Index 0 is the most
 * recently pushed (usually the deepest/most specific).
 */

/* Copy the ring message at ring-relative index idx (0 = most recent) into
 * buf, snprintf-style: writes up to cap bytes, NUL-terminates if room
 * remains, and returns the FULL trimmed message length (a return value
 * >= cap means truncation). Returns 0 (buf untouched) if idx is out of
 * range or cap <= 0 — NOT an error code, just "nothing there". */
int rdb_ocean_last_error(int idx, char *buf, int cap);

/* Flush the buffered log stream (stdout) so it is at least ordered
 * relative to the caller's own output. The ring above, not this stream,
 * is the channel to trust for the specific failure reason. */
void rdb_flush_logs(void);

/* ---- P2: D<->H refresh ----------------------------------------------
 * Every getter below returns (ptr, extents..., generation) WITHOUT
 * refreshing the host copy. generation is the handle's outer-step count
 * at the call. Compare it against what you last saw; on a mismatch, call
 * this before dereferencing ptr. Idempotent / cheap when already current.
 */
int rdb_ocean_refresh_host(rdb_ocean_handle handle);

/* ---- P2: getters ------------------------------------------------------
 * Raw pointers into the live host arrays — FULL extent, ghosts included
 * (see rdb_ocean_get_grid_info for nghost). Never refreshed by the
 * getter itself; see rdb_ocean_refresh_host above. All real values are
 * the build's working-precision float (see rdb_working_precision) —
 * NOT necessarily double.
 */

/* Layer thickness (m), cell-centred: (nx_total, ny_total, nz_ml). */
int rdb_ocean_get_h_layer_ptr(rdb_ocean_handle handle, void **ptr,
                                 int *nx, int *ny, int *nz, int *gen);
/* West-face x-velocity (m/s): (nx_total+1, ny_total, nz_ml). */
int rdb_ocean_get_u_face_x_layer_ptr(rdb_ocean_handle handle, void **ptr,
                                        int *nx, int *ny, int *nz, int *gen);
/* South-face y-velocity (m/s): (nx_total, ny_total+1, nz_ml). */
int rdb_ocean_get_v_face_y_layer_ptr(rdb_ocean_handle handle, void **ptr,
                                        int *nx, int *ny, int *nz, int *gen);
/* West-face x transport h*u (m^2/s), same shape as u_face_x_layer. */
int rdb_ocean_get_hu_ptr(rdb_ocean_handle handle, void **ptr, int *nx,
                            int *ny, int *nz, int *gen);
/* South-face y transport h*v (m^2/s), same shape as v_face_y_layer. */
int rdb_ocean_get_hv_ptr(rdb_ocean_handle handle, void **ptr, int *nx,
                            int *ny, int *nz, int *gen);
/* Vertical velocity at layer interfaces (m/s): (nx_total, ny_total, nz_ml+1),
 * k=1 bed .. k=nz_ml+1 surface. */
int rdb_ocean_get_w_interface_ptr(rdb_ocean_handle handle, void **ptr,
                                     int *nx, int *ny, int *nz, int *gen);
/* In-situ density (kg/m^3), same shape as h_layer. */
int rdb_ocean_get_rho_layer_ptr(rdb_ocean_handle handle, void **ptr,
                                   int *nx, int *ny, int *nz, int *gen);
/* Bathymetry (m, POSITIVE DOWN): (nx_total, ny_total). */
int rdb_ocean_get_b_ptr(rdb_ocean_handle handle, void **ptr, int *nx,
                           int *ny, int *gen);
/* Sea-surface height (m), diagnostic: (nx_total, ny_total). */
int rdb_ocean_get_bt_eta_ptr(rdb_ocean_handle handle, void **ptr,
                                int *nx, int *ny, int *gen);
/* East-face wind stress (N/m^2): (nx_total+1, ny_total). */
int rdb_ocean_get_tau_x_ptr(rdb_ocean_handle handle, void **ptr,
                               int *nx, int *ny, int *gen);
/* North-face wind stress (N/m^2): (nx_total, ny_total+1). */
int rdb_ocean_get_tau_y_ptr(rdb_ocean_handle handle, void **ptr,
                               int *nx, int *ny, int *gen);
/* Net surface heat flux (W/m^2): (nx_total, ny_total). */
int rdb_ocean_get_q_heat_ptr(rdb_ocean_handle handle, void **ptr,
                                int *nx, int *ny, int *gen);
/* Net surface salt flux: (nx_total, ny_total). */
int rdb_ocean_get_q_salt_ptr(rdb_ocean_handle handle, void **ptr,
                                int *nx, int *ny, int *gen);
/* Vertical viscosity (m^2/s): (nx_total, ny_total, nz_ml+1), interfaces. */
int rdb_ocean_get_kv_ptr(rdb_ocean_handle handle, void **ptr, int *nx,
                            int *ny, int *nz, int *gen);
/* Vertical heat diffusivity (m^2/s): (nx_total, ny_total, nz_ml+1). */
int rdb_ocean_get_kt_ptr(rdb_ocean_handle handle, void **ptr, int *nx,
                            int *ny, int *nz, int *gen);
/* Vertical salt (+ passive tracer) diffusivity (m^2/s): same shape as kv. */
int rdb_ocean_get_ks_ptr(rdb_ocean_handle handle, void **ptr, int *nx,
                            int *ny, int *nz, int *gen);
/* Wet mask at T points (1=wet, 0=land), read-only: (nx_total, ny_total). */
int rdb_ocean_get_wet_t_ptr(rdb_ocean_handle handle, void **ptr,
                               int *nx, int *ny, int *gen);

/* ---- P2: tracers by name --------------------------------------------- */

/* Number of registered tracers (S, T, + any passive tracers). */
int rdb_ocean_get_tracer_count(rdb_ocean_handle handle, int *count_out);

/* Name of the tracer at registry index idx (0-based), snprintf-style same
 * contract as rdb_ocean_last_error. */
int rdb_ocean_list_tracers(rdb_ocean_handle handle, int idx, char *buf,
                              int cap);

/* Raw h*Tr store (NOT concentration) for the tracer named `name`
 * (name_len bytes, need not be NUL-terminated): (nx_total, ny_total,
 * nz_ml). RDB_OCEAN_ERR_NOT_FOUND if no tracer matches. */
int rdb_ocean_get_tracer_ptr(rdb_ocean_handle handle, const char *name,
                                int name_len, void **ptr, int *nx, int *ny,
                                int *nz, int *gen);

/* ---- P2: narrow setters -----------------------------------------------
 * Each pushes ONLY the array(s) it mutates to device — never a broad
 * "push everything" that would rewind live device prognostics to a stale
 * host snapshot if called mid-run with no intervening read. All data
 * arrays here are double (converted to the build's working precision
 * internally); shapes match the corresponding getter above, restricted to
 * the PHYSICAL interior (no ghosts): pass nx_p = physical nx (etc.), and
 * face arrays get the "+1" in the staggered dimension as noted.
 */

/* Overwrite h_layer on the physical interior (nx_p, ny_p, nz_p) and
 * recompute rho_layer from it. */
int rdb_ocean_set_h(rdb_ocean_handle handle, const double *h_data,
                       int nx_p, int ny_p, int nz_p);
/* Overwrite u_face_x_layer: u_data is (nx_p+1, ny_p, nz_p). */
int rdb_ocean_set_u(rdb_ocean_handle handle, const double *u_data,
                       int nx_p, int ny_p, int nz_p);
/* Overwrite v_face_y_layer: v_data is (nx_p, ny_p+1, nz_p). */
int rdb_ocean_set_v(rdb_ocean_handle handle, const double *v_data,
                       int nx_p, int ny_p, int nz_p);
/* Overwrite bathymetry b (nx_p, ny_p), positive down; fills ghosts,
 * re-wraps the periodic/fold seam, re-derives bt_H_ref. */
int rdb_ocean_set_bathymetry(rdb_ocean_handle handle,
                                const double *b_data, int nx_p, int ny_p);
/* Overwrite wind stress: taux_data is (nx_p+1, ny_p), tauy_data is
 * (nx_p, ny_p+1). */
int rdb_ocean_set_wind(rdb_ocean_handle handle, const double *taux_data,
                          const double *tauy_data, int nx_p, int ny_p);
/* Overwrite net surface heat flux (nx_p, ny_p). */
int rdb_ocean_set_heat_flux(rdb_ocean_handle handle,
                               const double *q_data, int nx_p, int ny_p);
/* Overwrite net surface salt flux (nx_p, ny_p). */
int rdb_ocean_set_salt_flux(rdb_ocean_handle handle,
                               const double *q_data, int nx_p, int ny_p);
/* Set the tracer named `name` (concentration, its own units — degC/PSU/
 * etc, NOT h*Tr) on the physical interior (nx_p, ny_p, nz_p).
 * RDB_OCEAN_ERR_NOT_FOUND if no tracer matches. */
int rdb_ocean_set_tracer(rdb_ocean_handle handle, const char *name,
                            int name_len, const double *data, int nx_p,
                            int ny_p, int nz_p);

/* ---- P2: scalar diagnostics -------------------------------------------
 * Unlike the getters above, this refreshes the host itself — it hands
 * back a number, not a pointer a caller could otherwise defer syncing.
 */
int rdb_ocean_get_kinetic_energy(rdb_ocean_handle handle,
                                    double *ke_out);

/* ---- P7: diagnostics discoverability -----------------------------------
 * Two STATIC catalogs (no handle needed -- "what CAN I request", callable
 * before create()) plus the LIVE registered set on one handle ("what IS
 * this instance actually emitting"). Name getters share the
 * rdb_ocean_list_tracers snprintf-style contract: copy up to cap bytes,
 * NUL-terminate if room remains, return the FULL trimmed name length; 0
 * (buf untouched) on an out-of-range index or cap <= 0.
 */

/* Number of names in the DERIVED diagnostic catalog (vorticity_z, ke_total,
 * mld_density, ... -- opt-in via &ocean_diag_nml diags). */
int rdb_ocean_derived_catalog_size(void);
/* Name of derived-catalog entry idx (0-based). */
int rdb_ocean_derived_catalog_name(int idx, char *buf, int cap);

/* Number of names in the CANONICAL diagnostic catalog (SSH, temperature,
 * salinity, u, v, KE, ... -- some entries are gated by another namelist
 * group; this is what register_default_diags MAY register). */
int rdb_ocean_canonical_catalog_size(void);
/* Name of canonical-catalog entry idx (0-based). */
int rdb_ocean_canonical_catalog_name(int idx, char *buf, int cap);

/* Number of diagnostics REGISTERED on this live instance right now
 * (canonical + derived + anything &ocean_diag_nml diags added). */
int rdb_ocean_get_diag_count(rdb_ocean_handle handle, int *count_out);
/* Name of the registered diagnostic at index idx (0-based). */
int rdb_ocean_list_diags(rdb_ocean_handle handle, int idx, char *buf,
                            int cap);

/* ---- P7: in-memory diagnostic access -----------------------------------
 * The diag manager's OWN output_buffer for a registered diagnostic, by
 * name -- no NetCDF round-trip. output_buffer is pulled host-ward
 * unconditionally at every cadence fire, so (unlike the P2 state getters)
 * this needs no separate refresh call: by the time rdb_ocean_step
 * returns, every diagnostic that fired is already host-current. gen is a
 * per-variable fire count (not the global step count) -- a diagnostic on
 * a multi-hour cadence should not look "stale" between fires.
 */

/* Shape: (nx_total, ny_total, nz_ml) for a layered var at the default
 * LAYER vgrid, (nx_total, ny_total, 1) for a 2D var, or the remapped
 * shape when output_vgrid != LAYER. RDB_OCEAN_ERR_NOT_FOUND if `name`
 * is not currently REGISTERED on this instance (see
 * rdb_ocean_get_diag_count / rdb_ocean_list_diags). */
int rdb_ocean_get_diagnostic_ptr(rdb_ocean_handle handle,
                                    const char *name, int name_len,
                                    void **ptr, int *nx, int *ny, int *nz,
                                    int *gen);

#ifdef __cplusplus
}
#endif

#endif /* RDB_OCEAN_H */
