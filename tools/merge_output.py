#!/usr/bin/env python3
"""Merge per-rank Roundabout output NetCDF files into a single global file.

Roundabout writes one `<prefix>_rank_NNNNNN.nc` per MPI rank — no online gather.
This tool stitches them back into a global file as if the run had not been
partitioned. Both grid backends are supported and auto-detected:

  Structured     — global attrs `i_start`, `j_start`, `nx_global`, `ny_global`,
                   `nx_local`, `ny_local` define each rank's window into the
                   global `(y, x)` array. Coordinate variables `x`/`y` already
                   carry absolute positions.
  Unstructured   — per-face int64 variable `cell_global_id` maps each owned
                   triangle to its global cell index. Each rank writes only its
                   owned faces. Mesh topology (`mesh_node_x/y`,
                   `mesh_face_nodes`) is per-rank and uses local node indexing;
                   we don't try to stitch nodes (would need coord-based dedup
                   or the original UGRID file). Pass `--mesh-from <ugrid.nc>`
                   to copy the topology from the source UGRID input.

  Ocean diag     — per-rank files produced by the ocean diag manager.  Each
                   file has global attrs `i_start`, `j_start`, `nx_global`,
                   `ny_global`, `nx_local`, `ny_local`, `nghost`.  The spatial
                   arrays include `nghost` ghost cells on each side; the merge
                   trims them before stitching.  Time dims are per-variable
                   (`time_<name>`) and may be unlimited; z dims are per-variable
                   (`z_<name>`).  There are no `x`/`y` coordinate variables.

Usage
-----

    python tools/merge_output.py <run_dir> <merged.nc>
    python tools/merge_output.py <run_dir> <merged.nc> --prefix seamount_bench_full
    python tools/merge_output.py <run_dir> <merged.nc> --mesh-from input.nc
    python tools/merge_output.py --grid structured <run_dir> <merged.nc>

The default for unstructured is to skip topology and warn. Time-varying and
static face-centred fields are always merged.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings
from typing import List, Tuple

import netCDF4 as nc
import numpy as np


# ------------------------------------------------------------------ utilities


def find_rank_files(run_dir: str, prefix: str = "output") -> List[str]:
    pattern = os.path.join(run_dir, f"{prefix}_rank_*.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"error: no files matched '{pattern}'")
    return files


def detect_grid(path: str) -> str:
    with nc.Dataset(path) as ds:
        if "cell_global_id" in ds.variables:
            return "unstructured"
        if "i_start" in ds.ncattrs() and "nx_global" in ds.ncattrs():
            # Ocean diag files have nghost and no x/y coordinate variables.
            if "nghost" in ds.ncattrs() and "x" not in ds.variables:
                return "ocean_diag"
            return "structured"
    sys.exit(f"error: cannot detect grid type from {path} "
             f"(no cell_global_id var, no i_start/nx_global global attrs)")


def copy_attrs(src, dst, skip=()) -> None:
    for a in src.ncattrs():
        if a in skip:
            continue
        dst.setncattr(a, src.getncattr(a))


def make_var(out: nc.Dataset, src_var: nc.Variable, dims: Tuple[str, ...]) -> nc.Variable:
    """Create `src_var` in `out` with `dims`, preserving fill, deflate, attrs."""
    fill = getattr(src_var, "_FillValue", None)
    filt = src_var.filters() or {}
    chunks = None  # let netCDF4 auto-chunk for the global shape
    v = out.createVariable(
        src_var.name, src_var.dtype, dims,
        zlib=filt.get("zlib", False),
        complevel=filt.get("complevel", 0),
        shuffle=filt.get("shuffle", False),
        fill_value=fill,
        chunksizes=chunks,
    )
    copy_attrs(src_var, v, skip=("_FillValue",))
    return v


# --------------------------------------------------------------- structured


def merge_structured(rank_files: List[str], out_path: str) -> None:
    with nc.Dataset(rank_files[0]) as ref:
        nx_global = int(ref.getncattr("nx_global"))
        ny_global = int(ref.getncattr("ny_global"))
        ref_attrs = {a: ref.getncattr(a) for a in ref.ncattrs()}
        ref_dims = {d: len(ref.dimensions[d]) for d in ref.dimensions}
        ref_unlim = {d for d in ref.dimensions if ref.dimensions[d].isunlimited()}
        ref_vars = list(ref.variables.keys())

    print(f"[structured] {len(rank_files)} ranks → {nx_global} x {ny_global}")

    skip_attrs = {"i_start", "j_start", "nx_local", "ny_local"}

    with nc.Dataset(out_path, "w") as out:
        # Global attrs (drop per-rank window keys)
        for k, v in ref_attrs.items():
            if k in skip_attrs:
                continue
            out.setncattr(k, v)

        # Dimensions: rewrite x and y to global; pass through everything else.
        for d, n in ref_dims.items():
            if d == "x":
                out.createDimension(d, nx_global)
            elif d == "y":
                out.createDimension(d, ny_global)
            elif d in ref_unlim:
                out.createDimension(d, None)
            else:
                out.createDimension(d, n)

        # Define variables with global shape.
        with nc.Dataset(rank_files[0]) as ref:
            for name in ref_vars:
                make_var(out, ref.variables[name], ref.variables[name].dimensions)

        # ------- 1D global coordinates: stitch from every rank's slice -------
        x_global = np.zeros(nx_global, dtype=np.float64)
        y_global = np.zeros(ny_global, dtype=np.float64)
        x_seen = np.zeros(nx_global, dtype=bool)
        y_seen = np.zeros(ny_global, dtype=bool)

        for f in rank_files:
            with nc.Dataset(f) as ds:
                i0 = int(ds.getncattr("i_start")) - 1
                j0 = int(ds.getncattr("j_start")) - 1
                nxl = int(ds.getncattr("nx_local"))
                nyl = int(ds.getncattr("ny_local"))
                x_global[i0:i0 + nxl] = ds.variables["x"][:]
                y_global[j0:j0 + nyl] = ds.variables["y"][:]
                x_seen[i0:i0 + nxl] = True
                y_seen[j0:j0 + nyl] = True

        if not (x_seen.all() and y_seen.all()):
            sys.exit("error: rank files do not cover the full global grid "
                     "(gaps in x/y coverage)")

        out.variables["x"][:] = x_global
        out.variables["y"][:] = y_global

        # Time and any other purely non-spatial vars: take from rank 0.
        with nc.Dataset(rank_files[0]) as ref:
            for name, v in ref.variables.items():
                if name in ("x", "y"):
                    continue
                d = v.dimensions
                if "x" in d or "y" in d:
                    continue
                out.variables[name][:] = v[:]

        # ---- spatial fields: paste each rank's slice into global array ----
        # Identify spatial vars and walk them once per rank to keep memory low.
        with nc.Dataset(rank_files[0]) as ref:
            spatial_vars = [
                name for name, v in ref.variables.items()
                if "x" in v.dimensions and "y" in v.dimensions
            ]

        for name in spatial_vars:
            print(f"  merging {name}")
            for f in rank_files:
                with nc.Dataset(f) as ds:
                    i0 = int(ds.getncattr("i_start")) - 1
                    j0 = int(ds.getncattr("j_start")) - 1
                    nxl = int(ds.getncattr("nx_local"))
                    nyl = int(ds.getncattr("ny_local"))
                    src = ds.variables[name]
                    dst = out.variables[name]
                    # Build slicers indexed by dim name (avoids hardcoding order).
                    src_idx = []
                    dst_idx = []
                    for dim in src.dimensions:
                        if dim == "x":
                            src_idx.append(slice(None))
                            dst_idx.append(slice(i0, i0 + nxl))
                        elif dim == "y":
                            src_idx.append(slice(None))
                            dst_idx.append(slice(j0, j0 + nyl))
                        else:
                            src_idx.append(slice(None))
                            dst_idx.append(slice(None))
                    dst[tuple(dst_idx)] = src[tuple(src_idx)]


# --------------------------------------------------------- ocean diagnostics


def merge_ocean_diag(rank_files: List[str], out_path: str) -> None:
    """Merge per-rank ocean diag files into a global file.

    Ocean diag files differ from coastal structured output in three ways:

    1. Ghost cells — spatial arrays are (nx_local + 2*nghost) × (ny_local + 2*nghost);
       the merge strips the ghost band before pasting into the global array.
    2. Per-variable time dims — each variable has its own `time_<name>` unlimited
       dim and a matching 1-D time coord variable; there is no shared `time` dim.
    3. No x/y coordinate variables — `x` and `y` exist only as dimensions.
    """
    with nc.Dataset(rank_files[0]) as ref:
        nx_global = int(ref.getncattr("nx_global"))
        ny_global = int(ref.getncattr("ny_global"))
        nghost = int(ref.getncattr("nghost"))
        ref_attrs = {a: ref.getncattr(a) for a in ref.ncattrs()}
        ref_unlim = {d for d in ref.dimensions if ref.dimensions[d].isunlimited()}
        ref_vars = list(ref.variables.keys())
        # Collect all z dims (z_<name>) and their sizes from the reference file.
        ref_z_dims = {d: len(ref.dimensions[d])
                      for d in ref.dimensions if d.startswith("z_")}

    print(f"[ocean_diag] {len(rank_files)} ranks → {nx_global} x {ny_global} "
          f"(nghost={nghost})")

    skip_attrs = {"i_start", "j_start", "nx_local", "ny_local", "nghost"}

    with nc.Dataset(out_path, "w") as out:
        # Global attrs (drop per-rank window keys).
        for k, v in ref_attrs.items():
            if k in skip_attrs:
                continue
            out.setncattr(k, v)

        # Dimensions: rewrite x, y to global physical sizes; pass through
        # all per-var time dims (unlimited) and z dims unchanged.
        out.createDimension("x", nx_global)
        out.createDimension("y", ny_global)
        with nc.Dataset(rank_files[0]) as ref:
            for d, dim in ref.dimensions.items():
                if d in ("x", "y"):
                    continue
                if dim.isunlimited():
                    out.createDimension(d, None)
                else:
                    out.createDimension(d, len(dim))

        # Define all variables preserving dims and fill/compression settings.
        with nc.Dataset(rank_files[0]) as ref:
            for name in ref_vars:
                make_var(out, ref.variables[name], ref.variables[name].dimensions)

        # Identify spatial vars (those that have both 'x' and 'y' dims) vs
        # non-spatial (1-D time coord vars, which have only a time dim).
        with nc.Dataset(rank_files[0]) as ref:
            spatial_vars = [
                name for name, v in ref.variables.items()
                if "x" in v.dimensions and "y" in v.dimensions
            ]
            non_spatial_vars = [
                name for name in ref_vars if name not in spatial_vars
            ]

        # Non-spatial vars (the per-var 1-D time coordinate arrays): copy
        # from rank 0.  All ranks write the same time values since cadences
        # are globally synchronised.
        with nc.Dataset(rank_files[0]) as ref:
            for name in non_spatial_vars:
                out.variables[name][:] = ref.variables[name][:]

        # ---- spatial fields: strip ghosts, paste each rank's slice ----
        # The per-rank local array has shape (nx_local + 2*ng, ny_local + 2*ng)
        # in x/y.  We read only the physical interior [ng : ng+nx_local] before
        # pasting into the global array at [i_start-1 : i_start-1+nx_local].
        for name in spatial_vars:
            print(f"  merging {name}")
            for f in rank_files:
                with nc.Dataset(f) as ds:
                    i0 = int(ds.getncattr("i_start")) - 1   # 0-based global start
                    j0 = int(ds.getncattr("j_start")) - 1
                    nxl = int(ds.getncattr("nx_local"))
                    nyl = int(ds.getncattr("ny_local"))
                    ng = int(ds.getncattr("nghost"))

                    src = ds.variables[name]
                    dst = out.variables[name]

                    # Build slicers per dim name.
                    src_idx = []
                    dst_idx = []
                    for dim in src.dimensions:
                        if dim == "x":
                            # Strip ghost: physical cells are [ng : ng+nxl]
                            src_idx.append(slice(ng, ng + nxl))
                            dst_idx.append(slice(i0, i0 + nxl))
                        elif dim == "y":
                            src_idx.append(slice(ng, ng + nyl))
                            dst_idx.append(slice(j0, j0 + nyl))
                        else:
                            # Time or z dim: pass through entirely.
                            src_idx.append(slice(None))
                            dst_idx.append(slice(None))

                    dst[tuple(dst_idx)] = src[tuple(src_idx)]


# ------------------------------------------------------------- unstructured


def merge_unstructured(rank_files: List[str], out_path: str,
                       mesh_from: str | None) -> None:
    # Pass 1: scan all ranks, find global cell index range and detect base.
    cgid_min = None
    cgid_max = None
    for f in rank_files:
        with nc.Dataset(f) as ds:
            cgid = ds.variables["cell_global_id"][:]
            cgid_min = cgid.min() if cgid_min is None else min(cgid_min, cgid.min())
            cgid_max = cgid.max() if cgid_max is None else max(cgid_max, cgid.max())

    base = int(cgid_min)
    if base not in (0, 1):
        sys.exit(f"error: cell_global_id base is {base}, expected 0 or 1")
    ncells_global = int(cgid_max) - base + 1

    print(f"[unstructured] {len(rank_files)} ranks → {ncells_global} faces "
          f"(cell_global_id base={base})")

    topo_var_names = {
        "mesh2d", "mesh_node_x", "mesh_node_y", "mesh_face_nodes",
        "cell_global_id",
    }
    # face_x/face_y/face_area are face-centred and merge cleanly via cgid.

    with nc.Dataset(rank_files[0]) as ref, nc.Dataset(out_path, "w") as out:
        ref_attrs = {a: ref.getncattr(a) for a in ref.ncattrs()}
        ref_unlim = {d for d in ref.dimensions if ref.dimensions[d].isunlimited()}

        for k, v in ref_attrs.items():
            if k == "ncells_owned":
                out.setncattr("ncells_owned", ncells_global)
            else:
                out.setncattr(k, v)

        # Dimensions. `face` becomes the global count; node count depends on
        # whether we're copying topology from a source UGRID.
        out.createDimension("face", ncells_global)
        out.createDimension("three", 3)
        for d in ref_unlim:
            out.createDimension(d, None)

        node_count = None
        if mesh_from is not None:
            with nc.Dataset(mesh_from) as src:
                if "node" in src.dimensions:
                    node_count = len(src.dimensions["node"])
                elif "nMesh2_node" in src.dimensions:
                    node_count = len(src.dimensions["nMesh2_node"])
        if node_count is None:
            # No source mesh: keep node dim sized to the rank-0 file as a
            # placeholder; the topology variables are skipped below.
            node_count = (len(ref.dimensions["node"])
                          if "node" in ref.dimensions else 0)
        out.createDimension("node", node_count)

        # Define face-centred + non-topology variables. Topology is either
        # copied verbatim from the UGRID source (--mesh-from) or skipped.
        face_vars: List[str] = []
        for name, v in ref.variables.items():
            if name in topo_var_names:
                continue  # handled below (or skipped)
            make_var(out, v, v.dimensions)
            if "face" in v.dimensions:
                face_vars.append(name)

        # Time + scalar pass-throughs from rank 0.
        for name, v in ref.variables.items():
            if name in topo_var_names or "face" in v.dimensions:
                continue
            out.variables[name][:] = v[:]

        # Mesh topology, if asked for: copy verbatim from source UGRID.
        if mesh_from is not None:
            print(f"  copying mesh topology from {mesh_from}")
            _copy_mesh_topology(mesh_from, out)
        else:
            print("  skipping mesh topology (pass --mesh-from <ugrid.nc> to include)")

        # One-shot coverage check on cell_global_id (independent of variable).
        seen = np.zeros(ncells_global, dtype=bool)
        for f in rank_files:
            with nc.Dataset(f) as ds:
                cgid = ds.variables["cell_global_id"][:].astype(np.int64) - base
                seen[cgid] = True
        missing = int((~seen).sum())
        if missing:
            warnings.warn(
                f"{missing} of {ncells_global} global cells are not owned by "
                f"any rank (mesh / partition gap)"
            )

        # Scatter face-centred variables into a numpy buffer per variable, then
        # write the buffer to netCDF in one slab. Avoids the per-cell fancy-
        # index write penalty on the netCDF variable that dominates wall time.
        for name in face_vars:
            print(f"  merging {name}")
            src0_dims = ref.variables[name].dimensions
            src0_dtype = ref.variables[name].dtype

            if src0_dims == ("face",):
                buf = np.full(ncells_global, np.nan, dtype=src0_dtype) \
                    if np.issubdtype(src0_dtype, np.floating) \
                    else np.zeros(ncells_global, dtype=src0_dtype)
                for f in rank_files:
                    with nc.Dataset(f) as ds:
                        cgid = ds.variables["cell_global_id"][:].astype(np.int64) - base
                        buf[cgid] = ds.variables[name][:]
                out.variables[name][:] = buf

            elif src0_dims and src0_dims[-1] == "face" and src0_dims[0] == "time":
                # (time, face) — assemble the full (nt, ncells_global) slab.
                nt = ref.variables[name].shape[0]
                buf = np.full((nt, ncells_global), np.nan, dtype=src0_dtype) \
                    if np.issubdtype(src0_dtype, np.floating) \
                    else np.zeros((nt, ncells_global), dtype=src0_dtype)
                for f in rank_files:
                    with nc.Dataset(f) as ds:
                        cgid = ds.variables["cell_global_id"][:].astype(np.int64) - base
                        buf[:, cgid] = ds.variables[name][:]
                out.variables[name][:] = buf

            else:
                sys.exit(f"error: unexpected dim layout for face var "
                         f"'{name}': {src0_dims}")


def _copy_mesh_topology(mesh_path: str, out: nc.Dataset) -> None:
    """Copy mesh_node_x, mesh_node_y, mesh_face_nodes (and mesh2d if present)
    from a UGRID source into `out`. Tolerates SCHISM-style alternative names."""
    name_map = {
        "mesh2d":          ["mesh2d", "Mesh2"],
        "mesh_node_x":     ["mesh_node_x", "Mesh2_node_x"],
        "mesh_node_y":     ["mesh_node_y", "Mesh2_node_y"],
        "mesh_face_nodes": ["mesh_face_nodes", "Mesh2_face_nodes"],
    }
    with nc.Dataset(mesh_path) as src:
        for tgt, candidates in name_map.items():
            srcv = next((src.variables[c] for c in candidates if c in src.variables), None)
            if srcv is None:
                if tgt == "mesh2d":
                    continue  # optional
                sys.exit(f"error: {mesh_path} missing {tgt} (tried: {candidates})")
            # Map src dims to our canonical names: any node-sized dim → 'node'
            mapped_dims = []
            for d in srcv.dimensions:
                dlen = len(src.dimensions[d])
                if dlen == len(out.dimensions["node"]):
                    mapped_dims.append("node")
                elif dlen == len(out.dimensions["face"]):
                    mapped_dims.append("face")
                elif dlen == 3:
                    mapped_dims.append("three")
                else:
                    mapped_dims.append(d)
                    if d not in out.dimensions:
                        out.createDimension(d, dlen)
            v = make_var(out, srcv, tuple(mapped_dims))
            v[:] = srcv[:]


# -------------------------------------------------------------- entry point


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_dir", help="directory containing <prefix>_rank_*.nc files")
    p.add_argument("output", help="path for the merged NetCDF file")
    p.add_argument("--prefix", default="output",
                   help="filename prefix of the per-rank files "
                        "(default: 'output', giving 'output_rank_NNNNNN.nc'). "
                        "Use e.g. --prefix seamount_bench_full for ocean diag files.")
    p.add_argument("--grid", choices=("auto", "structured", "unstructured", "ocean_diag"),
                   default="auto", help="grid backend (default: auto-detect)")
    p.add_argument("--mesh-from", default=None,
                   help="(unstructured) UGRID NetCDF whose mesh topology to "
                        "copy into the merged file. If omitted, the merged "
                        "file contains only face-centred fields.")
    args = p.parse_args()

    if os.path.exists(args.output):
        sys.exit(f"error: refuse to overwrite existing {args.output}")

    files = find_rank_files(args.run_dir, prefix=args.prefix)
    grid = args.grid
    if grid == "auto":
        grid = detect_grid(files[0])
        print(f"auto-detected grid: {grid}")

    if grid == "structured":
        if args.mesh_from is not None:
            print("warning: --mesh-from has no effect for structured grids")
        merge_structured(files, args.output)
    elif grid == "ocean_diag":
        if args.mesh_from is not None:
            print("warning: --mesh-from has no effect for ocean_diag files")
        merge_ocean_diag(files, args.output)
    else:
        merge_unstructured(files, args.output, args.mesh_from)

    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
