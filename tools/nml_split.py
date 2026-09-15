#!/usr/bin/env python3
"""
Transform every `&ocean_setup_nml ... /` block into per-group sub-namelists
and rewrite `&ocean_diag_nml` keys to drop the `ocean_diag_` prefix.

In-place edits, one file at a time.  Keeps comments + ordering as much as
the simple line-by-line model allows.
"""

import sys
import re
from pathlib import Path

# knob → (group, new_short_name) mapping.
# Source-of-truth for the PR-B2 split.  Keys are the OLD knob names that
# appear in `&ocean_setup_nml` blocks today; values are
# (new sub-nml group, new short key) tuples.
SETUP_REMAP = {
    # coriolis
    "ocean_coriolis_form":               ("coriolis", "form"),
    # thermo
    "ocean_enable_thermodynamics":       ("thermo", "enable_thermodynamics"),
    # tracers
    "ocean_enable_ideal_age":            ("tracers", "enable_ideal_age"),
    # bt
    "ocean_use_bt_cont_type":            ("bt", "use_cont_type"),
    "ocean_bt_cont_corr_bounds":         ("bt", "cont_corr_bounds"),
    "ocean_bt_upstream_h_face":          ("bt", "upstream_h_face"),
    "ocean_bt_correction_h_weighted":    ("bt", "correction_h_weighted"),
    "ocean_bt_correction_visc_rem":      ("bt", "correction_visc_rem"),
    "ocean_bt_correction_bc_pgf":        ("bt", "correction_bc_pgf"),
    "ocean_bt_substep_drag":             ("bt", "substep_drag"),
    "ocean_debug_bt_budget":             ("bt", "debug_budget"),
    "n_inner":                           ("bt", "n_inner"),
    "auto_n_inner":                      ("bt", "auto_n_inner"),
    "cfl_bt_safety":                     ("bt", "cfl_bt_safety"),
    # pgf
    "ocean_pgf_form":                    ("pgf", "form"),
    "ocean_gprime_gfs":                  ("pgf", "gprime_gfs"),
    "ocean_gprime_gint":                 ("pgf", "gprime_gint"),
    "ocean_gfs_scale":                   ("pgf", "gfs_scale"),
    "ocean_maxvel":                      ("pgf", "maxvel"),
    # bdrag
    "ocean_bdrag_form":                  ("bdrag", "form"),
    "ocean_bdrag_cd":                    ("bdrag", "cd"),
    "ocean_bdrag_r":                     ("bdrag", "r"),
    "ocean_bdrag_hbbl":                  ("bdrag", "hbbl"),
    "ocean_bdrag_bg_vel":                ("bdrag", "bg_vel"),
    "ocean_bdrag_bbl_thick_min":         ("bdrag", "bbl_thick_min"),
    "ocean_bdrag_bed_factor":            ("bdrag", "bed_factor"),
    # hvisc
    "ocean_lateral_closure":             ("hvisc", "lateral_closure"),
    "ocean_c_smag":                      ("hvisc", "c_smag"),
    "ocean_c_leith":                     ("hvisc", "c_leith"),
    "ocean_kh_vel_scale":                ("hvisc", "kh_vel_scale"),
    "ocean_ah_bg":                       ("hvisc", "ah_bg"),
    "ocean_ah_max":                      ("hvisc", "ah_max"),
    "ocean_smag_ah":                     ("hvisc", "smag_ah"),
    "ocean_smag_bi_const":               ("hvisc", "smag_bi_const"),
    "ocean_nu_4_bg":                     ("hvisc", "nu_4_bg"),
    "ocean_nu_4_max":                    ("hvisc", "nu_4_max"),
    "ocean_nu_h":                        ("hvisc", "nu_h"),
    "ocean_nu_4":                        ("hvisc", "nu_4"),
    # vmix
    "vmix_use_closure":                  ("vmix", "use_closure"),
    "vmix_use_kpp":                      ("vmix", "use_kpp"),
    "ocean_direct_stress":               ("vmix", "direct_stress"),
    "ocean_hmix_stress":                 ("vmix", "hmix_stress"),
    "ocean_kv_ml_invz2":                 ("vmix", "kv_ml_invz2"),
    "ocean_hmix_fixed":                  ("vmix", "hmix_fixed"),
    "ocean_harmonic_visc":               ("vmix", "harmonic_visc"),
    "ocean_dt_therm_ratio":              ("vmix", "dt_therm_ratio"),
    # continuity
    "ocean_continuity_h_min":            ("continuity", "h_min"),
    "ocean_continuity_ppm_limit_pos":    ("continuity", "ppm_limit_pos"),
    # topo
    "topo_config":                       ("topo", "topo_config"),
    "ocean_max_depth":                   ("topo", "max_depth"),
    "ocean_edge_depth":                  ("topo", "edge_depth"),
    "ocean_slope_scale":                 ("topo", "slope_scale"),
    "wind_config":                       ("topo", "wind_config"),
    "taux_magnitude":                    ("topo", "taux_magnitude"),
    "coriolis_beta":                     ("topo", "coriolis_beta"),
    "coriolis_y_ref":                    ("topo", "coriolis_y_ref"),
    # ic
    "ic_config":                         ("ic", "ic_config"),
    "ocean_alpha_T":                     ("ic", "alpha_T"),
    "ocean_rho_0":                       ("ic", "rho_0"),
    "eady_dT_dy":                        ("ic", "eady_dT_dy"),
    "eady_dT_dz":                        ("ic", "eady_dT_dz"),
    "eady_T_ref":                        ("ic", "eady_T_ref"),
    "eady_pert_amp":                     ("ic", "eady_pert_amp"),
    "eady_pert_seed":                    ("ic", "eady_pert_seed"),
    "ga_eta_amp":                        ("ic", "ga_eta_amp"),
    "ga_length_scale":                   ("ic", "ga_length_scale"),
    "ga_x_center":                       ("ic", "ga_x_center"),
    "ga_y_center":                       ("ic", "ga_y_center"),
}

DIAG_REMAP = {
    "ocean_diag_enabled":    "enabled",
    "ocean_diag_filename":   "filename",
    "ocean_diag_dt_out":     "dt_out",
    "ocean_diag_vgrid":      "vgrid",
    "ocean_diag_z_levels":   "z_levels",
    "ocean_diag_n_z_levels": "n_z_levels",
}

# Fixed group ordering for emitted sub-namelists.
GROUP_ORDER = ["coriolis", "thermo", "tracers", "topo", "ic", "pgf", "bdrag",
               "hvisc", "vmix", "continuity", "bt"]


def parse_assignments(block_body):
    """
    Parse `key = value [! comment]` lines from a namelist block body.
    Returns a list of (key, raw_line, leading_blank_or_comment_lines).
    Preserves blank lines + comment-only lines as `pre` chunks before
    each assignment so the rewrite can keep them attached.
    """
    out = []
    pre = []
    for raw in block_body.splitlines():
        # Match `<key> = <value>` (key is alpha+digit+underscore).
        m = re.match(r'^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=', raw)
        if m:
            key = m.group(2)
            out.append((key, raw, pre))
            pre = []
        else:
            pre.append(raw)
    # Trailing pre = trailing blanks/comments at end of block.
    return out, pre


def split_setup_block(body):
    """
    Take an &ocean_setup_nml block body and split into per-group blocks.
    Returns dict {group_name: [(key, raw_line, pre_lines), ...]} plus a
    list of orphan lines that didn't map (warn the user).
    """
    assignments, trailing = parse_assignments(body)
    per_group = {g: [] for g in GROUP_ORDER}
    orphans = []
    for key, raw, pre in assignments:
        if key in SETUP_REMAP:
            group, new_key = SETUP_REMAP[key]
            # Rewrite the key in the raw line.
            new_raw = re.sub(
                r'^(\s*)' + re.escape(key) + r'(\s*=)',
                lambda m: m.group(1) + new_key + m.group(2),
                raw,
                count=1,
            )
            per_group[group].append((new_key, new_raw, pre))
        else:
            orphans.append((key, raw))
    return per_group, orphans, trailing


def rewrite_diag_block(body):
    """Rewrite &ocean_diag_nml keys to drop the `ocean_diag_` prefix."""
    lines = []
    for raw in body.splitlines():
        new = raw
        for old_key, new_key in DIAG_REMAP.items():
            new = re.sub(
                r'^(\s*)' + re.escape(old_key) + r'(\s*=)',
                lambda m, nk=new_key: m.group(1) + nk + m.group(2),
                new,
                count=1,
            )
        lines.append(new)
    return '\n'.join(lines)


def render_setup_replacement(per_group, trailing, indent="   "):
    """Render the per-group sub-namelists into one text block."""
    out_blocks = []
    for group in GROUP_ORDER:
        items = per_group[group]
        if not items:
            continue
        chunk = [f"&ocean_{group}_nml"]
        for key, raw, pre in items:
            for p in pre:
                if p.strip():    # skip leading blanks at start of group
                    chunk.append(p)
            chunk.append(raw)
        chunk.append("/")
        out_blocks.append('\n'.join(chunk))
    if trailing:
        # Drop trailing blanks; comments at the very end get lost.  Most
        # of our nmls don't have those.
        pass
    return '\n\n'.join(out_blocks)


def process_file(path):
    text = path.read_text()
    orig = text
    notes = []

    # ---- &ocean_setup_nml → per-group sub-nmls ----
    setup_re = re.compile(
        r'(&ocean_setup_nml\b[ \t]*\n)(.*?)(^\s*/\s*$)',
        re.DOTALL | re.MULTILINE,
    )
    m = setup_re.search(text)
    if m:
        body = m.group(2)
        per_group, orphans, trailing = split_setup_block(body)
        if orphans:
            notes.append(f"  ⚠ {path.name} orphan keys: {[k for k, _ in orphans]}")
        replacement = render_setup_replacement(per_group, trailing)
        text = text[:m.start()] + replacement + text[m.end():]

    # ---- &ocean_diag_nml: rewrite keys, keep block intact ----
    diag_re = re.compile(
        r'(&ocean_diag_nml\b[ \t]*\n)(.*?)(^\s*/\s*$)',
        re.DOTALL | re.MULTILINE,
    )
    m = diag_re.search(text)
    if m:
        new_body = rewrite_diag_block(m.group(2))
        text = text[:m.start()] + m.group(1) + new_body + '\n' + m.group(3) + text[m.end():]

    if text != orig:
        path.write_text(text)
        notes.insert(0, f"  ✓ updated {path}")
    else:
        notes.insert(0, f"  (no ocean nmls in {path})")
    return notes


def main():
    if len(sys.argv) < 2:
        print("usage: nml_split.py <file.nml> [<file.nml> ...]", file=sys.stderr)
        sys.exit(1)
    for p in sys.argv[1:]:
        for note in process_file(Path(p)):
            print(note)


if __name__ == '__main__':
    main()
