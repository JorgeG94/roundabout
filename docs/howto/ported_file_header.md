# Header for a file ported from CICE / Icepack

Every roundabout source file that is a port of CICE or Icepack code carries the
BSD-3-Clause header below. The full notice and the pinned upstream versions are
in [`NOTICE`](../../NOTICE); the licence text is in
[`LICENSES/BSD-3-Clause-CICE.txt`](../../LICENSES/BSD-3-Clause-CICE.txt).

## Template

Put it immediately after the module's `!!` summary line (FORD renders it), or as
a plain `!` block above `module`:

```fortran
! SPDX-License-Identifier: BSD-3-Clause
!
! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_mechred.
! Translated to roundabout's do-concurrent/GPU form and modified by the
! roundabout developers; see NOTICE for the full notice.
!
! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
! Redistribution and use in source and binary forms, with or without
! modification, are permitted under the BSD-3-Clause conditions reproduced in
! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
! DISCLAIMED (full disclaimer in NOTICE).
!
! Original authors (from the upstream header): William H. Lipscomb (LANL),
! Elizabeth C. Hunke (LANL).
! Physics: Lipscomb et al. (2007, JGR 112 C03S91); Thorndike et al. (1975);
! Rothrock (1975); Hibler (1980); Flato and Hibler (1995).
```

Replace the `Upstream:` module (and the CICE hash for a CICE-side module), the
authors and the papers per file. The example above is the ridging module.

## Rules

- Name the upstream **module and version**, never upstream file paths or line
  numbers. Everywhere else (code, docs, knobs, tests) cite the papers only.
- Carry over the upstream original-authors line where the source names authors.
  This attributes; it never implies endorsement.
- A wholly new file that implements a published algorithm written from the
  paper, with no side-by-side transliteration of CICE code, gets no BSD block:
  cite the paper only. If code was read side by side and transliterated, it is
  a port and gets the header.
- In `docs/CLOSURE_MATRIX.md` and `docs/CAPABILITIES_AND_LIMITATIONS.md`, each
  ported scheme row says "ported from Icepack (BSD-3)".
- The binary-redistribution condition means `NOTICE` must accompany any release
  tarball and the rendered FORD docs (see `ford.md`).

## Publication acknowledgement

The CICE Consortium Distribution Policy asks that publications using CICE or its
tools acknowledge the Consortium. Suggested text for papers and posters that use
roundabout's CICE-derived ice:

> Sea-ice physics in Roundabout is ported from CICE6/Icepack (CICE Consortium;
> CICE main 2d271f5, Icepack 9a1f3e4; doi:10.5281/zenodo.1205674,
> doi:10.5281/zenodo.1213462). "The CICE Consortium is supported by the U.S.
> Department of Energy, Department of Defense, Department of Commerce, National
> Science Foundation, Environment and Climate Change Canada, Danish
> Meteorological Institute, Australian Community Climate and Earth System
> Simulator National Research Infrastructure, and Polish Academy of Sciences."

The policy also asks users to tell the CICE Consortium about publications and
model use.
