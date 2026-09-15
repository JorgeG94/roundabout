#!/usr/bin/env python3
"""Rewrite `do concurrent` blocks as `!$omp target teams distribute parallel do`.

    do concurrent (j=1:ny, i=1:nx) local(a, b)
        ...
    end do

becomes

    !$omp target teams distribute parallel do collapse(2) private(a, b)
    do j = 1, ny
    do i = 1, nx
        ...
    end do
    end do
    !$omp end target teams distribute parallel do

Rules:
  - Preserve header order: leftmost induction → outermost loop
  - `local(...)` → `private(...)`; a specifier wrapped in a portability
    macro (`DO_LOCALITY(local(a, b))`, MOM6-style) is unwrapped first
  - Drop `collapse(N)` when N == 1
  - Match the `end do` of the block by depth-tracking forward
  - Skip occurrences inside comments / string literals
  - Join `&` continuations on the header
  - Strip the `pure` attribute from any procedure that ends up containing an
    `!$omp target` region.  An executable target region launches a device
    kernel + moves data — side effects a `pure` procedure may not have, so
    NVHPC/ifx reject `pure` + `!$omp target`.  (`do concurrent` itself is
    side-effect-free and is fine in a pure procedure, which is why the
    OpenACC/`-stdpar` source keeps `pure`; only this OpenMP-target variant
    needs it removed.)  Declarative directives like `!$omp declare target`
    are pure-compatible and are left alone.

Targets (--target):
  gpu (default)  `!$omp target teams distribute parallel do` (offload).
  cpu            `!$omp parallel do` (host multicore) — drops `target teams
                 distribute`, whose host-fallback path is flaky on some
                 compilers (notably GNU). The `pure`-strip cascade keys on the
                 host `!$omp parallel do` opening directive instead of
                 `!$omp target` (a host parallel region is likewise rejected
                 inside a PURE procedure). Build with
                 -DRDB_PARALLEL_BACKEND=openmp -DRDB_ENABLE_GPU=OFF.

Usage:
    python tools/dc_to_omp.py                  # dry-run, scans src/ (gpu)
    python tools/dc_to_omp.py --write          # apply in place (gpu)
    python tools/dc_to_omp.py --target cpu --write src   # host multicore
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from functools import partial
from pathlib import Path

from _locality_macro import unwrap_locality_macros
from _parallel import pmap, resolve_jobs
from _progress import Progress, track

DC_RE = re.compile(r"\bdo\s+concurrent\b", re.IGNORECASE)
DO_RE = re.compile(r"^\s*(?:[a-zA-Z_][a-zA-Z0-9_]*\s*:\s*)?do\b", re.IGNORECASE)
END_DO_RE = re.compile(r"^\s*end\s*do\b", re.IGNORECASE)
KNOWN_CLAUSES = ("local_init", "local", "shared", "default", "reduce")

# An *opening* `!$omp target` directive (not `!$omp end target`, and not a
# declarative `!$omp declare target`).
OMP_TARGET_RE = re.compile(r"^\s*!\$omp\s+target\b", re.IGNORECASE)
# The CPU-target opening worksharing directive (`!$omp parallel do`, not its
# `!$omp end parallel do`). Like a target region, a host parallel region is a
# side-effecting construct that the compiler rejects inside a PURE procedure.
OMP_PARALLEL_DO_RE = re.compile(r"^\s*!\$omp\s+parallel\s+do\b", re.IGNORECASE)
# A procedure-opening statement (subroutine/function with a name), e.g.
# `pure subroutine foo(`, `pure real(wp) function bar(`.
PROC_OPEN_RE = re.compile(
    r"^\s*[^!]*\b(?:subroutine|function)\b\s+[a-zA-Z]\w*", re.IGNORECASE
)
PROC_END_RE = re.compile(r"^\s*end\s*(?:subroutine|function)\b", re.IGNORECASE)
# The `pure` prefix token plus trailing blanks (leaves indentation intact;
# `\bpure\b` does not match inside `impure`).
PURE_RE = re.compile(r"\bpure\b[ \t]*", re.IGNORECASE)


def strip_strings(line: str) -> str:
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(" " if ch != quote else ch)
            if ch == quote:
                quote = None
        else:
            if ch in ("'", '"'):
                quote = ch
            out.append(ch)
    return "".join(out)


def code_part(line: str) -> str:
    s = strip_strings(line)
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def split_top_level(s: str, sep: str = ",") -> list[str]:
    out, buf, depth = [], [], 0
    for ch in s:
        if ch == "(":
            depth += 1; buf.append(ch)
        elif ch == ")":
            depth -= 1; buf.append(ch)
        elif ch == sep and depth == 0:
            out.append("".join(buf).strip()); buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf).strip())
    return out


@dataclass
class Header:
    indent: str
    indices: list[tuple[str, str, str, str | None]]  # (var, lo, hi, step)
    privates: list[str]
    has_mask: bool
    mask: str | None = None  #: lowered to `if (mask) then` inside the loop nest


def parse_triplet(item: str) -> tuple[str, str, str, str | None] | None:
    if "=" not in item:
        return None
    name, _, rng = item.partition("=")
    name = name.strip()
    parts = split_top_level(rng, ":")
    if len(parts) < 2:
        return None
    lo = parts[0].strip()
    hi = parts[1].strip()
    step = parts[2].strip() if len(parts) >= 3 else None
    return (name, lo, hi, step)


COMPOUND_END_DO_RE = re.compile(r"^\s*end\s*do\b", re.IGNORECASE)


def _split_positions(s: str) -> list[int]:
    """Indices of top-level `;` in `s` (paren depth 0)."""
    pos, depth = [], 0
    for k, ch in enumerate(s):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == ";" and depth == 0:
            pos.append(k)
    return pos


def normalize_compound_lines(lines: list[str]) -> tuple[list[str], list[int]]:
    """Split `;`-compound lines that the block machinery needs to own outright.

    Fortran allows several statements per line. MOM6 uses that everywhere:

        if (use_temperature) then ; do concurrent (j=js:je, i=is:ie)
        do k=1,nz ; do concurrent (i=is:ie)
        do concurrent (i=is:ie) ; press(i,j) = 0.0 ; enddo
        endif ; enddo

    The emitter replaces WHOLE LINES — a directive plus one `do` per index in
    place of the header, and the closers in place of the `end do`. Any
    statement sharing one of those lines would be destroyed by the rewrite. So
    rather than refuse such loops (MOM6 has ~180), give each statement its own
    line first and let the line-oriented machinery downstream be exactly right.

    Only two kinds of line are touched, so the diff stays proportional to the
    loops actually converted instead of reflowing whole files:
      * a line whose code holds a `do concurrent` plus anything else, and
      * a multi-statement line containing an `end do` / `enddo`.

    Never touched: comment and directive lines, preprocessor lines, and any
    line involved in a `&` continuation — statement boundaries there are not
    what a `;` split would suggest.

    Returns the new lines plus an `origin` map (new index -> original 0-based
    line number) so diagnostics can still name the line the reader has.
    """
    out: list[str] = []
    origin: list[int] = []
    in_continuation = False
    for n, raw in enumerate(lines):
        stripped = raw.lstrip()
        code = code_part(raw)
        bare = strip_strings(code)
        was_cont, continues = in_continuation, bare.rstrip().endswith("&")
        in_continuation = continues

        # A line that ENDS in `&` is still safe to split: every `;`-separated
        # piece before the last is a complete statement, and the `&` stays
        # attached to the final piece. Only a line that IS a continuation of a
        # previous one is off limits — there a `;` may sit inside a continued
        # expression, where it is not a statement boundary.
        # `do K=2,nz ; do concurrent(I=is-1:ie) &` is exactly this shape, and
        # skipping it silently deleted the `do K=2,nz` when the header line was
        # replaced.
        cuts = [] if was_cont else _split_positions(bare)
        if (not cuts or not stripped or stripped[0] in "!#"):
            out.append(raw); origin.append(n); continue

        # Slice the ORIGINAL text at those positions (string literals intact).
        bounds = [-1] + cuts + [len(code)]
        pieces = [code[a + 1:b].strip() for a, b in zip(bounds, bounds[1:])]
        pieces = [x for x in pieces if x]
        if len(pieces) < 2:
            out.append(raw); origin.append(n); continue

        has_dc = any(DC_RE.search(strip_strings(x)) for x in pieces)
        has_end_do = any(END_DO_RE.match(x) for x in pieces)
        if not (has_dc or has_end_do):
            out.append(raw); origin.append(n); continue

        bang = strip_strings(raw).find("!")
        comment = raw[bang:].rstrip() if bang >= 0 else ""
        indent = raw[: len(raw) - len(stripped)]
        for k, piece in enumerate(pieces):
            tail = f"  {comment}" if (comment and k == len(pieces) - 1) else ""
            out.append(f"{indent}{piece}{tail}")
            origin.append(n)
    return out, origin


# A hand-written OpenMP COMPUTE region: `!$omp target teams ...` (loop,
# distribute parallel do, ...). Data-mapping directives (`target data`,
# `target enter/exit data`, `target update`) are NOT compute regions — a
# compute construct nested inside those is perfectly legal.
HANDWRITTEN_TARGET_RE = re.compile(r"^!\$omp\s+target\s+teams\b", re.IGNORECASE)


def handwritten_target_extent(lines: list[str]) -> set[int]:
    """Indices covered by an existing `!$omp target teams` region's loop nest.

    A `do concurrent` inside one of these must NOT be given a directive of its
    own: OpenMP forbids a target region inside a target region, and the loop is
    already covered by the enclosing construct. MOM6's GPU branch has ~50 such
    loops — a `!$omp target teams loop` over `j`, with `do concurrent (i=...)`
    inside — so this is the common case, not a corner.

    The construct's `end` directive is OPTIONAL for a loop-associated one, so
    the extent is found by walking the associated loop nest to its close rather
    than by looking for `!$omp end`.
    """
    covered: set[int] = set()
    n = 0
    while n < len(lines):
        s = code_part(lines[n]).strip()
        if not HANDWRITTEN_TARGET_RE.match(lines[n].strip()):
            n += 1
            continue
        m = n
        while m < len(lines) and lines[m].rstrip().endswith("&"):
            m += 1
        depth, k, started = 0, m + 1, False
        while k < len(lines):
            for st in statements(code_part(lines[k])):
                if not st:
                    continue
                if END_DO_RE.match(st):
                    depth -= 1
                elif DO_RE.match(st):
                    depth += 1
                    started = True
            covered.add(k)
            if started and depth <= 0:
                break
            k += 1
        n = m + 1
    return covered


def _src_line(idx: int, line_delta: int, origin: list[int] | None) -> int:
    """Map a buffer index back to a 1-based ORIGINAL source line number.

    Two shifts stand between the two: `line_delta` lines added by conversions
    already applied to this buffer, and the compound-line normalisation that
    ran before any of them (recorded in `origin`). A skip message is the only
    handle on what needs a manual patch, so one naming a line the loop is not
    on is worse than none.
    """
    j = idx - line_delta
    if origin is not None and 0 <= j < len(origin):
        return origin[j] + 1
    return j + 1


def find_dc_block(lines: list[str], start: int, line_delta: int = 0,
                  label: str = "",
                  origin: list[int] | None = None) -> tuple[Header, int, int] | None:
    """Locate the next CONVERTIBLE `do concurrent` block at or after `start`.

    Returns (header, header_first_line, header_last_line), or None when the
    file holds no further convertible block.

    A header we cannot convert (mask, unbalanced parens, no parseable index
    triplet) is SKIPPED and the scan continues past it — it must not return
    None, because the caller reads None as "end of file" and stops. Conflating
    the two once meant a single masked loop silenced every remaining loop in
    the file: on MOM6, 57 masked headers produced 5 skip messages and
    abandoned several hundred perfectly convertible loops downstream.
    """
    i = start
    while i < len(lines):
        c = code_part(lines[i])
        # Detect on the STRING-STRIPPED text: a character literal that
        # mentions the construct (`call check(..., "do concurrent reduce(+) ..."`,
        # a test assertion message) is not a loop. Matching it derails the whole
        # scan of that file.
        m = DC_RE.search(strip_strings(c))
        if not m:
            i += 1
            continue
        # Join continuations on the header.
        joined_parts = []
        first = i
        j = i
        while j < len(lines):
            piece = code_part(lines[j]).rstrip()
            if piece.endswith("&"):
                joined_parts.append(piece[:-1])
                j += 1
                if j < len(lines):
                    nxt = lines[j].lstrip()
                    if nxt.startswith("&"):
                        lines[j] = lines[j].replace("&", " ", 1)
            else:
                joined_parts.append(piece)
                break
        joined = " ".join(joined_parts)
        # Unwrap CPP-hidden locality specifiers (MOM6's `DO_LOCALITY(local(...))`)
        # BEFORE clause parsing. We read unpreprocessed source, so a wrapped
        # `local(...)` would otherwise be invisible and the emitted OpenMP loop
        # would silently drop its `private(...)` — a wrong-answer bug.
        joined = unwrap_locality_macros(joined)
        last = j

        # Parse header.
        m2 = DC_RE.search(joined)
        rest = joined[m2.end():].lstrip()
        if not rest.startswith("("):
            sys.stderr.write(
                f"  {label}skip: unparseable do-concurrent header at line {_src_line(first, line_delta, origin)}\n")
            i = last + 1
            continue
        depth, end = 0, -1
        for k, ch in enumerate(rest):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    end = k
                    break
        if end < 0:
            sys.stderr.write(
                f"  {label}skip: unbalanced parens in header at line {_src_line(first, line_delta, origin)}\n")
            i = last + 1
            continue
        inside = rest[1:end]
        trailing = rest[end + 1:].strip()

        indices: list[tuple[str, str, str, str | None]] = []
        masks: list[str] = []
        for it in split_top_level(inside, ","):
            t = parse_triplet(it)
            if t is not None:
                indices.append(t)
            elif it.strip():
                masks.append(it.strip())
        has_mask = bool(masks)

        privates: list[str] = []
        t = trailing
        while t:
            tlow = t.lower()
            matched = None
            for c in KNOWN_CLAUSES:
                if tlow.startswith(c) and (len(t) == len(c) or t[len(c)] in "( "):
                    matched = c
                    break
            if not matched:
                break
            after = t[len(matched):].lstrip()
            if after.startswith("("):
                d, jj = 0, -1
                for k, ch in enumerate(after):
                    if ch == "(":
                        d += 1
                    elif ch == ")":
                        d -= 1
                        if d == 0:
                            jj = k
                            break
                if jj < 0:
                    break
                payload = after[1:jj]
                if matched == "local":
                    privates.extend(v.strip() for v in split_top_level(payload, ","))
                t = after[jj + 1:].lstrip()
            else:
                t = after

        # Code BEFORE the construct on the header line would be destroyed by
        # the replacement just as trailing code would. The normaliser splits
        # such lines, so reaching here means it declined to (a continuation it
        # could not safely cut) — refuse rather than delete the statement.
        prefix = code_part(lines[first])[:
            DC_RE.search(strip_strings(code_part(lines[first]))).start()]
        if prefix.strip():
            sys.stderr.write(
                f"  {label}skip: header at line {_src_line(first, line_delta, origin)} "
                f"carries leading code ({prefix.strip()!r}); manual handling needed\n")
            i = last + 1
            continue

        # Anything left in `t` is real code the header carries after its
        # clauses — MOM6's `do concurrent (j=js:je, I=is-1:ie) ; if (...) then`.
        # The emitter replaces the whole header line, so converting would
        # DELETE that statement. Refuse.
        if t.strip():
            sys.stderr.write(
                f"  {label}skip: header at line {_src_line(first, line_delta, origin)} carries "
                f"trailing statements ({t.strip()!r}); manual handling needed\n")
            i = last + 1
            continue

        if not indices:
            sys.stderr.write(
                f"  {label}skip: no parseable index triplet at line {_src_line(first, line_delta, origin)}\n")
            i = last + 1
            continue

        indent = lines[first][: len(lines[first]) - len(lines[first].lstrip())]
        return (Header(indent=indent, indices=indices, privates=privates,
                       has_mask=has_mask,
                       mask=" .and. ".join(f"({m})" for m in masks) if masks else None),
                first, last)
    return None


def statements(code: str) -> list[str]:
    r"""Split one line's code into its `;`-separated Fortran statements.

    Fortran allows several statements per line separated by `;`, and MOM6 uses
    that constantly (`endif ; enddo`, `if (m==1) then ; do concurrent (...) ;
    press(i,j) = 0.0 ; enddo`). A line-ANCHORED `^\s*do` / `^\s*end\s*do`
    matcher miscounts those lines in both directions: the `enddo` in
    `endif ; enddo` is invisible so depth never closes, and the `do` in
    `... then ; do ...` is invisible so depth never opens.

    The second is the dangerous one. An unopened `do` whose `enddo` IS at line
    start makes the matcher decrement past the loop it was tracking and return
    an `end do` belonging to an OUTER construct — structurally broken output,
    emitted silently. Splitting into statements first makes both exact.

    Split on the string-stripped text (length-preserving) so a `;` inside a
    character literal never splits a statement.
    """
    return [s.strip() for s in split_top_level(strip_strings(code), ";")]


def find_matching_end_do(lines: list[str], header_last: int) -> int | None:
    """Scan forward from the line after the header for the matching `end do`.

    Tracks `do` depth statement-by-statement (any `do` opens, `end do`
    closes). Returns None when there is no match, or when the closing
    `end do` shares its line with other statements: the caller replaces that
    whole line with the loop footer, so anything else on it would be lost.
    """
    depth = 1
    i = header_last + 1
    while i < len(lines):
        sts = statements(code_part(lines[i]))
        for n, st in enumerate(sts):
            if END_DO_RE.match(st):
                depth -= 1
                if depth == 0:
                    # Must be the only statement on the line, or the emitted
                    # footer would swallow its neighbours.
                    if any(s for k, s in enumerate(sts) if k != n):
                        return None
                    return i
            elif DO_RE.match(st):
                depth += 1
        i += 1
    return None


# Worksharing directive emitted per target. CPU drops `target teams distribute`
# (the host-fallback path that is flaky on some compilers, notably GNU) for a
# plain host `parallel do`.
DIRECTIVE_KW = {
    "gpu": "target teams distribute parallel do",
    "cpu": "parallel do",
}


def emit_replacement(h: Header, target: str = "gpu", nested: bool = False,
                     extra_privates: list[str] | None = None
                     ) -> tuple[list[str], list[str]]:
    """Return (header_lines, footer_lines).

    `nested` emits the loop nest with NO directive. A `do concurrent` sitting
    inside another one is already covered by the outer construct's
    parallelism, and OpenMP forbids a `target` region inside a `target`
    region, so the inner loop becomes a plain sequential nest. Its `local(...)`
    variables are hoisted into the OUTER directive's `private(...)` by the
    caller — the outer iterations are the parallel unit, so per-thread copies
    there give the inner loop exactly the privacy `local` asked for.
    """
    n = len(h.indices)
    kw = DIRECTIVE_KW[target]
    privates = list(h.privates) + list(extra_privates or [])
    seen: set[str] = set()
    privates = [x for x in privates
                if not (x.lower() in seen or seen.add(x.lower()))]
    parts = [f"!$omp {kw}"]
    if n > 1:
        parts.append(f"collapse({n})")
    if privates:
        parts.append(f"private({', '.join(privates)})")
    directive = h.indent + " ".join(parts)
    do_lines = []
    for (var, lo, hi, step) in h.indices:
        rng = f"{lo}, {hi}" + (f", {step}" if step else "")
        do_lines.append(f"{h.indent}do {var} = {rng}")
    end_lines = [f"{h.indent}end do" for _ in h.indices]
    if not nested:
        end_lines.append(f"{h.indent}!$omp end {kw}")
    if h.mask:
        # `do concurrent (i=lo:hi, MASK)` runs the body only where MASK holds.
        # OpenMP has no mask clause, so it becomes a guard inside the nest.
        #
        # Equivalent for any LEGAL do concurrent: the standard requires each
        # iteration to be independent, so the body cannot write what the mask
        # reads. (Were it to, the mask is conceptually evaluated for all indices
        # before the first iteration, and a guard evaluated per-iteration would
        # differ — but such a loop is already invalid as `do concurrent`.)
        do_lines.append(f"{h.indent}if ({h.mask}) then")
        end_lines.insert(0, f"{h.indent}end if")
    return ((do_lines if nested else [directive] + do_lines), end_lines)


def transform_lines(lines: list[str], target: str = "gpu",
                    label: str = "", nested: bool = False
                    ) -> tuple[list[str], int] | tuple[list[str], int, list[str]]:
    """Transform every `do concurrent` block. Returns (new_lines, n_changed).

    Nested loops are handled by recursing into each block's body BEFORE
    emitting its directive: an inner `do concurrent` becomes a plain `do` nest
    (no directive — OpenMP forbids a target region inside a target region) and
    surrenders its `local(...)` variables to the outer directive's
    `private(...)`. With `nested=True` the call itself emits no directives and
    returns those hoisted names as a third element.
    """
    # Give every header and every `end do` a line of its own before scanning:
    # the emitter rewrites whole lines, so a statement sharing one would be lost.
    if any(DC_RE.search(strip_strings(code_part(l))) for l in lines):
        out, origin = normalize_compound_lines(lines)
    else:
        out, origin = list(lines), list(range(len(lines)))
    # Loops already inside a hand-written `!$omp target teams` region must not
    # get a directive of their own. Computed once, on the normalised buffer,
    # before any rewriting shifts the indices.
    covered = handwritten_target_extent(out) if not nested else set()
    hoisted: list[str] = []
    changed = 0
    i = 0
    # Lines added to `out` so far. Reported line numbers subtract it so every
    # warning points at the original source, not the partially rewritten buffer.
    delta = 0
    while i < len(out):
        block = find_dc_block(out, i, delta, label, origin)
        if block is None:
            break
        h, first, last = block
        end_idx = find_matching_end_do(out, last)
        if end_idx is None:
            sys.stderr.write(f"  {label}warn: no matching end do for header "
                             f"at line {_src_line(first, delta, origin)}\n")
            i = last + 1
            continue
        # Convert the body first: nested loops lose their directives and give
        # up their locality names, which this header must then declare.
        body, body_changed, body_privates = transform_lines(
            out[last + 1:end_idx], target, label, nested=True)
        in_handwritten = (first - delta) in covered
        if in_handwritten and (h.privates or body_privates):
            # Nothing to hoist them into — editing MOM6's own directive is not
            # this tool's business. Say so instead of dropping them silently.
            names = ", ".join(h.privates + body_privates)
            sys.stderr.write(
                f"  {label}warn: loop at line {_src_line(first, delta, origin)} is inside a "
                f"hand-written !$omp target region and has locality ({names}); "
                f"add these to that directive's private() by hand\n")
        new_header, new_footer = emit_replacement(
            h, target, nested=nested or in_handwritten,
            extra_privates=body_privates)
        if nested:
            # No directive here either; pass every private further out.
            hoisted.extend(h.privates)
            hoisted.extend(body_privates)
        out = (out[:first] + new_header + body + new_footer + out[end_idx + 1:])
        changed += 1 + body_changed
        new_len = len(new_header) + len(body) + len(new_footer)
        delta += new_len - (end_idx - first + 1)
        i = first + new_len
    if nested:
        return (out, changed, hoisted)
    return (out, changed)


PROC_NAME_RE = re.compile(r"\b(?:subroutine|function)\b\s+([a-zA-Z]\w*)",
                          re.IGNORECASE)


def _parse_procedures(lines: list[str]) -> list[dict]:
    """List every subroutine/function as {open, end, name}.

    Stack-matches procedure-opening statements to their `end subroutine` /
    `end function` so nested (contained) procedures pair correctly.
    """
    procs: list[dict] = []
    stack: list[tuple[int, str]] = []
    for idx, raw in enumerate(lines):
        c = code_part(raw)
        if PROC_END_RE.match(c):
            if stack:
                open_idx, name = stack.pop()
                procs.append({"open": open_idx, "end": idx, "name": name})
        elif PROC_OPEN_RE.match(c):
            m = PROC_NAME_RE.search(c)
            stack.append((idx, m.group(1).lower() if m else ""))
    return procs


def seed_target_strips(lines: list[str], procs: list[dict],
                       open_re: re.Pattern = OMP_TARGET_RE) -> set[int]:
    """Open-line indices of procedures that DIRECTLY contain an opening
    worksharing directive (`!$omp target` on gpu, `!$omp parallel do` on cpu) —
    those launch a kernel / spawn a team (a side effect) and so can't stay
    `pure`.  The nearest enclosing procedure of each such directive is the one
    that loses purity."""
    seed: set[int] = set()
    for idx, raw in enumerate(lines):
        if not open_re.match(raw):
            continue
        enclosing = None
        for p in procs:
            if p["open"] < idx < p["end"]:
                if enclosing is None or (p["end"] - p["open"]) < (
                    enclosing["end"] - enclosing["open"]):
                    enclosing = p
        if enclosing is not None:
            seed.add(enclosing["open"])
    return seed


IDENT_RE = re.compile(r"[A-Za-z_]\w*")


def _convert_task(path: Path, target: str = "gpu") -> tuple:
    """Phase-1 unit of work: convert one file and seed its local impure names.

    Self-contained (reads its own file, shares nothing), so it is safe to run in
    a worker process. Returns the `unit` tuple plus the names this file
    contributes to the global impure set.
    """
    open_re = OMP_PARALLEL_DO_RE if target == "cpu" else OMP_TARGET_RE
    try:
        text = path.read_text()
    except (UnicodeDecodeError, OSError):
        # Dangling symlink / unreadable file: skip it rather than abort the run.
        return (path, [], 0, [], set(), True, set())
    lines, n_loops = transform_lines(text.splitlines(), target, f"{path}: ")
    procs = _parse_procedures(lines)
    strip_idx = seed_target_strips(lines, procs, open_re)
    # Seed with every procedure that is non-pure after conversion: the
    # target-stripped ones AND any procedure that already lacks `pure`
    # (legitimately impure, or stripped by a prior run — re-running must still
    # cascade to ITS pure callers). In valid Fortran a pure procedure never
    # calls an already-impure one, so the already-impure names never trigger a
    # spurious strip.
    local_impure = {p["name"] for p in procs
                    if p["open"] in strip_idx
                    or not PURE_RE.search(lines[p["open"]])}
    return (path, lines, n_loops, procs, strip_idx,
            text.endswith("\n"), local_impure)


def _write_task(unit: tuple, write: bool = False) -> tuple:
    """Phase-3 unit of work: apply this file's pure-strips and write it back."""
    path, lines, n_loops, _procs, strip_idx, trailing_nl = unit
    n_pure = 0
    for idx in strip_idx:
        if PURE_RE.search(lines[idx]):
            lines[idx] = PURE_RE.sub("", lines[idx], count=1)
            n_pure += 1
    if (n_loops or n_pure) and write:
        path.write_text("\n".join(lines) + ("\n" if trailing_nl else ""))
    return path, n_loops, n_pure


def walk(roots: list[Path]) -> list[Path]:
    files: list[Path] = []
    for r in roots:
        if r.is_file():
            files.append(r); continue
        for p in sorted(r.rglob("*.F90")):
            files.append(p)
        for p in sorted(r.rglob("*.f90")):
            files.append(p)
    return files


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="apply changes in place (default: dry-run)")
    ap.add_argument(
        "--target", choices=["gpu", "cpu"], default="gpu",
        help="gpu (default): rewrite do concurrent as !$omp target teams "
             "distribute parallel do. cpu: !$omp parallel do (host multicore).")
    ap.add_argument(
        "-j", "--jobs", "-t", dest="jobs", type=int, default=1, metavar="N",
        help="run the per-file phases across N worker processes (0 = one per "
             "available core). The whole-program pure cascade stays serial; "
             "output is order-independent either way.")
    ap.add_argument("paths", nargs="*", default=["src"])
    args = ap.parse_args()

    roots = [Path(p) for p in args.paths]
    for r in roots:
        if not r.exists():
            print(f"error: {r} does not exist", file=sys.stderr)
            return 1
    jobs = resolve_jobs(args.jobs)

    # ---- Phase 1: per-file DC→target conversion + seed the impure set ----
    # `procs` are parsed from the POST-conversion lines so body spans line up.
    # `impure` is the GLOBAL set of procedure names that lose `pure` — seeded
    # with every procedure that directly contains an `!$omp target`.
    units = []  # [path, lines, n_loops, procs, strip_idx, trailing_nl]
    impure: set[str] = set()
    for (path, lines, n_loops, procs, strip_idx, trailing_nl,
         local_impure) in pmap(partial(_convert_task, target=args.target),
                               walk(roots), jobs, "[1/3] convert"):
        impure |= local_impure
        units.append([path, lines, n_loops, procs, strip_idx, trailing_nl])

    # ---- Phase 2: WHOLE-PROGRAM cascade to a fixpoint ----
    # A `pure` procedure that references an impure procedure (in ANY file —
    # a pure procedure may only call pure procedures) must lose `pure` too.
    # Iterate across all files until the impure set stops growing, so the
    # impurity propagates up cross-module call chains (the impure callee
    # lives in one file, its caller in another module in a different file).
    #
    # Run as a worklist over a reverse index rather than re-sweeping every
    # file. The old form rebuilt each procedure body and ran one regex per
    # impure name over it, on every sweep — O(sweeps x procs x |impure|), and
    # |impure| grows into the thousands, which is what made this the slow
    # stage. A `\b`-anchored case-insensitive search for a Fortran identifier
    # hits iff that identifier occurs as a token in the body, so a lowercase
    # token set is an exact substitute for the regex; the fixpoint is
    # unchanged, it is just reached by propagating along the call graph
    # instead of rediscovering it each pass.
    candidates = []            # [unit_idx, open_idx, name_lc] still strippable
    referenced_by: dict[str, list[int]] = {}
    for u_i, (_p, lines, _nl, procs, strip_idx, _t) in enumerate(
            track(units, "[2/3] index")):
        for p in procs:
            if p["open"] in strip_idx:
                continue
            if not PURE_RE.search(lines[p["open"]]):
                continue  # already non-pure / nothing to strip
            body = "\n".join(code_part(lines[k])
                             for k in range(p["open"] + 1, p["end"]))
            c_i = len(candidates)
            candidates.append([u_i, p["open"], p["name"].lower()])
            for tok in {m.group(0).lower() for m in IDENT_RE.finditer(body)}:
                referenced_by.setdefault(tok, []).append(c_i)

    cascade = Progress("[2/3] pure cascade", total=None)
    cascade.draw(force=True)
    stripped = [False] * len(candidates)
    impure_lc = {c.lower() for c in impure if c}
    work = list(impure_lc)
    while work:
        callee = work.pop()
        for c_i in referenced_by.get(callee, ()):
            if stripped[c_i]:
                continue
            u_i, open_idx, name_lc = candidates[c_i]
            stripped[c_i] = True
            units[u_i][4].add(open_idx)
            cascade.advance(suffix=f"{len(impure_lc)} impure")
            if name_lc and name_lc not in impure_lc:
                impure_lc.add(name_lc)
                work.append(name_lc)
    cascade.done(f"{len(impure_lc)} impure procedures")

    # ---- Phase 3: apply pure-strips + write ----
    total = 0
    total_pure = 0
    files_changed = 0
    for path, n_loops, n_pure in pmap(partial(_write_task, write=args.write),
                                      units, jobs, "[3/3] write"):
        if n_loops == 0 and n_pure == 0:
            continue
        total += n_loops
        total_pure += n_pure
        files_changed += 1
        verb = "rewrote" if args.write else "would rewrite"
        extra = f"  (-{n_pure} pure)" if n_pure else ""
        print(f"  {verb} {n_loops:3d}  {path}{extra}")
    mode = "WRITE" if args.write else "DRY RUN"
    print(f"\n[{mode}] {total} loops, {total_pure} pure attrs removed "
          f"across {files_changed} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
