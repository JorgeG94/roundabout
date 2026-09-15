"""Unwrap CPP locality-specifier macros in `do concurrent` headers.

Not every Fortran compiler accepts F2018 locality specifiers (`local`,
`local_init`, `shared`, `default(none)`) on a `do concurrent` header, so a
portable codebase hides them behind a CPP macro that expands to the specifier
where it is supported and to nothing where it is not. MOM6 does exactly this,
in `src/framework/do_concurrent_compat.h`:

    #ifdef HAVE_FC_DO_CONCURRENT_LOCAL
    #define DO_LOCALITY(X) X
    #else
    #define DO_LOCALITY(X) ;
    #endif

used as

    do concurrent (i=is:ie, j=js:je) DO_LOCALITY(local(a, b))

The translators parse **unpreprocessed** source, so they see the macro call,
not the specifier. Without unwrapping, `local(a, b)` is invisible and the
emitted `!$omp` loop silently loses its `private(a, b)` — the variables become
shared across threads. That is a wrong-answer bug, not a cosmetic one, so
unwrapping happens before any clause parsing.

The macro name is configurable: pass `names=`, or set the environment variable
DCPORT_LOCALITY_MACROS to a comma-separated list. Default: DO_LOCALITY.

The false-branch expansion (a bare `;`, an empty statement that keeps the line
syntactically valid) is also removed, so a header preprocessed by the
*unsupported* path parses identically.
"""
from __future__ import annotations

import os
import re

DEFAULT_MACROS = ("DO_LOCALITY",)


def configured_macros() -> tuple[str, ...]:
    """Locality-macro names, from DCPORT_LOCALITY_MACROS or the default."""
    env = os.environ.get("DCPORT_LOCALITY_MACROS")
    if not env:
        return DEFAULT_MACROS
    names = tuple(n.strip() for n in env.split(",") if n.strip())
    return names or DEFAULT_MACROS


def unwrap_locality_macros(text: str, names: tuple[str, ...] | None = None) -> str:
    """Replace `NAME(<args>)` with `<args>` for every configured macro name.

    Balanced-paren matching, so a nested call such as
    `DO_LOCALITY(local(a, b) shared(c))` unwraps to `local(a, b) shared(c)`.
    Applied repeatedly to a fixpoint so nested macro calls resolve. Also drops
    the bare `;` the false branch expands to. Case-sensitive: CPP macro names
    are, even though Fortran identifiers are not.
    """
    if names is None:
        names = configured_macros()
    for name in names:
        pat = re.compile(r"\b" + re.escape(name) + r"\s*\(")
        while True:
            m = pat.search(text)
            if not m:
                break
            depth, end = 0, -1
            for k in range(m.end() - 1, len(text)):
                ch = text[k]
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        end = k
                        break
            if end < 0:
                # Unbalanced (a continuation the caller failed to join).
                # Leave it alone rather than corrupt the line.
                break
            text = text[: m.start()] + text[m.end(): end] + text[end + 1:]
    # The `#define DO_LOCALITY(X) ;` branch: an empty statement carrying no
    # meaning for us. Strip it ONLY at end of text, which is where a locality
    # specifier would have stood. Fortran also uses `;` as a statement
    # separator, and MOM6 leans on that heavily (`do concurrent (i=is:ie) ;
    # press(i,j) = 0.0 ; enddo`). Removing THAT `;` would silently weld the
    # loop body onto the header and hide it from the trailing-statement guard
    # in the caller, so the strip has to be anchored.
    text = re.sub(r";\s*$", "", text)
    return text
