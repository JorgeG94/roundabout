#!/usr/bin/env python3
"""Commit-msg guard: no Claude session links in commit messages.

Session URLs (`https://claude.ai/code/session_...`) and the matching
`Claude-Session:` trailer are per-conversation handles.  They are useless to
anyone reading `git log` (nobody else can open them), they never expire out of
the history, and they leak which assistant session produced a change into a
permanent, public record.  The maintainer's rule is simple: they do not belong
in a commit message.

`Co-Authored-By:` is explicitly NOT affected — attribution trailers are fine,
it is only the session handle that is banned.

Wired as a `local` hook with `stages: [commit-msg]`, so pre-commit hands us the
path of the prepared message file (`.git/COMMIT_EDITMSG`) as argv[1].

Run the embedded cases with:

    python3 tools/check_commit_msg_no_session_links.py --self-test
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# A session URL in any form: bare, in angle brackets, http/https, with or
# without a trailing path.  We match the distinctive middle of the URL so the
# scheme and host prefix ("www.", …) do not matter.
SESSION_URL_RE = re.compile(r"claude\.ai/code/session", re.IGNORECASE)

# The trailer form, which can appear without a URL at all.
SESSION_TRAILER_RE = re.compile(r"^\s*Claude-Session\s*:", re.IGNORECASE)

# `git commit --verbose` appends the staged diff below this marker; git strips
# everything from it downwards before the message is stored, so neither it nor
# the comment lines are part of the commit.  Scanning them would fail a commit
# whose *diff* legitimately touches this very file.
SCISSORS = "# ------------------------ >8 ------------------------"


def offending_lines(message: str) -> list[tuple[int, str]]:
    """Return `(1-based line number, line)` for every banned line.

    Comment lines and everything below the `--verbose` scissors marker are
    skipped: git discards them, so they never reach the stored message.
    """
    hits: list[tuple[int, str]] = []
    for lineno, line in enumerate(message.splitlines(), start=1):
        if line.rstrip() == SCISSORS:
            break
        if line.startswith("#"):
            continue
        if SESSION_URL_RE.search(line) or SESSION_TRAILER_RE.match(line):
            hits.append((lineno, line.strip()))
    return hits


def check_file(path: Path) -> int:
    try:
        message = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        print(f"check_commit_msg_no_session_links: cannot read {path}: {exc}")
        return 1

    hits = offending_lines(message)
    if not hits:
        return 0

    print("check_commit_msg_no_session_links: FAIL")
    for lineno, line in hits:
        print(f"  - line {lineno}: {line}")
    print(
        "\nCommit messages must not carry Claude session links "
        "(`claude.ai/code/session...` URLs or a `Claude-Session:` trailer) — "
        "they are per-conversation handles nobody else can open, and they are "
        "permanent once committed.  Drop those lines and commit again.\n"
        "`Co-Authored-By:` trailers are fine and are not what this rejects."
    )
    return 1


def self_test() -> int:
    """Tiny in-file test suite — no pytest, no third-party imports."""
    cases: list[tuple[str, str, bool]] = [
        # (name, message, should_fail)
        (
            "clean message with Co-Authored-By",
            "feat: add a thing\n\nBody text.\n\n"
            "Co-Authored-By: Someone <noreply@example.com>\n",
            False,
        ),
        (
            "Claude-Session trailer",
            "fix: a bug\n\nClaude-Session: https://claude.ai/code/session_abc123\n",
            True,
        ),
        (
            "Claude-Session trailer without a URL",
            "fix: a bug\n\nClaude-Session: session_abc123\n",
            True,
        ),
        (
            "lowercase trailer spelling",
            "fix: a bug\n\nclaude-session: whatever\n",
            True,
        ),
        (
            "bare session URL in the body",
            "docs: notes\n\nSee https://claude.ai/code/session_xyz for context.\n",
            True,
        ),
        (
            "uppercase host",
            "docs: notes\n\nSee HTTPS://CLAUDE.AI/CODE/SESSION_XYZ\n",
            True,
        ),
        (
            "session URL only in a comment line",
            "chore: tidy\n\n# https://claude.ai/code/session_abc\n",
            False,
        ),
        (
            "session URL only below the --verbose scissors",
            "chore: tidy\n\n"
            f"{SCISSORS}\n"
            "diff --git a/x b/x\n"
            "+Claude-Session: https://claude.ai/code/session_abc\n",
            False,
        ),
        (
            "unrelated claude.ai link is allowed",
            "docs: link the docs\n\nSee https://claude.ai/code for the CLI.\n",
            False,
        ),
        (
            "Co-Authored-By naming Claude is allowed",
            "feat: x\n\nCo-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>\n",
            False,
        ),
    ]

    failures = 0
    for name, message, should_fail in cases:
        got = bool(offending_lines(message))
        if got != should_fail:
            failures += 1
            print(
                f"  FAIL {name}: expected "
                f"{'rejection' if should_fail else 'acceptance'}, got the opposite"
            )

    if failures:
        print(f"check_commit_msg_no_session_links self-test: {failures} FAILED")
        return 1
    print(f"check_commit_msg_no_session_links self-test: OK ({len(cases)} cases)")
    return 0


def main(argv: list[str]) -> int:
    if "--self-test" in argv:
        return self_test()
    if not argv:
        print(
            "usage: check_commit_msg_no_session_links.py <commit-msg-file> "
            "| --self-test"
        )
        return 1
    return max(check_file(Path(a)) for a in argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
