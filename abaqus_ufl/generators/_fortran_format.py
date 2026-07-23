"""Shared Fortran line-wrapping utilities."""

import re


# Token regex for Fortran literals and operators
# Matches: quoted strings, scientific literals (1.23d-05, 4.5e+10), power operator (**)
_LITERAL_RE = re.compile(r"'[^']*'|\"[^\"]*\"|\d+\.?\d*[dDeE][+-]?\d+|\*\*")


def _is_safe_break(line, pos):
    """
    Return True if splitting *after* `pos` (i.e. at pos+1) is safe.

    Unsafe splits:
      - Inside a quoted string
      - Inside a numeric literal with exponent (e.g. 1.234d-05)
      - Between the two stars of **
    """
    # Check against literal tokens
    for m in _LITERAL_RE.finditer(line):
        start, end = m.span()
        # If pos is strictly inside the token, or pos+1 splits the token, it's unsafe
        if start <= pos < end:
            return False
        # Also unsafe if break is exactly at token start (would split token)
        if pos + 1 == start:
            # This is okay - breaking before token is fine
            pass
        # Unsafe if break is inside token
        if start < pos + 1 <= end:
            return False
    return True


def wrap_fortran_line(line, max_col=72):
    """
    Wrap a single Fortran fixed-format line at column max_col.

    Comment lines (column 1 is C, c, !, or *) are NOT wrapped.
    Continuation lines use '     &' in columns 1-6.

    Never splits inside numeric literals (e.g. 1.234d-05), quoted strings,
    or the ** operator. Prefers breaks after commas, then spaces, then
    other operators at top-level.

    Raises RuntimeError if no safe break point exists.
    """
    # Comment or empty — pass through
    if not line or line[0] in ('C', 'c', '!', '*'):
        return line.rstrip()

    if len(line) <= max_col:
        return line.rstrip()

    result_lines = []
    remaining = line

    while len(remaining) > max_col:
        limit = max_col

        # Search backwards from limit for a safe break character
        break_chars = ' ,+-*/'
        break_pos = -1
        for pos in range(min(limit, len(remaining)) - 1, 5, -1):
            if remaining[pos] in break_chars and _is_safe_break(remaining, pos):
                break_pos = pos + 1  # break AFTER the character
                break

        if break_pos <= 6:
            # No safe break found — this is a genuine error
            raise RuntimeError(
                f"Cannot wrap Fortran line safely at column {max_col}: "
                f"{remaining[:max_col+20]!r}..."
            )

        result_lines.append(remaining[:break_pos].rstrip())
        remaining = '     &' + remaining[break_pos:]

    if remaining:
        result_lines.append(remaining.rstrip())

    return '\n'.join(result_lines)


def wrap_fortran_lines(lines, max_col=72):
    """Wrap a sequence of Fortran lines, returning a new list."""
    return [wrap_fortran_line(line, max_col) for line in lines]
