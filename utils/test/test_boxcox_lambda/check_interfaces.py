#!/usr/bin/env python3
"""Check error-last interfaces and direct calls in the available SUMMA source tree.

This focused style check supplements the compiled tests; it is not a Fortran
compiler and does not resolve renamed imports or indirect procedure calls.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

ERROR_NAMES = {
    "compute_metric": "err",
    "prepare_series": "err",
    "apply_transformation": "err",
    "write_evaluation": "err",
    "define_evaluation_target": "ierr",
}
PATTERN = re.compile(
    r"\b(subroutine|call)\s+(" + "|".join(ERROR_NAMES) + r")\s*\(", re.I
)


def without_comment(line: str) -> str:
    quote = None
    i = 0
    while i < len(line):
        char = line[i]
        if quote:
            if char == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    i += 2
                    continue
                quote = None
        elif char in "'\"":
            quote = char
        elif char == "!":
            return line[:i]
        i += 1
    return line


def statements(path: Path):
    pending = ""
    first = 1
    for number, line in enumerate(path.read_text().splitlines(), 1):
        code = without_comment(line).strip()
        if not code:
            continue
        if not pending:
            first = number
        elif code.startswith("&"):
            code = code[1:].lstrip()
        continued = code.endswith("&")
        pending += " " + (code[:-1] if continued else code)
        if not continued:
            yield first, pending
            pending = ""
    if pending:
        raise ValueError(f"{path}:{first}: unterminated continuation")


def arguments(text: str, start: int) -> list[str]:
    """Split one parenthesized list, respecting nested arrays, calls and strings."""
    depth = 0
    quote = None
    result = []
    field = start
    i = start
    while i < len(text):
        char = text[i]
        if quote:
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    i += 2
                    continue
                quote = None
        elif char in "'\"":
            quote = char
        elif char in "([":
            depth += 1
        elif char == ")" and depth == 0:
            return result + [text[field:i].strip()]
        elif char in ")]":
            depth -= 1
        elif char == "," and depth == 0:
            result.append(text[field:i].strip())
            field = i + 1
        i += 1
    raise ValueError("Unbalanced argument list")


def audit(root: Path) -> tuple[int, int]:
    definitions = calls = 0
    seen = set()
    failures = []
    paths = sorted({p for tree in (root / "build/source", root / "utils/test")
                    for p in tree.rglob("*") if p.suffix.lower() == ".f90"})
    for path in paths:
        for line, statement in statements(path):
            for match in PATTERN.finditer(statement):
                kind, name = (s.lower() for s in match.groups())
                args = arguments(statement, match.end())
                error_name = ERROR_NAMES[name]
                if kind == "subroutine":
                    definitions += 1
                    seen.add(name)
                    valid = [s.lower() for s in args[-3:]] == [
                        "boxcox_lambda", error_name, "message"
                    ]
                    requirement = "signature must end with boxcox_lambda," + error_name + ",message"
                else:
                    calls += 1
                    valid = len(args) >= 2 and all(
                        re.match(r"^" + key + r"\s*=", arg, re.I)
                        for key, arg in zip((error_name, "message"), args[-2:])
                    )
                    requirement = "call must end with keyword " + error_name + "=...,message=..."
                if not valid:
                    failures.append(f"{path.relative_to(root)}:{line}: {name}: {requirement}")
    if seen != set(ERROR_NAMES):
        failures.append("Missing definitions: " + ", ".join(sorted(set(ERROR_NAMES) - seen)))
    if calls == 0:
        failures.append("No affected calls found")
    if failures:
        raise ValueError("\n".join(failures))
    return definitions, calls


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    try:
        definitions, calls = audit(args.root.resolve())
    except (OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(f"PASS: error-last signatures and keyword calls ({definitions} definitions, {calls} call sites).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
