#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


ALIASES = {
    "TRAI": "TRA1",
    "TRA21": "TRA2",
}


def tokenize(line: str) -> list[str]:
    return [t for t in re.split(r"[,\s]+", line.strip()) if t]


def clean_numeric_token(token: str) -> str:
    token = token.replace("−", "-")
    return re.sub(r"[^0-9eE+\-\.]", "", token)


def parse_number(token: str):
    t = clean_numeric_token(token)
    if not t or t in {"+", "-", ".", "+.", "-."}:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def to_f5_field(token: str) -> str:
    value = parse_number(token)
    if value is None:
        raise ValueError(f"Cannot parse numeric token '{token}'")

    # Prefer decimal literals if they fit in 5 chars; otherwise use implied-decimal style.
    for dec in (4, 3, 2, 1, 0):
        s = f"{value:.{dec}f}".rstrip("0").rstrip(".")
        if s in {"", "-0"}:
            s = "0"
        if len(s) <= 5:
            return s.rjust(5)

    scaled = int(round(value * 100.0))
    s = str(scaled)
    if len(s) <= 5:
        return s.rjust(5)

    raise ValueError(f"Value '{token}' does not fit F5 field")


def to_scaled_field(token: str, scale: int) -> str:
    value = parse_number(token)
    if value is None:
        raise ValueError(f"Cannot parse numeric token '{token}'")
    s = str(int(round(value * scale)))
    if len(s) > 5:
        raise ValueError(f"Scaled value '{s}' does not fit F5 field")
    return s.rjust(5)


def build_line(cmd: str, head: str, values: list[str], scales: list[int] | None = None) -> str:
    fields = []
    for i, tok in enumerate(values[:14]):
        if scales and i < len(scales):
            fields.append(to_scaled_field(tok, scales[i]))
        else:
            fields.append(to_f5_field(tok))
    line = f"{cmd:<4}{head}{''.join(fields)}"
    return line[:80].rstrip()


def convert_line(line: str, line_no: int, warnings: list[str]) -> str | None:
    stripped = line.strip()
    if not stripped:
        return ""
    if stripped.startswith(("#", "!")):
        return stripped

    tokens = tokenize(line)
    if not tokens:
        return ""

    raw_cmd = tokens[0].upper()
    cmd = ALIASES.get(raw_cmd, raw_cmd)

    if cmd == "REMO":
        warnings.append(f"line {line_no}: skipped REMO (not supported by current parser)")
        return None

    if cmd in {"TRA1", "TRA2"}:
        if len(tokens) < 2:
            raise ValueError(f"line {line_no}: {cmd} requires an airfoil ID token")
        airfoil_id = tokens[1]
        if len(airfoil_id) > 4:
            warnings.append(f"line {line_no}: airfoil id '{airfoil_id}' truncated to 4 chars")
            airfoil_id = airfoil_id[-4:]
        head = f"  {airfoil_id:>4}"
        if cmd == "TRA1":
            return build_line(cmd, head, tokens[2:], scales=[100] * 14)
        # TRA2 scaling mirrors book-style fixed examples (A4,3I1,I3,14F5.2 + legacy scale factors).
        tra2_scales = [100, 100, 100, 1000, 1000, 100, 100, 100, 1000, 1000, 100, 100, 100, 100]
        return build_line(cmd, head, tokens[2:], scales=tra2_scales)

    if cmd == "ALFA":
        if len(tokens) < 2:
            raise ValueError(f"line {line_no}: ALFA requires angle count followed by angle values")
        nupu = int(round(float(clean_numeric_token(tokens[1]))))
        head = f"000{nupu:>3}"
        return build_line(cmd, head, tokens[2:], scales=[100] * 14)

    if cmd == "ENDE":
        return "ENDE"

    warnings.append(f"line {line_no}: unsupported command '{raw_cmd}' copied as-is")
    return stripped[:80]


def convert_file(src: Path, dst: Path) -> list[str]:
    warnings: list[str] = []
    out_lines: list[str] = []
    for i, line in enumerate(src.read_text(encoding="utf-8").splitlines(), start=1):
        converted = convert_line(line, i, warnings)
        if converted is None:
            continue
        out_lines.append(converted)

    dst.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return warnings


def main():
    parser = argparse.ArgumentParser(
        description="Convert Eppler-style free-format deck lines to fixed-column lines."
    )
    parser.add_argument("input", help="Path to free-format input deck")
    parser.add_argument("output", help="Path to write converted fixed-format deck")
    args = parser.parse_args()

    warnings = convert_file(Path(args.input), Path(args.output))
    if warnings:
        print("Warnings:")
        for w in warnings:
            print(f"  - {w}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
