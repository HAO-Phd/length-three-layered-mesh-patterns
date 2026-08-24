#!/usr/bin/env python3
"""Deterministic lookup for length-three mesh-pattern classes on layered hosts.

The classification data are read from ``equidistribution_audit.json``, the
same 146-row audit used to produce Table 4 of the manuscript.  A raw mesh
pattern is supplied by its underlying permutation and its shaded cells.  The
script computes the active-cell signature, finds its Table 4 row, lists all
equivalent signatures, and can export every raw mesh diagram in that row.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


Cell = tuple[int, int]
Permutation = tuple[int, ...]
ALL_CELLS: tuple[Cell, ...] = tuple((x, y) for x in range(4) for y in range(4))
SIGNATURE_RE = re.compile(
    r"^r=\((?P<blocks>[^)]*)\);z=\((?P<z>[^)]*)\);E=\{(?P<external>[^}]*)\}$"
)


@dataclass(frozen=True, order=True)
class Signature:
    blocks: tuple[int, ...]
    z: tuple[int, ...]
    external: tuple[int, ...]

    def audit_key(self) -> str:
        blocks = ",".join(map(str, self.blocks))
        z_values = ",".join(map(str, self.z))
        external = ",".join(map(str, self.external))
        return f"r=({blocks});z=({z_values});E={{{external}}}"

    def table_label(self) -> str:
        blocks = "".join(map(str, self.blocks))
        z_values = "".join(map(str, self.z))
        external = ",".join(map(str, self.external))
        return f"({blocks};{z_values};{{{external}}})"


def _integer_tuple(text: str) -> tuple[int, ...]:
    text = text.strip()
    if not text:
        return ()
    return tuple(int(part.strip()) for part in text.split(","))


def parse_signature(label: str) -> Signature:
    match = SIGNATURE_RE.fullmatch(label)
    if match is None:
        raise ValueError(f"invalid audit signature: {label!r}")
    signature = Signature(
        _integer_tuple(match.group("blocks")),
        _integer_tuple(match.group("z")),
        _integer_tuple(match.group("external")),
    )
    if len(signature.blocks) != len(signature.z):
        raise ValueError(f"block/z length mismatch in {label!r}")
    return signature


def load_audit(path: Path) -> tuple[dict[int, dict], dict[str, int]]:
    try:
        records = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise SystemExit(f"classification data not found: {path}") from exc
    if not isinstance(records, list):
        raise ValueError("the audit JSON root must be a list")

    rows: dict[int, dict] = {}
    signature_to_class: dict[str, int] = {}
    for record in records:
        class_id = int(record["distribution_id"])
        members = list(record["members"])
        if class_id in rows:
            raise ValueError(f"duplicate distribution class {class_id}")
        if int(record["signature_count"]) != len(members):
            raise ValueError(f"signature count mismatch in class {class_id}")
        for member in members:
            parse_signature(member)
            if member in signature_to_class:
                raise ValueError(f"signature occurs in two rows: {member}")
            signature_to_class[member] = class_id
        rows[class_id] = record

    if sorted(rows) != list(range(1, 147)):
        raise ValueError("the audit must contain Table 4 classes 1 through 146")
    if len(signature_to_class) != 644:
        raise ValueError(
            f"the audit must contain 644 signatures, found {len(signature_to_class)}"
        )
    return rows, signature_to_class


def parse_permutation(text: str) -> Permutation:
    values = tuple(int(value) for value in re.findall(r"\d+", text))
    if len(values) == 1 and len(str(values[0])) == 3:
        values = tuple(int(char) for char in str(values[0]))
    if len(values) != 3 or sorted(values) != [1, 2, 3]:
        raise argparse.ArgumentTypeError(
            "the underlying permutation must be 123, 132, 213, 231, 312, or 321"
        )
    return values


def parse_cells(text: str) -> frozenset[Cell]:
    if not text.strip():
        return frozenset()
    pair_re = re.compile(r"\(?\s*(\d+)\s*[,/]\s*(\d+)\s*\)?")
    cells = [(int(x), int(y)) for x, y in pair_re.findall(text)]
    residue = pair_re.sub("", text)
    residue = re.sub(r"[\s,;{}\[\]()]", "", residue)
    if residue or not cells:
        raise argparse.ArgumentTypeError(
            "shaded cells must look like '0,0;1,1;3,1'"
        )
    for cell in cells:
        if not (0 <= cell[0] <= 3 and 0 <= cell[1] <= 3):
            raise argparse.ArgumentTypeError(
                f"cell {cell} is outside the length-three grid [0,3]^2"
            )
    return frozenset(cells)


def parse_mask(text: str) -> int:
    token = text.strip().lower().replace("_", "")
    base = 16 if token.startswith("0x") or re.search(r"[a-f]", token) else 10
    try:
        value = int(token, base)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid shading mask: {text!r}") from exc
    if not 0 <= value < (1 << 16):
        raise argparse.ArgumentTypeError("the shading mask must be between 0 and 0xffff")
    return value


def cells_to_mask(cells: Iterable[Cell]) -> int:
    mask = 0
    for x, y in cells:
        mask |= 1 << (4 * x + y)
    return mask


def mask_to_cells(mask: int) -> frozenset[Cell]:
    return frozenset(cell for bit, cell in enumerate(ALL_CELLS) if mask & (1 << bit))


def cells_text(cells: Iterable[Cell]) -> str:
    return ";".join(f"{x},{y}" for x, y in sorted(cells))


def decreasing_blocks(tau: Permutation) -> tuple[int, ...] | None:
    blocks: list[int] = []
    start = 0
    while start < len(tau):
        end = start + 1
        while end < len(tau) and tau[end - 1] == tau[end] + 1:
            end += 1
        block = tau[start:end]
        expected = tuple(range(start + len(block), start, -1))
        if block != expected:
            return None
        blocks.append(len(block))
        start = end
    return tuple(blocks)


def permutation_from_blocks(blocks: Sequence[int]) -> Permutation:
    permutation: list[int] = []
    offset = 0
    for size in blocks:
        permutation.extend(range(offset + size, offset, -1))
        offset += size
    return tuple(permutation)


def active_cell_data(
    blocks: Sequence[int],
) -> tuple[frozenset[Cell], tuple[tuple[Cell, ...], ...], tuple[Cell, ...]]:
    channels: list[tuple[Cell, ...]] = []
    offset = 0
    for size in blocks:
        channels.append(
            tuple((offset + h, offset + size - h) for h in range(size + 1))
        )
        offset += size
    boundaries: list[Cell] = [(0, 0)]
    offset = 0
    for size in blocks:
        offset += size
        boundaries.append((offset, offset))
    active = set(boundaries)
    for channel in channels:
        active.update(channel)
    return frozenset(active), tuple(channels), tuple(boundaries)


def signature_from_pattern(
    tau: Permutation, shading: frozenset[Cell]
) -> tuple[Signature | None, frozenset[Cell], frozenset[Cell]]:
    blocks = decreasing_blocks(tau)
    if blocks is None:
        return None, frozenset(), shading
    active, channels, boundaries = active_cell_data(blocks)
    z_values = tuple(sum(cell in shading for cell in channel) for channel in channels)
    external = tuple(i for i, cell in enumerate(boundaries) if cell in shading)
    signature = Signature(blocks, z_values, external)
    return signature, shading & active, shading - active


def raw_count_for_signature(signature: Signature) -> int:
    active, channels, _ = active_cell_data(signature.blocks)
    count = 1 << (16 - len(active))
    for channel, z_value in zip(channels, signature.z):
        count *= math.comb(len(channel), z_value)
    return count


def class_signatures(record: dict) -> tuple[Signature, ...]:
    return tuple(parse_signature(member) for member in record["members"])


def class_raw_count(record: dict) -> int:
    return sum(raw_count_for_signature(sig) for sig in class_signatures(record))


def lookup_pattern(
    tau: Permutation,
    shading: frozenset[Cell],
    rows: dict[int, dict],
    signature_to_class: dict[str, int],
) -> dict:
    signature, active_shaded, inactive_shaded = signature_from_pattern(tau, shading)
    base = {
        "underlying_permutation": "".join(map(str, tau)),
        "shading_mask": f"0x{cells_to_mask(shading):04x}",
        "shaded_cells": [list(cell) for cell in sorted(shading)],
    }
    if signature is None:
        return {
            **base,
            "zero_class": True,
            "table4_class": None,
            "message": "Underlying patterns 231 and 312 form the zero class, which is not numbered in Table 4.",
            "raw_patterns_in_class": 2 * (1 << 16),
        }

    key = signature.audit_key()
    try:
        class_id = signature_to_class[key]
    except KeyError as exc:
        raise RuntimeError(f"signature missing from the canonical audit: {key}") from exc
    record = rows[class_id]
    equivalent = class_signatures(record)
    return {
        **base,
        "zero_class": False,
        "table4_class": class_id,
        "block_word": "".join(map(str, signature.blocks)),
        "signature": signature.table_label(),
        "audit_key": key,
        "active_shaded_cells": [list(cell) for cell in sorted(active_shaded)],
        "inactive_shaded_cells": [list(cell) for cell in sorted(inactive_shaded)],
        "raw_patterns_for_signature": raw_count_for_signature(signature),
        "raw_patterns_in_class": class_raw_count(record),
        "equivalent_signatures": [sig.table_label() for sig in equivalent],
    }


def class_summary(class_token: str, rows: dict[int, dict]) -> dict:
    if class_token.lower() in {"zero", "z", "0"}:
        return {
            "zero_class": True,
            "table4_class": None,
            "message": "Underlying patterns 231 and 312 form the zero class, which is not numbered in Table 4.",
            "raw_patterns_in_class": 2 * (1 << 16),
        }
    try:
        class_id = int(class_token)
    except ValueError as exc:
        raise SystemExit("--class-id must be an integer from 1 to 146, or 'zero'") from exc
    if class_id not in rows:
        raise SystemExit("--class-id must be an integer from 1 to 146, or 'zero'")
    record = rows[class_id]
    signatures = class_signatures(record)
    return {
        "zero_class": False,
        "table4_class": class_id,
        "signature_count": len(signatures),
        "raw_patterns_in_class": class_raw_count(record),
        "equivalent_signatures": [sig.table_label() for sig in signatures],
    }


def raw_rows_for_signature(signature: Signature, class_id: int) -> Iterator[dict]:
    active, channels, boundaries = active_cell_data(signature.blocks)
    inactive = tuple(cell for cell in ALL_CELLS if cell not in active)
    channel_choices = [
        tuple(itertools.combinations(channel, z_value))
        for channel, z_value in zip(channels, signature.z)
    ]
    external_cells = {boundaries[i] for i in signature.external}
    permutation = "".join(map(str, permutation_from_blocks(signature.blocks)))

    for chosen_channels in itertools.product(*channel_choices):
        active_shading = set(external_cells)
        for chosen in chosen_channels:
            active_shading.update(chosen)
        for inactive_mask in range(1 << len(inactive)):
            shading = set(active_shading)
            for bit, cell in enumerate(inactive):
                if inactive_mask & (1 << bit):
                    shading.add(cell)
            yield {
                "table4_class": class_id,
                "permutation": permutation,
                "shading_mask": f"0x{cells_to_mask(shading):04x}",
                "shaded_cells": cells_text(shading),
                "signature": signature.table_label(),
            }


def raw_rows_for_class(class_token: str, rows: dict[int, dict]) -> Iterator[dict]:
    if class_token.lower() in {"zero", "z", "0"}:
        for permutation in ("231", "312"):
            for mask in range(1 << 16):
                shading = mask_to_cells(mask)
                yield {
                    "table4_class": "zero",
                    "permutation": permutation,
                    "shading_mask": f"0x{mask:04x}",
                    "shaded_cells": cells_text(shading),
                    "signature": "zero",
                }
        return
    class_id = int(class_token)
    for signature in class_signatures(rows[class_id]):
        yield from raw_rows_for_signature(signature, class_id)


def export_raw_patterns(records: Iterable[dict], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    if path.suffix.lower() == ".csv":
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "table4_class",
                    "permutation",
                    "shading_mask",
                    "shaded_cells",
                    "signature",
                ),
            )
            writer.writeheader()
            for record in records:
                writer.writerow(record)
                count += 1
    elif path.suffix.lower() in {".jsonl", ".ndjson"}:
        with path.open("w", encoding="utf-8") as handle:
            for record in records:
                handle.write(json.dumps(record, ensure_ascii=False) + "\n")
                count += 1
    else:
        raise SystemExit("--export-raw requires a .csv, .jsonl, or .ndjson path")
    return count


def print_human(summary: dict) -> None:
    if "underlying_permutation" in summary:
        print(f"Underlying permutation: {summary['underlying_permutation']}")
        print(f"Shading mask: {summary['shading_mask']}")
        print(
            "Shaded cells: "
            + (cells_text(tuple(map(tuple, summary["shaded_cells"]))) or "none")
        )
    if summary["zero_class"]:
        print("Distribution class: zero (not numbered in Table 4)")
        print(summary["message"])
        print(f"Raw mesh patterns in class: {summary['raw_patterns_in_class']}")
        return

    print(f"Table 4 class: {summary['table4_class']}")
    if "signature" in summary:
        print(f"Active signature: {summary['signature']}")
        print(f"Audit key: {summary['audit_key']}")
        inactive = tuple(map(tuple, summary["inactive_shaded_cells"]))
        print(f"Inactive shaded cells ignored: {cells_text(inactive) or 'none'}")
        print(
            "Raw mesh patterns represented by this signature: "
            f"{summary['raw_patterns_for_signature']}"
        )
    print(f"Raw mesh patterns in this class: {summary['raw_patterns_in_class']}")
    signatures = summary["equivalent_signatures"]
    print(f"Equivalent signatures ({len(signatures)}):")
    for signature in signatures:
        print(f"  {signature}")


def self_test(rows: dict[int, dict], signature_to_class: dict[str, int]) -> dict:
    layered_total = sum(class_raw_count(record) for record in rows.values())
    if layered_total != 4 * (1 << 16):
        raise AssertionError(
            f"raw expansion should contain 262144 layered diagrams, found {layered_total}"
        )

    seen: set[str] = set()
    for tau in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (3, 2, 1)):
        for mask in range(1 << 16):
            signature, _, _ = signature_from_pattern(tau, mask_to_cells(mask))
            assert signature is not None
            key = signature.audit_key()
            if key not in signature_to_class:
                raise AssertionError(f"unmapped signature: {key}")
            seen.add(key)
    if len(seen) != 644:
        raise AssertionError(f"raw diagrams should realize 644 signatures, found {len(seen)}")

    example = lookup_pattern(
        (1, 3, 2),
        frozenset({(0, 0), (1, 1), (1, 2), (3, 1)}),
        rows,
        signature_to_class,
    )
    if example["table4_class"] != 19 or example["signature"] != "(12;01;{0,1})":
        raise AssertionError("the manuscript's worked example does not resolve to Class 19")

    return {
        "status": "ok",
        "table4_classes": len(rows),
        "signatures": len(signature_to_class),
        "layered_raw_patterns": layered_total,
        "zero_class_raw_patterns": 2 * (1 << 16),
        "all_raw_patterns": 6 * (1 << 16),
        "worked_example_class": example["table4_class"],
        "worked_example_signature": example["signature"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Look up length-three mesh-pattern distribution classes on layered permutations."
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=Path(__file__).with_name("equidistribution_audit.json"),
        help="path to equidistribution_audit.json (default: beside this script)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--perm", type=parse_permutation, help="underlying permutation, e.g. 132")
    mode.add_argument(
        "--class-id", help="Table 4 class 1--146, or 'zero', to inspect directly"
    )
    mode.add_argument("--self-test", action="store_true", help="run the exhaustive consistency audit")
    shading = parser.add_mutually_exclusive_group()
    shading.add_argument(
        "--shade",
        type=parse_cells,
        help="shaded cells, e.g. '0,0;1,1;1,2;3,1' (default: none)",
    )
    shading.add_argument(
        "--mask",
        type=parse_mask,
        help="16-bit shading mask; bit 4*x+y represents cell (x,y)",
    )
    parser.add_argument(
        "--export-raw",
        type=Path,
        help="export every raw pattern in the resulting class to .csv or .jsonl",
    )
    parser.add_argument("--json", action="store_true", help="print the lookup result as JSON")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    rows, signature_to_class = load_audit(args.data)

    if args.self_test:
        if args.shade is not None or args.mask is not None or args.export_raw is not None:
            parser.error("--self-test cannot be combined with shading or export options")
        result = self_test(rows, signature_to_class)
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            for key, value in result.items():
                print(f"{key}: {value}")
        return 0

    if args.class_id is not None:
        if args.shade is not None or args.mask is not None:
            parser.error("--shade and --mask are available only with --perm")
        summary = class_summary(args.class_id, rows)
        export_token = args.class_id
    else:
        if args.mask is not None:
            shading = mask_to_cells(args.mask)
        elif args.shade is not None:
            shading = args.shade
        else:
            shading = frozenset()
        summary = lookup_pattern(args.perm, shading, rows, signature_to_class)
        export_token = "zero" if summary["zero_class"] else str(summary["table4_class"])

    if args.export_raw is not None:
        exported = export_raw_patterns(raw_rows_for_class(export_token, rows), args.export_raw)
        summary["exported_raw_patterns"] = exported
        summary["export_path"] = str(args.export_raw.resolve())

    if args.json:
        print(json.dumps(summary, indent=2, ensure_ascii=False))
    else:
        print_human(summary)
        if args.export_raw is not None:
            print(f"Exported raw patterns: {summary['exported_raw_patterns']}")
            print(f"Export path: {summary['export_path']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error
