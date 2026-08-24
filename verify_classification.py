"""Exact reproducibility audit for the length-three layered classification.

The script uses integer/rational arithmetic only.  Run

    python verify_classification.py --write-certificates

once to regenerate the enriched orbit map, explicit coordinate witnesses, and
lower-bound certificate, then run it without arguments to verify the committed
files.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True, order=True)
class Signature:
    r: tuple[int, ...]
    z: tuple[int, ...]
    E: frozenset[int]

    @property
    def s(self) -> int:
        return len(self.r)

    @property
    def f(self) -> int:
        return self.s + 1 - len(self.E)

    @property
    def degree_multiset(self) -> tuple[int, ...]:
        return tuple(sorted(r - z for r, z in zip(self.r, self.z) if z <= r))

    def rc(self) -> "Signature":
        return Signature(
            self.r[::-1], self.z[::-1], frozenset(self.s - e for e in self.E)
        )

    def text(self) -> str:
        r = ",".join(map(str, self.r))
        z = ",".join(map(str, self.z))
        E = ",".join(map(str, sorted(self.E)))
        return f"r=({r});z=({z});E={{{E}}}"


def compositions(n: int):
    if n == 0:
        yield ()
        return
    for mask in range(1 << (n - 1)):
        parts: list[int] = []
        previous = 0
        for i in range(n - 1):
            if mask & (1 << i):
                parts.append(i + 1 - previous)
                previous = i + 1
        parts.append(n - previous)
        yield tuple(parts)


def signature_count(k: int) -> int:
    return sum(
        2 ** (len(r) + 1) * math.prod(part + 2 for part in r)
        for r in compositions(k)
    )


def all_signatures() -> list[Signature]:
    result: list[Signature] = []
    for r in compositions(3):
        s = len(r)
        for z in itertools.product(*(range(part + 2) for part in r)):
            for mask in range(1 << (s + 1)):
                E = frozenset(i for i in range(s + 1) if mask & (1 << i))
                result.append(Signature(r, tuple(z), E))
    return sorted(result)


def parse_signature(text: str) -> Signature:
    match = re.fullmatch(r"r=\(([^)]*)\);z=\(([^)]*)\);E=\{([^}]*)\}", text)
    if match is None:
        raise ValueError(f"invalid signature: {text}")
    r = tuple(int(x) for x in match.group(1).split(",") if x)
    z = tuple(int(x) for x in match.group(2).split(",") if x)
    E = frozenset(int(x) for x in match.group(3).split(",") if x)
    return Signature(r, z, E)


def kernel(r: int, z: int, ell: int) -> int:
    if ell < r:
        return 0
    if z == r + 1:
        return int(ell == r)
    return math.comb(ell - z, r - z)


# Coefficients in B=(1,x,C(x,2),C(x,3),1_{x=1},1_{x=2},1_{x=3}).
KERNEL_BASIS: dict[tuple[int, int], dict[int, int]] = {
    (1, 0): {1: 1},
    (1, 1): {0: 1},
    (1, 2): {4: 1},
    (2, 0): {2: 1},
    (2, 1): {0: -1, 1: 1},
    (2, 2): {0: 1, 4: -1},
    (2, 3): {5: 1},
    (3, 0): {3: 1},
    (3, 1): {0: 1, 1: -1, 2: 1},
    (3, 2): {0: -2, 1: 1, 4: 1},
    (3, 3): {0: 1, 4: -1, 5: -1},
    (3, 4): {6: 1},
}


def basis_value(letter: int, x: int) -> int:
    if letter == 0:
        return 1
    if letter == 1:
        return x
    if letter == 2:
        return math.comb(x, 2)
    if letter == 3:
        return math.comb(x, 3)
    return int(x == letter - 3)


def supports(signature: Signature, m: int):
    for phi in itertools.combinations(range(m), signature.s):
        if 0 in signature.E and phi[0] != 0:
            continue
        if signature.s in signature.E and phi[-1] != m - 1:
            continue
        if any(
            j in signature.E and phi[j] != phi[j - 1] + 1
            for j in range(1, signature.s)
        ):
            continue
        yield phi


def occurrence(signature: Signature, host: tuple[int, ...]) -> int:
    total = 0
    for phi in supports(signature, len(host)):
        weight = 1
        for j, position in enumerate(phi):
            weight *= kernel(signature.r[j], signature.z[j], host[position])
        total += weight
    return total


def tensor(signature: Signature, m: int) -> dict[tuple[int, ...], int]:
    result: dict[tuple[int, ...], int] = defaultdict(int)
    for phi in supports(signature, m):
        choices = [
            list(KERNEL_BASIS[(signature.r[j], signature.z[j])].items())
            for j in range(signature.s)
        ]
        for selected in itertools.product(*choices):
            word = [0] * m
            coefficient = 1
            for position, (letter, local_coefficient) in zip(phi, selected):
                word[position] = letter
                coefficient *= local_coefficient
            result[tuple(word)] += coefficient
    return {word: coefficient for word, coefficient in result.items() if coefficient}


def tensor_key(signature: Signature, maximum_m: int) -> tuple:
    return tuple(
        tuple(sorted(tensor(signature, m).items()))
        for m in range(1, maximum_m + 1)
    )


# Explicit representatives for the 43 coordinate edges in Table 2.  Each tuple
# is (final class, source local orbit, target local orbit, U-index,
# source signature, target signature).
COORDINATE_WITNESS_DATA = (
    (3, 1, 2, 1, "r=(1,1,1);z=(0,0,0);E={0,1}", "r=(1,1,1);z=(0,0,0);E={0,3}"),
    (11, 1, 2, 1, "r=(1,1,1);z=(0,0,1);E={0,1}", "r=(1,1,1);z=(0,1,0);E={0,3}"),
    (20, 1, 2, 1, "r=(1,1,1);z=(0,0,2);E={0,1}", "r=(1,1,1);z=(0,2,0);E={0,3}"),
    (12, 1, 2, 2, "r=(1,1,1);z=(0,0,1);E={0,2}", "r=(1,1,1);z=(0,0,1);E={2,3}"),
    (12, 1, 3, 3, "r=(1,1,1);z=(0,0,1);E={0,2}", "r=(1,1,1);z=(0,1,0);E={0,1}"),
    (12, 2, 3, 4, "r=(1,1,1);z=(0,0,1);E={2,3}", "r=(1,1,1);z=(0,1,0);E={0,1}"),
    (22, 1, 2, 2, "r=(1,1,1);z=(0,0,2);E={0,3}", "r=(1,1,1);z=(0,0,2);E={2,3}"),
    (22, 1, 3, 3, "r=(1,1,1);z=(0,0,2);E={0,3}", "r=(1,1,1);z=(0,2,0);E={0,1}"),
    (22, 2, 3, 4, "r=(1,1,1);z=(0,0,2);E={2,3}", "r=(1,1,1);z=(0,2,0);E={0,1}"),
    (43, 1, 2, 2, "r=(1,1,1);z=(0,1,2);E={0,3}", "r=(1,1,1);z=(1,0,2);E={2,3}"),
    (43, 1, 3, 3, "r=(1,1,1);z=(0,1,2);E={0,3}", "r=(1,1,1);z=(0,2,1);E={0,1}"),
    (43, 2, 3, 4, "r=(1,1,1);z=(1,0,2);E={2,3}", "r=(1,1,1);z=(0,2,1);E={0,1}"),
    (13, 1, 2, 8, "r=(1,1,1);z=(0,0,1);E={0}", "r=(1,1,1);z=(0,1,0);E={0}"),
    (21, 1, 2, 8, "r=(1,1,1);z=(0,0,2);E={0,2}", "r=(1,1,1);z=(0,2,0);E={0,2}"),
    (23, 1, 2, 8, "r=(1,1,1);z=(0,0,2);E={0}", "r=(1,1,1);z=(0,2,0);E={0}"),
    (44, 1, 2, 8, "r=(1,1,1);z=(0,1,2);E={0}", "r=(1,1,1);z=(0,2,1);E={0}"),
    (35, 1, 2, 5, "r=(1,1,1);z=(0,1,1);E={0,1}", "r=(1,1,1);z=(1,0,1);E={0,1}"),
    (37, 1, 2, 6, "r=(1,1,1);z=(0,1,1);E={1,2}", "r=(1,1,1);z=(1,0,1);E={0,2}"),
    (38, 1, 2, 6, "r=(1,1,1);z=(0,1,1);E={1}", "r=(1,1,1);z=(1,0,1);E={0}"),
    (65, 1, 2, 6, "r=(1,1,1);z=(0,2,2);E={2,3}", "r=(1,1,1);z=(2,0,2);E={0,3}"),
    (90, 1, 2, 6, "r=(1,1,1);z=(1,2,2);E={2,3}", "r=(1,1,1);z=(2,1,2);E={0,3}"),
    (102, 1, 2, 6, "r=(1,1,1);z=(2,2,2);E={2,3}", "r=(1,1,1);z=(2,2,2);E={0,3}"),
    (42, 1, 2, 5, "r=(1,1,1);z=(0,1,2);E={0,1}", "r=(1,1,1);z=(1,0,2);E={0,1}"),
    (42, 1, 3, 1, "r=(1,1,1);z=(0,1,2);E={0,1}", "r=(1,1,1);z=(0,2,1);E={0,2}"),
    (42, 2, 3, 7, "r=(1,1,1);z=(1,0,2);E={0,1}", "r=(1,1,1);z=(0,2,1);E={0,2}"),
    (46, 1, 2, 6, "r=(1,1,1);z=(0,1,2);E={1,3}", "r=(1,1,1);z=(1,0,2);E={0,3}"),
    (46, 1, 3, 11, "r=(1,1,1);z=(0,1,2);E={1,3}", "r=(1,1,1);z=(0,2,1);E={2,3}"),
    (46, 2, 3, 7, "r=(1,1,1);z=(1,0,2);E={0,3}", "r=(1,1,1);z=(0,2,1);E={2,3}"),
    (48, 1, 2, 9, "r=(1,1,1);z=(0,1,2);E={3}", "r=(1,1,1);z=(1,0,2);E={3}"),
    (63, 1, 2, 9, "r=(1,1,1);z=(0,2,2);E={1,3}", "r=(1,1,1);z=(2,0,2);E={1,3}"),
    (67, 1, 2, 9, "r=(1,1,1);z=(0,2,2);E={3}", "r=(1,1,1);z=(2,0,2);E={3}"),
    (92, 1, 2, 9, "r=(1,1,1);z=(1,2,2);E={3}", "r=(1,1,1);z=(2,1,2);E={3}"),
    (53, 1, 2, 7, "r=(1,1,1);z=(1,0,2);E={0,2}", "r=(1,1,1);z=(0,2,1);E={1,2}"),
    (55, 1, 2, 7, "r=(1,1,1);z=(1,0,2);E={0}", "r=(1,1,1);z=(0,2,1);E={2}"),
    (59, 1, 2, 3, "r=(1,1,1);z=(0,2,2);E={0,3}", "r=(1,1,1);z=(0,2,2);E={0,1}"),
    (59, 1, 3, 2, "r=(1,1,1);z=(0,2,2);E={0,3}", "r=(1,1,1);z=(2,0,2);E={2,3}"),
    (59, 2, 3, 10, "r=(1,1,1);z=(0,2,2);E={0,1}", "r=(1,1,1);z=(2,0,2);E={2,3}"),
    (87, 1, 2, 3, "r=(1,1,1);z=(1,2,2);E={0,3}", "r=(1,1,1);z=(1,2,2);E={0,1}"),
    (87, 1, 3, 2, "r=(1,1,1);z=(1,2,2);E={0,3}", "r=(1,1,1);z=(2,1,2);E={1,3}"),
    (87, 2, 3, 10, "r=(1,1,1);z=(1,2,2);E={0,1}", "r=(1,1,1);z=(2,1,2);E={1,3}"),
    (79, 1, 2, 1, "r=(1,1,1);z=(1,1,2);E={0,1}", "r=(1,1,1);z=(1,2,1);E={0,2}"),
    (81, 1, 2, 1, "r=(1,1,1);z=(1,1,2);E={0}", "r=(1,1,1);z=(1,2,1);E={2}"),
    (80, 1, 2, 11, "r=(1,1,1);z=(1,1,2);E={0,3}", "r=(1,1,1);z=(1,2,1);E={2,3}"),
)


def coordinate_permutation(index: int, m: int) -> tuple[int, ...]:
    if m < 3:
        return tuple(range(m))
    permutations = {
        1: (0, *range(2, m), 1),
        2: (*range(1, m - 1), 0, m - 1),
        3: (0, m - 1, *range(1, m - 1)),
        4: (m - 2, m - 1, *range(0, m - 2)),
        5: (1, 0, *range(2, m)),
        6: (m - 2, *range(0, m - 2), m - 1),
        7: (*range(1, m), 0),
        8: (0, m - 1, *range(m - 2, 0, -1)),
        9: (*range(m - 2, -1, -1), m - 1),
        10: (*range(2, m), 0, 1),
        11: (*range(0, m - 2), m - 1, m - 2),
    }
    result = permutations[index]
    assert tuple(sorted(result)) == tuple(range(m))
    return result


def relabel_tensor_entries(entries, m: int, index: int):
    """Express a tensor evaluated at U_index(L) in the coordinates of L."""
    permutation = coordinate_permutation(index, m)
    result = {}
    for word, coefficient in entries.items():
        relabelled = [0] * m
        for output_position, input_position in enumerate(permutation):
            relabelled[input_position] = word[output_position]
        result[tuple(relabelled)] = coefficient
    return result


def build_coordinate_witness_certificate():
    witnesses = []
    for class_id, source_orbit, target_orbit, index, source, target in (
        COORDINATE_WITNESS_DATA
    ):
        witnesses.append(
            {
                "distribution_class": class_id,
                "source_local_orbit": source_orbit,
                "target_local_orbit": target_orbit,
                "coordinate_map": f"U_{index}",
                "source_signature": source,
                "target_signature": target,
                "identity": f"P_source(L) = P_target(U_{index}(L))",
            }
        )
    return {
        "description": (
            "Explicit signature witnesses and orientations for all 43 "
            "layer-coordinate edges in Table 2."
        ),
        "identity_convention": (
            "For each record and every positive composition L, the source "
            "statistic at L equals the target statistic at U_i(L)."
        ),
        "witness_count": len(witnesses),
        "witnesses": witnesses,
    }


def verify_coordinate_witnesses(certificate, distributions_raw) -> None:
    assert certificate["witness_count"] == 43
    assert len(certificate["witnesses"]) == 43
    by_class = {item["distribution_id"]: item for item in distributions_raw}
    seen = set()
    tensor_cache = {}

    def cached_tensor(signature: Signature, m: int):
        key = (signature, m)
        if key not in tensor_cache:
            tensor_cache[key] = tensor(signature, m)
        return tensor_cache[key]

    for witness in certificate["witnesses"]:
        class_id = witness["distribution_class"]
        source_label = witness["source_local_orbit"]
        target_label = witness["target_local_orbit"]
        index = int(witness["coordinate_map"].removeprefix("U_"))
        edge = (class_id, source_label, target_label, index)
        assert edge not in seen
        seen.add(edge)

        item = by_class[class_id]
        local_orbits = {
            orbit["local_label"]: orbit for orbit in item["local_orbits"]
        }
        source = parse_signature(witness["source_signature"])
        target = parse_signature(witness["target_signature"])
        assert source.text() in local_orbits[source_label]["members"]
        assert target.text() in local_orbits[target_label]["members"]
        assert source.s == target.s == 3
        assert witness["identity"] == f"P_source(L) = P_target(U_{index}(L))"

        for m in range(1, 31):
            target_tensor = cached_tensor(target, m)
            relabelled = relabel_tensor_entries(target_tensor, m, index)
            assert cached_tensor(source, m) == relabelled

def linear_representation(signature: Signature):
    """Return the seven exact letter matrices for the full tensor series."""
    dimension = signature.s + 1
    matrices = [
        [[0 for _ in range(dimension)] for _ in range(dimension)]
        for _ in range(7)
    ]
    for state in range(dimension):
        if state not in signature.E:
            matrices[0][state][state] += 1
    for j, (r, z) in enumerate(zip(signature.r, signature.z), start=1):
        for letter, coefficient in KERNEL_BASIS[(r, z)].items():
            matrices[letter][j - 1][j] += coefficient
    initial = [Fraction(0) for _ in range(dimension)]
    final = [Fraction(0) for _ in range(dimension)]
    initial[0] = Fraction(1)
    final[-1] = Fraction(1)
    return matrices, initial, final


def multiply_row(row: list[Fraction], matrix: list[list[int]]) -> list[Fraction]:
    return [
        sum((row[i] * matrix[i][j] for i in range(len(row))), Fraction(0))
        for j in range(len(row))
    ]


class RowSpace:
    def __init__(self) -> None:
        self._rows: dict[int, list[Fraction]] = {}

    def reduce(self, row: list[Fraction]) -> list[Fraction]:
        value = list(row)
        for pivot in sorted(self._rows):
            if value[pivot]:
                factor = value[pivot]
                base = self._rows[pivot]
                value = [x - factor * y for x, y in zip(value, base)]
        return value

    def add(self, row: list[Fraction]) -> bool:
        value = self.reduce(row)
        pivot = next((i for i, x in enumerate(value) if x), None)
        if pivot is None:
            return False
        scale = value[pivot]
        value = [x / scale for x in value]
        for old_pivot, old_row in list(self._rows.items()):
            if old_row[pivot]:
                factor = old_row[pivot]
                self._rows[old_pivot] = [
                    x - factor * y for x, y in zip(old_row, value)
                ]
        self._rows[pivot] = value
        return True

    @property
    def rows(self) -> list[list[Fraction]]:
        return [self._rows[pivot] for pivot in sorted(self._rows)]


def block_diagonal(left: list[list[int]], right: list[list[int]]):
    n, m = len(left), len(right)
    result = [[0 for _ in range(n + m)] for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            result[i][j] = left[i][j]
    for i in range(m):
        for j in range(m):
            result[n + i][n + j] = right[i][j]
    return result


def series_equal(left: Signature, right: Signature) -> bool:
    left_matrices, left_initial, left_final = linear_representation(left)
    right_matrices, right_initial, right_final = linear_representation(right)
    matrices = [
        block_diagonal(left_matrices[i], right_matrices[i]) for i in range(7)
    ]
    initial = left_initial + [-x for x in right_initial]
    final = left_final + right_final

    space = RowSpace()
    space.add(initial)
    while True:
        old_dimension = len(space.rows)
        for row in list(space.rows):
            for matrix in matrices:
                space.add(multiply_row(row, matrix))
        if len(space.rows) == old_dimension:
            break
    return all(
        sum((x * y for x, y in zip(row, final)), Fraction(0)) == 0
        for row in space.rows
    )


def standardize(values: tuple[int, ...]) -> tuple[int, ...]:
    order = {value: i + 1 for i, value in enumerate(sorted(values))}
    return tuple(order[value] for value in values)


def layered_permutation(parts: tuple[int, ...]) -> tuple[int, ...]:
    result: list[int] = []
    offset = 0
    for part in parts:
        result.extend(range(offset + part, offset, -1))
        offset += part
    return tuple(result)


def active_cells(r: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    partial = [0]
    for part in r:
        partial.append(partial[-1] + part)
    cells: set[tuple[int, int]] = set()
    for j, part in enumerate(r):
        for h in range(part + 1):
            cells.add((partial[j] + h, partial[j + 1] - h))
    cells.update((value, value) for value in partial)
    return frozenset(cells)


def signature_from_shading(
    r: tuple[int, ...], shading: frozenset[tuple[int, int]]
) -> Signature:
    partial = [0]
    for part in r:
        partial.append(partial[-1] + part)
    channels = [
        {
            (partial[j] + h, partial[j + 1] - h)
            for h in range(part + 1)
        }
        for j, part in enumerate(r)
    ]
    z = tuple(len(shading & channel) for channel in channels)
    E = frozenset(j for j, value in enumerate(partial) if (value, value) in shading)
    return Signature(r, z, E)


def mesh_occurrence_count(
    pattern: tuple[int, ...],
    shading: frozenset[tuple[int, int]],
    permutation: tuple[int, ...],
) -> int:
    n, k = len(permutation), len(pattern)
    count = 0
    for chosen in itertools.combinations(range(n), k):
        selected = tuple(permutation[i] for i in chosen)
        if standardize(selected) != pattern:
            continue
        positions = (0,) + tuple(i + 1 for i in chosen) + (n + 1,)
        values = (0,) + tuple(sorted(selected)) + (n + 1,)
        selected_positions = {i + 1 for i in chosen}
        valid = True
        for a, b in shading:
            for position, value in enumerate(permutation, start=1):
                if position in selected_positions:
                    continue
                if (
                    positions[a] < position < positions[a + 1]
                    and values[b] < value < values[b + 1]
                ):
                    valid = False
                    break
            if not valid:
                break
        count += int(valid)
    return count


def verify_direct_mesh_formula(maximum_n: int = 6) -> int:
    checked = 0
    for r in compositions(3):
        cells = sorted(active_cells(r))
        pattern = layered_permutation(r)
        for mask in range(1 << len(cells)):
            shading = frozenset(
                cells[i] for i in range(len(cells)) if mask & (1 << i)
            )
            signature = signature_from_shading(r, shading)
            for n in range(1, maximum_n + 1):
                for host in compositions(n):
                    direct = mesh_occurrence_count(
                        pattern, shading, layered_permutation(host)
                    )
                    assert direct == occurrence(signature, host)
            checked += 1
    return checked


def distribution_polynomial(signature: Signature, n: int) -> list[list[int]]:
    polynomial = Counter(occurrence(signature, host) for host in compositions(n))
    return [[exponent, polynomial[exponent]] for exponent in sorted(polynomial)]


def load_json(name: str):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def enrich_orbit_certificate(distributions_raw, pointwise_raw):
    pointwise_members = {
        item["pointwise_id"]: sorted(item["members"]) for item in pointwise_raw
    }
    enriched = []
    for original in distributions_raw:
        item = dict(original)
        item["representative_signature"] = min(item["members"])
        local_orbits = []
        for local_label, pointwise_ids in enumerate(item["reversal_orbits"], start=1):
            members = sorted(
                {
                    signature
                    for pointwise_id in pointwise_ids
                    for signature in pointwise_members[pointwise_id]
                }
            )
            local_orbits.append(
                {
                    "local_label": local_label,
                    "pointwise_ids": pointwise_ids,
                    "representative_signature": members[0],
                    "members": members,
                }
            )
        item["local_orbits"] = local_orbits
        enriched.append(item)
    return enriched


def build_lower_bound_certificate(distributions_raw):
    representatives = {
        item["distribution_id"]: parse_signature(min(item["members"]))
        for item in distributions_raw
    }
    polynomials = {
        class_id: {
            n: distribution_polynomial(signature, n) for n in range(8)
        }
        for class_id, signature in representatives.items()
    }
    separations = []
    histogram: Counter[int] = Counter()
    maximum = 0
    for left, right in itertools.combinations(sorted(representatives), 2):
        smallest_n = next(
            n for n in range(8) if polynomials[left][n] != polynomials[right][n]
        )
        histogram[smallest_n] += 1
        maximum = max(maximum, smallest_n)
        separations.append(
            {
                "left_class": left,
                "right_class": right,
                "smallest_n": smallest_n,
                "left_polynomial": polynomials[left][smallest_n],
                "right_polynomial": polynomials[right][smallest_n],
            }
        )
    return {
        "description": (
            "Exact first separating size and the two distribution polynomials "
            "for every pair of the 146 nonzero class representatives.  A "
            "polynomial is stored as [exponent, coefficient] pairs."
        ),
        "arithmetic": "exact integers",
        "class_representatives": {
            str(class_id): signature.text()
            for class_id, signature in representatives.items()
        },
        "pair_count": len(separations),
        "maximum_smallest_separating_size": maximum,
        "first_witness_histogram": {
            str(n): histogram[n] for n in sorted(histogram)
        },
        "separations": separations,
    }


def write_json(name: str, data) -> None:
    (ROOT / name).write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def verify_certificates(write_certificates: bool) -> None:
    signature_counts = [signature_count(k) for k in range(1, 11)]
    assert signature_counts[:4] == [12, 88, 644, 4712]
    assert all(
        signature_counts[k - 1]
        == 8 * signature_counts[k - 2] - 5 * signature_counts[k - 3]
        for k in range(3, 11)
    )
    print("arbitrary-length signature enumeration:", signature_counts[:4])

    signatures = all_signatures()
    assert len(signatures) == 644

    for (r, z), expansion in KERNEL_BASIS.items():
        for ell in range(1, 20):
            expanded = sum(
                coefficient * basis_value(letter, ell)
                for letter, coefficient in expansion.items()
            )
            assert expanded == kernel(r, z, ell)
    print("kernel basis expansions: OK")

    pointwise_raw = load_json("pointwise_groups.json")
    distributions_raw = load_json("equidistribution_audit.json")
    pointwise_groups = {
        item["pointwise_id"]: [parse_signature(x) for x in item["members"]]
        for item in pointwise_raw
    }
    distributions = {
        item["distribution_id"]: [parse_signature(x) for x in item["members"]]
        for item in distributions_raw
    }
    assert len(pointwise_groups) == 523
    assert len(distributions) == 146
    assert {
        signature for members in pointwise_groups.values() for signature in members
    } == set(signatures)
    assert {
        signature for members in distributions.values() for signature in members
    } == set(signatures)

    for members in pointwise_groups.values():
        representative = members[0]
        for signature in members[1:]:
            assert series_equal(representative, signature)
    representative_keys = {
        pointwise_id: tensor_key(members[0], 4)
        for pointwise_id, members in pointwise_groups.items()
    }
    assert len(set(representative_keys.values())) == 523
    progression = [
        len({tensor_key(signature, m) for signature in signatures})
        for m in range(1, 5)
    ]
    assert progression == [6, 45, 223, 523]
    print("all-layer pointwise comparison: 644 -> 523; finite check", progression)

    signature_to_pointwise = {
        signature: pointwise_id
        for pointwise_id, members in pointwise_groups.items()
        for signature in members
    }
    fixed, rc_orbits = 0, set()
    for pointwise_id, members in pointwise_groups.items():
        rc_ids = {signature_to_pointwise[signature.rc()] for signature in members}
        assert len(rc_ids) == 1
        rc_id = next(iter(rc_ids))
        assert {signature.rc() for signature in members} == set(
            pointwise_groups[rc_id]
        )
        fixed += int(rc_id == pointwise_id)
        rc_orbits.add(tuple(sorted((pointwise_id, rc_id))))
    assert fixed == 49 and len(rc_orbits) == 286

    unique_support = [signature for signature in signatures if signature.f <= 1]
    assert Counter(signature.f for signature in unique_support) == Counter({0: 56, 1: 190})
    assert len(
        {(signature.f, signature.degree_multiset) for signature in unique_support}
    ) == 28
    print("quotient checks: 523 -> 286; unique-support families -> 28")

    enriched = enrich_orbit_certificate(distributions_raw, pointwise_raw)
    coordinate_witnesses = build_coordinate_witness_certificate()
    verify_coordinate_witnesses(coordinate_witnesses, enriched)
    print("coordinate witnesses: 43 explicit identities checked through m=30")

    lower_bound = build_lower_bound_certificate(enriched)
    assert lower_bound["pair_count"] == math.comb(146, 2)
    assert lower_bound["maximum_smallest_separating_size"] == 7

    if write_certificates:
        write_json("equidistribution_audit.json", enriched)
        write_json("coordinate_witnesses.json", coordinate_witnesses)
        write_json("lower_bound_separation.json", lower_bound)
        print("certificates regenerated")
    else:
        committed_coordinate = load_json("coordinate_witnesses.json")
        committed_lower = load_json("lower_bound_separation.json")
        assert committed_coordinate == coordinate_witnesses
        assert committed_lower == lower_bound
        assert distributions_raw == enriched
        print("orbit, coordinate, and lower-bound certificates: exact match")

    checked = verify_direct_mesh_formula(6)
    assert checked == 1600
    print("direct mesh/signature check: 1,600 effective shadings through n=6")
    print("REPRODUCIBILITY AUDIT COMPLETE")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-certificates",
        action="store_true",
        help=(
            "regenerate equidistribution_audit.json, coordinate_witnesses.json, "
            "and lower_bound_separation.json"
        ),
    )
    args = parser.parse_args()
    verify_certificates(args.write_certificates)


if __name__ == "__main__":
    main()
