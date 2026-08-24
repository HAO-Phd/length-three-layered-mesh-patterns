from __future__ import annotations

import itertools
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
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
    def D(self) -> tuple[int, ...]:
        return tuple(sorted(r - z for r, z in zip(self.r, self.z) if z <= r))

    def rc(self) -> "Signature":
        s = self.s
        return Signature(self.r[::-1], self.z[::-1], frozenset(s - e for e in self.E))

    def text(self) -> str:
        rs = ",".join(map(str, self.r))
        zs = ",".join(map(str, self.z))
        es = ",".join(map(str, sorted(self.E)))
        return f"r=({rs});z=({zs});E={{{es}}}"


def compositions(n: int):
    if n == 0:
        yield ()
        return
    for mask in range(1 << (n - 1)):
        out = []
        last = 0
        for i in range(n - 1):
            if mask & (1 << i):
                out.append(i + 1 - last)
                last = i + 1
        out.append(n - last)
        yield tuple(out)


def all_signatures() -> list[Signature]:
    result = []
    for r in compositions(3):
        s = len(r)
        for z in itertools.product(*(range(x + 2) for x in r)):
            for mask in range(1 << (s + 1)):
                E = frozenset(i for i in range(s + 1) if mask & (1 << i))
                result.append(Signature(r, tuple(z), E))
    return sorted(result)


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


def basis_value(i: int, x: int) -> int:
    if i == 0:
        return 1
    if i == 1:
        return x
    if i == 2:
        return math.comb(x, 2)
    if i == 3:
        return math.comb(x, 3)
    return int(x == i - 3)


def supports(sig: Signature, m: int):
    s = sig.s
    for phi in itertools.combinations(range(m), s):
        if 0 in sig.E and phi[0] != 0:
            continue
        if s in sig.E and phi[-1] != m - 1:
            continue
        if any(j in sig.E and phi[j] != phi[j - 1] + 1 for j in range(1, s)):
            continue
        yield phi


def occurrence(sig: Signature, L: tuple[int, ...]) -> int:
    total = 0
    for phi in supports(sig, len(L)):
        weight = 1
        for j, pos in enumerate(phi):
            weight *= kernel(sig.r[j], sig.z[j], L[pos])
        total += weight
    return total


def tensor(sig: Signature, m: int) -> dict[tuple[int, ...], int]:
    answer: dict[tuple[int, ...], int] = defaultdict(int)
    for phi in supports(sig, m):
        choices = [list(KERNEL_BASIS[(sig.r[j], sig.z[j])].items()) for j in range(sig.s)]
        for selected in itertools.product(*choices):
            word = [0] * m
            coeff = 1
            for pos, (letter, c) in zip(phi, selected):
                word[pos] = letter
                coeff *= c
            answer[tuple(word)] += coeff
    return {w: c for w, c in answer.items() if c}


def tensor_key(sig: Signature, max_m: int) -> tuple:
    return tuple(tuple(sorted(tensor(sig, m).items())) for m in range(1, max_m + 1))


def parse_signature(text: str) -> Signature:
    match = re.fullmatch(r"r=\(([^)]*)\);z=\(([^)]*)\);E=\{([^}]*)\}", text)
    if not match:
        raise ValueError(text)
    r = tuple(int(x) for x in match.group(1).split(",") if x)
    z = tuple(int(x) for x in match.group(2).split(",") if x)
    E = frozenset(int(x) for x in match.group(3).split(",") if x)
    return Signature(r, z, E)


def read_certificates():
    pointwise_raw = json.loads((ROOT / "pointwise_groups.json").read_text(encoding="utf-8"))
    dist_raw = json.loads((ROOT / "equidistribution_audit.json").read_text(encoding="utf-8"))
    pointwise = {
        item["pointwise_id"]: frozenset(parse_signature(x) for x in item["members"])
        for item in pointwise_raw
    }
    distributions = {
        item["distribution_id"]: frozenset(parse_signature(x) for x in item["members"])
        for item in dist_raw
    }
    return pointwise_raw, dist_raw, pointwise, distributions


def standardize(values: tuple[int, ...]) -> tuple[int, ...]:
    order = {v: i + 1 for i, v in enumerate(sorted(values))}
    return tuple(order[v] for v in values)


def layered_permutation(L: tuple[int, ...]) -> tuple[int, ...]:
    result = []
    offset = 0
    for ell in L:
        result.extend(range(offset + ell, offset, -1))
        offset += ell
    return tuple(result)


def underlying(r: tuple[int, ...]) -> tuple[int, ...]:
    return layered_permutation(r)


def active_cells(r: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    result = set()
    R = [0]
    for x in r:
        R.append(R[-1] + x)
    for j, x in enumerate(r):
        for h in range(x + 1):
            result.add((R[j] + h, R[j + 1] - h))
    result.update((x, x) for x in R)
    return frozenset(result)


def signature_from_shading(r: tuple[int, ...], shading: frozenset[tuple[int, int]]) -> Signature:
    R = [0]
    for x in r:
        R.append(R[-1] + x)
    channels = []
    for j, x in enumerate(r):
        channels.append({(R[j] + h, R[j + 1] - h) for h in range(x + 1)})
    z = tuple(len(shading & channel) for channel in channels)
    E = frozenset(j for j, x in enumerate(R) if (x, x) in shading)
    return Signature(r, z, E)


def mesh_occurrence_count(
    tau: tuple[int, ...], shading: frozenset[tuple[int, int]], pi: tuple[int, ...]
) -> int:
    n = len(pi)
    k = len(tau)
    count = 0
    for inds0 in itertools.combinations(range(n), k):
        selected = tuple(pi[i] for i in inds0)
        if standardize(selected) != tau:
            continue
        inds = (0,) + tuple(i + 1 for i in inds0) + (n + 1,)
        vals = (0,) + tuple(sorted(selected)) + (n + 1,)
        selected_positions = set(i + 1 for i in inds0)
        valid = True
        for a, b in shading:
            for j, value in enumerate(pi, 1):
                if j in selected_positions:
                    continue
                if inds[a] < j < inds[a + 1] and vals[b] < value < vals[b + 1]:
                    valid = False
                    break
            if not valid:
                break
        count += valid
    return count


def verify_active_and_signature_formula(max_n: int = 6):
    recovered = defaultdict(set)
    for r in compositions(3):
        tau = underlying(r)
        for L in compositions(4):
            pi = layered_permutation(L)
            for inds0 in itertools.combinations(range(4), 3):
                selected = tuple(pi[i] for i in inds0)
                if standardize(selected) != tau:
                    continue
                unused = next(i for i in range(4) if i not in inds0)
                inds = (0,) + tuple(i + 1 for i in inds0) + (5,)
                vals = (0,) + tuple(sorted(selected)) + (5,)
                j = unused + 1
                value = pi[unused]
                a = next(a for a in range(4) if inds[a] < j < inds[a + 1])
                b = next(b for b in range(4) if vals[b] < value < vals[b + 1])
                recovered[r].add((a, b))
        assert recovered[r] == set(active_cells(r)), (r, recovered[r], active_cells(r))

    checked = 0
    for r in compositions(3):
        cells = sorted(active_cells(r))
        tau = underlying(r)
        for mask in range(1 << len(cells)):
            shading = frozenset(cells[i] for i in range(len(cells)) if mask & (1 << i))
            sig = signature_from_shading(r, shading)
            for n in range(1, max_n + 1):
                for L in compositions(n):
                    pi = layered_permutation(L)
                    direct = mesh_occurrence_count(tau, shading, pi)
                    via_signature = occurrence(sig, L)
                    assert direct == via_signature, (r, shading, sig, L, direct, via_signature)
            checked += 1
    return checked


def distributions_through(signatures: list[Signature], max_n: int):
    result = {sig: [] for sig in signatures}
    for n in range(max_n + 1):
        comps = list(compositions(n))
        for sig in signatures:
            result[sig].append(tuple(sorted(Counter(occurrence(sig, L) for L in comps).items())))
    return result


def verify_lower_bound(distribution_rows: dict[int, frozenset[Signature]], max_n: int = 12):
    representatives = {row: min(members) for row, members in distribution_rows.items()}
    series = distributions_through(list(representatives.values()), max_n)
    first_witness = {}
    for a, b in itertools.combinations(sorted(representatives), 2):
        sa, sb = representatives[a], representatives[b]
        witness = next((n for n in range(max_n + 1) if series[sa][n] != series[sb][n]), None)
        assert witness is not None, (a, b)
        first_witness[a, b] = witness
    return max(first_witness.values()), Counter(first_witness.values())


def main():
    signatures = all_signatures()
    assert len(signatures) == 644
    assert len(set(signatures)) == 644

    for (r, z), expansion in KERNEL_BASIS.items():
        for ell in range(1, 20):
            value = sum(c * basis_value(i, ell) for i, c in expansion.items())
            assert value == kernel(r, z, ell), ((r, z), ell, value, kernel(r, z, ell))
    print("kernel expansions: OK")

    pointwise_raw, dist_raw, pointwise_cert, dist_cert = read_certificates()
    assert len(pointwise_cert) == 523
    assert frozenset().union(*pointwise_cert.values()) == frozenset(signatures)
    assert len(dist_cert) == 146
    assert frozenset().union(*dist_cert.values()) == frozenset(signatures)
    print("certificate coverage: 644 signatures, 523 pointwise groups, 146 rows")

    current = {sig: () for sig in signatures}
    progression = []
    for m in range(1, 16):
        current = {sig: current[sig] + (tuple(sorted(tensor(sig, m).items())),) for sig in signatures}
        progression.append(len(set(current.values())))
    computed_groups = defaultdict(set)
    for sig, key in current.items():
        computed_groups[key].add(sig)
    computed_partition = {frozenset(group) for group in computed_groups.values()}
    certificate_partition = set(pointwise_cert.values())
    assert computed_partition == certificate_partition
    print("pointwise tensor progression m=1..15:", progression)
    print("pointwise certificate: exact match")

    fixed = 0
    seen = set()
    sig_to_pid = {sig: pid for pid, members in pointwise_cert.items() for sig in members}
    for pid, members in pointwise_cert.items():
        rc_members = frozenset(sig.rc() for sig in members)
        rc_pid = sig_to_pid[next(iter(rc_members))]
        assert pointwise_cert[rc_pid] == rc_members
        if rc_pid == pid:
            fixed += 1
        seen.add(tuple(sorted((pid, rc_pid))))
    assert fixed == 49 and len(seen) == 286
    print("reverse-complement quotient: 49 fixed, 286 orbits")

    unique = [sig for sig in signatures if sig.f <= 1]
    assert Counter(sig.f for sig in unique) == Counter({0: 56, 1: 190})
    assert len({(sig.f, sig.D) for sig in unique}) == 28
    print("unique-support reduction: 56 + 190 signatures -> 28 classes")

    checked = verify_active_and_signature_formula(6)
    assert checked == 1600
    print("active cells and direct mesh/signature comparison: 1600 shadings through n=6 OK")

    max_witness, witness_hist = verify_lower_bound(dist_cert, 12)
    print("lower-bound separation: all 146 representatives separated; maximum first n =", max_witness)
    print("lower-bound first-witness histogram:", dict(sorted(witness_hist.items())))

    ancillary_keys = set().union(*(item.keys() for item in dist_raw))
    print("equidistribution_audit.json keys:", sorted(ancillary_keys))
    print("AUDIT CORE COMPLETE")


if __name__ == "__main__":
    main()
