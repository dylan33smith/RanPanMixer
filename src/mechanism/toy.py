"""The toy cohort: small enough that every output path can be enumerated exactly.

This module is the falsification target for Theorem 1, and it is the reason the
sampler is written against something checkable rather than against chr21. At
scale `tv_empirical` is biased UPWARD by the plug-in estimator, so a large
measured value there is the estimator misbehaving and a small one proves nothing.
Only exact enumeration on a support this size can falsify the theorem.

⚠ IT READS NOTHING. No graph, no VCF, no cohort artifact, no PanMixer path. The
allele matrix below is hand-written, so a defect in preprocessing cannot mask a
defect in the mechanism, and the gate keeps working on a machine with no `data/`.

WHAT EACH FEATURE OF THE TOY IS FOR. Every row and column earns its place; a toy
that exercises only the easy path is worse than no toy, because it reports green.

| feature                         | position | what it makes testable                     |
|---------------------------------|----------|--------------------------------------------|
| D1 and D2 differ EVERYWHERE     | all      | the sharp pair -- see `sharp_pair()`       |
| D3 == D1 except one cell        | v4       | two state paths, ONE output path           |
| a missing donor allele          | v4, v5   | `missing_policy` = RENORMALISE             |
| a three-allele position         | v2       | `n_alleles` bounds, non-binary support     |
| a ZERO genetic-distance step    | v2 -> v3 | P(switch) = 0, so Y_D is strictly smaller  |
| a very LARGE step               | v5 -> v6 | switching is near-uniform                  |
| unequal support                 | v4, v5   | non-uniform beta_v = w_v / sum w_u         |

THE ZERO-DISTANCE STEP IS THE LOAD-BEARING ONE. chr21 has 618 interior map
intervals with dcM = 0 holding about 1,650 kept variants, so this is real data,
not a contrivance. It has two consequences the toy exists to pin down: the output
support `Y_D` genuinely shrinks (paths that switch donors across it have
probability zero, and the mechanism must never give them mass at any `tau`), and
it is what rules out reading `missing_policy` as a locally renormalised
transition -- across a forced-stay step, a donor that becomes unavailable has a
transition row summing to zero, so there is no chain to renormalise. See
`hmm.py`.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

# Hand-written. Rows are donor haplotypes, columns are chain positions, -1 is
# "no called allele" -- the same encoding as the real `haplotypes.npy`, so the
# mechanism cannot tell a toy cohort from a chromosome.
_HAPLOTYPES = np.array(
    [[0, 1, 2, 0, 1, 0, 1],    # D1  fully called
     [1, 0, 0, 1, 0, 1, 0],    # D2  fully called, differs from D1 at EVERY position
     [0, 1, 2, 0, -1, 0, 1],   # D3  == D1 except no call at v4
     [1, 0, 1, 1, 0, -1, 0]],  # D4  == D2 except v2, and no call at v5
    dtype=np.int16,
)

# Base-pair positions. Strictly increasing, as VCF record order requires.
_POSITIONS = np.array([1_000, 1_200, 2_200, 3_200, 3_600, 200_000, 200_200],
                      dtype=np.int64)

# centiMorgans, as a genetic map would supply them -- NOT derived from the bp
# positions. v2 -> v3 is a plateau (dcM = 0 across 1,000 bp), which is what the
# chr21 map does in 618 places; v5 -> v6 is a 0.02 cM gap, far enough that
# switching goes uniform.
_GENETIC_POS = np.array([0.0, 0.00002, 0.00012, 0.00012, 0.00016, 0.02, 0.02002],
                        dtype=np.float64)

# Declared alleles per position, as `num_alleles` = n_alt + 1 from the VCF. Held
# separately from the matrix because the real thing is: the graph declares the
# alleles, and which of them the cohort happens to carry is a different question.
_N_ALLELES = np.array([2, 2, 3, 2, 2, 2, 2], dtype=np.int32)


@dataclass(frozen=True)
class Toy:
    """A cohort the mechanism cannot distinguish from a real one.

    Field names and dtypes deliberately match the cohort artifacts -- the same
    code path consumes both, so there is no test mode to drift out of sync with
    the real one. That is the same reason `--root` exists in `cli.py`.
    """

    haplotypes: np.ndarray    # (K, T) int16, -1 = no called allele
    positions: np.ndarray     # (T,) int64   bp
    genetic_pos: np.ndarray   # (T,) float64 cM
    n_alleles: np.ndarray     # (T,) int32   declared allele count per position
    name: str = "toy4x7"

    @property
    def K(self) -> int:
        return self.haplotypes.shape[0]

    @property
    def T(self) -> int:
        return self.haplotypes.shape[1]

    @property
    def support(self) -> np.ndarray:
        """support(v): called haplotypes per position. The `w_v` denominator."""
        return (self.haplotypes != -1).sum(axis=0).astype(np.int32)

    @property
    def allowed(self) -> np.ndarray:
        """(K, T) bool: which donors the chain may occupy, per position.

        This is `missing_policy` = RENORMALISE expressed as data. A donor with no
        called allele has no `phi_v`, so it is not in the state set there.
        """
        return self.haplotypes != -1


def toy_cohort() -> Toy:
    """The toy, with its invariants checked at construction.

    The checks are here rather than in a test because a toy that has quietly
    stopped exercising a feature still reports green, and a green gate on a
    degenerate fixture is the failure mode this whole file exists to avoid.
    """
    t = Toy(haplotypes=_HAPLOTYPES.copy(), positions=_POSITIONS.copy(),
            genetic_pos=_GENETIC_POS.copy(), n_alleles=_N_ALLELES.copy())

    if not np.all(np.diff(t.positions) > 0):
        raise ValueError("toy positions are not strictly increasing")
    if not np.all(np.diff(t.genetic_pos) >= 0):
        raise ValueError("toy cM is not monotonic non-decreasing")
    called = t.haplotypes[t.allowed]
    if called.size and called.max() >= t.n_alleles.max():
        bad = np.flatnonzero((t.haplotypes >= t.n_alleles[None, :]).any(axis=0))
        raise ValueError(f"observed allele index >= n_alleles at positions {bad}")
    if (t.support == 0).any():
        raise ValueError("a toy position has support 0; the real chain excludes those")

    # The features the table at the top of this file claims. Each is a real gate.
    if not (t.haplotypes[0] != t.haplotypes[1]).all():
        raise ValueError("D1 and D2 no longer differ everywhere; sharp_pair() is blunted")
    if not (np.diff(t.genetic_pos) == 0).any():
        raise ValueError("the zero-distance step is gone; forced-stay is untested")
    if not (t.haplotypes == -1).any():
        raise ValueError("no missing donor allele; RENORMALISE is untested")
    if len(set(t.support.tolist())) < 2:
        raise ValueError("support is uniform; non-uniform beta_v is untested")
    if t.n_alleles.max() < 3:
        raise ValueError("no multi-allelic position; allele bounds are untested")
    return t


def sharp_pair(toy: Toy) -> tuple[np.ndarray, np.ndarray]:
    """The input pair that SATURATES the pointwise likelihood-ratio bound.

    Why this pair exists and why it matters. The theorem's pointwise statement is
    that for any two inputs p, q the ratio `Q_p(y)/Q_q(y)` spans a factor of at
    most `exp(epsilon_tau)` across the support, because

        Q_p(y)/Q_q(y) = exp(eta_tau * (u(p,y) - u(q,y))) * Z_q/Z_p

    and the `Z` factor is common to every `y`, so it cancels out of the SPREAD.
    The spread is therefore exactly `exp(eta_tau * range_y(u_p - u_q))`, and the
    range is at most 2 because each utility is in [0,1].

    Take p = D1's alleles and q = D2's alleles. They differ at every position, and
    both "copy D1 throughout" and "copy D2 throughout" are in the support. So
    `u_p - u_q` is +1 at one output and -1 at another, the range is exactly 2, and
    the spread EQUALS `exp(epsilon_tau)` -- the bound is attained, not approached.

    ⚠ THIS IS WHY THE PAIR IS CONSTRUCTED AND NOT SAMPLED. A saturating pair makes
    the gate sharp: any inflation of the tilt exponent fails it immediately. The
    `tv_empirical <= tau` gate does NOT have that property -- measured on this
    toy, the worst exact TV reaches only about half of `tau`, so TV alone cannot
    detect a doubled `eta_tau`, which is the single most likely calibration error
    in the project (`epsilon_tau = 2 * eta_tau`, and reaching for the wrong one is
    a one-character mistake). Keep both gates; only this one is sharp.
    """
    return toy.haplotypes[0].copy(), toy.haplotypes[1].copy()


def toy_inputs(toy: Toy, *, kind: str = "called") -> np.ndarray:
    """Every admissible input path, enumerated.

    `kind="called"` gives the fully-called inputs: the target has a declared
    allele at every position. `kind="all"` also admits -1, which is what a real
    `path.npy` carries wherever the join found no record for the target -- 14.9%
    of chr21 chain positions. Those inputs are admissible and must be covered,
    because their maximum achievable utility is below 1 and that is exactly the
    regime the deferred decision R1 is about.
    """
    if kind == "called":
        per_pos = [range(int(n)) for n in toy.n_alleles]
    elif kind == "all":
        per_pos = [[-1, *range(int(n))] for n in toy.n_alleles]
    else:
        raise ValueError(f"kind must be 'called' or 'all', not {kind!r}")
    return np.array(list(itertools.product(*per_pos)), dtype=np.int16)
