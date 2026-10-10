"""Brute force: every output path, every probability, exactly. The falsifier.

This module is the reference the mechanism is checked against, and it is written
to be WRONG IN DIFFERENT WAYS than the mechanism is. It re-derives the transition
probability from the published closed form and materialises the dense K x K
matrix; `hmm.py` never materialises anything and uses the rank-one decomposition.
It enumerates all K^T state paths with a Python loop; `ffbs.py` does one O(TK)
forward pass. Neither imports the other.

That independence is the whole value. `support.py` says it already about a
cheaper check: PanMixer's equivalent assertion compares two expressions that are
equal by construction, so it can never fire. Two agreeing implementations are
evidence only when they could have disagreed.

WHAT THIS CAN FALSIFY, with no sampler in existence:

1. `log_z_p` in [0, eta_tau] -- computed here by direct summation over the
   support, not from a forward pass, so it tests the DEFINITION of the release
   distribution rather than an implementation of it.
2. `TV(Q_p, Q_q) <= tau` for every pair of admissible inputs. Exact, so a
   violation is a violation.
3. The pointwise ratio spread `<= exp(epsilon_tau)`, which is sharp -- see
   `toy.sharp_pair`. This is the gate that catches a miscalibrated tilt.
4. `tau = 0` reproduces `baseline_model` exactly.

Points 1-4 need only the definitions of `R_D` and `u_path`. So the toy gate is
buildable BEFORE the sampler, and the sampler is then written against a target
that is already known to be self-consistent. Checking that `ffbs` draws from
`Q_p` is a fifth, separate gate that does need the sampler.

⚠ COST. K^T state paths, so this is for the toy and nothing else. At K = 88,
T = 305,886 the enumeration has 10^576,000 terms; the whole point of `ffbs` is
that the same distribution is sampled exactly in 30M operations.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

# The published Methods, `paper/s41467-026-77591-0_reference.pdf` PDF page 10:
#   Ne = 10,000, r = 1.26, Delta_x in centiMorgans, d = Delta_x * Ne * r,
#   p_stay = e^(-d/n) + (1 - e^(-d/n)) / n.
# Stated again here rather than imported from `hmm.py` so that the two modules
# can disagree about it. ⚠ NOT PanMixer's constant: its EFFECTIVE_N is 1/10,000,
# the reciprocal, and its "distance" is a difference of VCF row indices.
NE = 10_000.0
R = 1.26


def dense_transition(delta_cm: float, n_states: int) -> np.ndarray:
    """The full n x n Li-Stephens transition matrix, written from the closed form.

    `delta_cm` is the genetic distance to the next position, in centiMorgans.
    Off-diagonal entries are the per-state switch probability and the diagonal is
    the published `p_stay`; the two are filled independently and the row sums are
    checked, so a factor-of-n error in either cannot pass silently. That error is
    the one `stay_switch` in `docs/terms.md` warns about by name.
    """
    if delta_cm < 0:
        raise ValueError(f"negative genetic distance {delta_cm}; the map is not monotonic")
    d = delta_cm * NE * R
    decay = math.exp(-d / n_states)
    per_state_switch = (1.0 - decay) / n_states
    p_stay = decay + per_state_switch
    A = np.full((n_states, n_states), per_state_switch, dtype=np.float64)
    np.fill_diagonal(A, p_stay)
    rows = A.sum(axis=1)
    if not np.allclose(rows, 1.0, atol=1e-12):
        raise AssertionError(f"transition rows sum to {rows} at delta_cm={delta_cm}")
    return A


@dataclass(frozen=True)
class ExactBaseline:
    """`R_D` enumerated: the output support and its probability, exhaustively.

    `paths` are OUTPUT paths -- allele vectors, not donor-state vectors. That
    distinction is not pedantry. The released object is a path of alleles, and
    two different donor-state sequences can emit the same one (the toy's D3 is D1
    with one cell removed, precisely so this is exercised). Theorem 1 is a
    statement about the distribution over what is RELEASED, so the collapse has
    to happen before any distance is measured; measuring over state paths instead
    would answer a question nobody asked and would be quietly conservative.
    """

    paths: np.ndarray        # (M, T) int16  -- Y_D, the output support
    prob: np.ndarray         # (M,) float64  -- R_D(y), sums to 1
    n_state_paths: int       # how many donor-state sequences had positive mass
    n_enumerated: int        # K^T, for the record

    @property
    def M(self) -> int:
        return len(self.prob)


def enumerate_baseline(haplotypes: np.ndarray, genetic_pos: np.ndarray) -> ExactBaseline:
    """Enumerate `R_D` over output paths, from the cohort alone.

    ⚠ NO TARGET ARGUMENT, and that is structural, not stylistic. `R_D` and the
    output support `Y_D` are functions of `G` and `D` only (standing constraint
    1); a target-dependent support is the §4.4 disjoint-support argument that
    voids the theorem. The support this returns is the one every `tau` shares,
    and the toy gate asserts it does not move when the target does.

    RENORMALISE is applied as an INDICATOR over states, not as a renormalised
    transition row -- see `hmm.py` for why the locally-renormalised reading is
    not even well defined. So `R_D` here is the Li-Stephens chain CONDITIONED on
    never occupying a donor with no called allele, and the conditioning constant
    is divided out below.
    """
    K, T = haplotypes.shape
    allowed = haplotypes != -1
    if not allowed.any(axis=0).all():
        dead = np.flatnonzero(~allowed.any(axis=0))
        raise ValueError(
            f"positions {dead} have no called donor at all, so the chain has "
            f"nothing to emit. The real chain excludes support == 0 positions "
            f"upstream, in genmap.chain_mask -- this cohort should not reach here."
        )

    A = [dense_transition(float(genetic_pos[t + 1] - genetic_pos[t]), K)
         for t in range(T - 1)]
    rho = np.full(K, 1.0 / K)     # uniform over donors: a cohort-only choice

    acc: dict[bytes, float] = {}
    key_to_path: dict[bytes, np.ndarray] = {}
    n_state_paths = 0
    for z in itertools.product(range(K), repeat=T):
        if not all(allowed[z[t], t] for t in range(T)):
            continue                                  # RENORMALISE: state excluded
        p = rho[z[0]]
        for t in range(1, T):
            p *= A[t - 1][z[t - 1], z[t]]
        if p <= 0.0:
            continue                                  # e.g. a forced-stay step
        n_state_paths += 1
        y = np.array([haplotypes[z[t], t] for t in range(T)], dtype=np.int16)
        k = y.tobytes()
        acc[k] = acc.get(k, 0.0) + p
        key_to_path.setdefault(k, y)

    if not acc:
        raise ValueError("R_D has empty support: every path was excluded")

    keys = sorted(acc)
    prob = np.array([acc[k] for k in keys], dtype=np.float64)
    prob /= prob.sum()                                # the conditioning constant
    paths = np.stack([key_to_path[k] for k in keys])
    return ExactBaseline(paths=paths, prob=prob, n_state_paths=n_state_paths,
                         n_enumerated=K ** T)


def utility_vector(base: ExactBaseline, target: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """`u(p, y)` for every `y` in the support, as a vector over the support.

    `phi_v` is an indicator under model C: 1 if the output carries the target's
    allele at `v`, else 0. A position where the TARGET has no call (-1) scores 0
    for every donor -- the potential is 1 everywhere there, so the position
    contributes no tilt at all. That is the inert budget of deferred decision R1,
    and it is why `u` cannot reach 1 for such an input.
    """
    if beta.shape != target.shape:
        raise ValueError(f"beta has shape {beta.shape}, target {target.shape}")
    match = (base.paths == target[None, :]) & (target[None, :] != -1)
    return match @ beta


def release_distribution(base: ExactBaseline, u: np.ndarray, eta: float
                         ) -> tuple[np.ndarray, float]:
    """`Q^tau_p` over the support, and `log Z_p`, by direct summation.

    Returns `(Q, log_z_p)` with `log_z_p = log sum_y R_D(y) exp(eta * u(p,y))`.
    Theory pins it to [0, eta_tau] for every input: `u` in [0,1] gives
    `Z_p` in [1, e^eta]. This is the cheapest correctness statement in the whole
    project and it is computed here WITHOUT a forward pass, so it checks the
    definition rather than an implementation of it.
    """
    if eta < 0:
        raise ValueError(f"eta_tau must be non-negative, got {eta}")
    if u.min() < -1e-12 or u.max() > 1 + 1e-12:
        raise AssertionError(
            f"u_path left [0,1]: [{u.min():.6g}, {u.max():.6g}]. An unbounded "
            f"utility does not weaken the bound -- it removes it (constraint 2)."
        )
    w = base.prob * np.exp(eta * u)
    z = w.sum()
    return w / z, float(np.log(z))


def tv(q1: np.ndarray, q2: np.ndarray) -> float:
    """Total variation distance, exact on a shared enumerated support."""
    return 0.5 * float(np.abs(q1 - q2).sum())


def ratio_spread(q1: np.ndarray, q2: np.ndarray, *, floor: float = 0.0) -> float:
    """`max_y Q_p(y)/Q_q(y)` divided by `min_y` of the same ratio.

    THE SHARP GATE. Theory: the ratio is
    `exp(eta_tau (u_p(y) - u_q(y))) * Z_q/Z_p`, and the `Z` factor is the same for
    every `y`, so the spread is `exp(eta_tau * range_y(u_p - u_q)) <=
    exp(2 * eta_tau) = exp(epsilon_tau)`. Nothing is estimated and nothing
    cancels approximately.

    It is strictly more powerful than `TV <= tau` on this toy, which is worth
    stating plainly: the worst exact TV here reaches only **27.4% of `tau` at
    `tau` = 0.01, rising to 43.1% at `tau` = 0.9** (measured), so a tilt applied
    with `epsilon_tau` in place of `eta_tau` -- twice the intended exponent, and a
    one-token error -- passes the TV gate at every `tau` we test and fails this one
    by a factor of `exp(epsilon_tau)`.

    Support entries below `floor` are dropped; both distributions share the
    enumerated support of `R_D`, so with `floor = 0` nothing is dropped and the
    statement covers all of `Y_D`.
    """
    keep = (q1 > floor) & (q2 > floor)
    if not keep.any():
        raise ValueError("no shared support left to compare")
    r = q1[keep] / q2[keep]
    return float(r.max() / r.min())


def dense_forward(haplotypes: np.ndarray, genetic_pos: np.ndarray,
                  potentials: np.ndarray) -> tuple[np.ndarray, float]:
    """An O(TK^2) forward pass with materialised matrices. The reference for `ffbs`.

    `potentials` is (T, K) and must already carry BOTH the RENORMALISE indicator
    and the tilt, because that is what the fast path multiplies in. Returns the
    per-position normalised alphas and `sum_t log c_t` -- which is the log
    normaliser of the potentials-and-transitions product, NOT `log_z_p`. Turning
    it into `log_z_p` means subtracting the baseline's own conditioning constant;
    `ffbs.py` does that and explains why forgetting to is a trap.

    Exists only so the rank-one decomposition can be checked against a dense
    matrix-vector product. `docs/plan.md` asks for exactly this equivalence.
    """
    K, T = haplotypes.shape
    if potentials.shape != (T, K):
        raise ValueError(f"potentials must be (T, K) = {(T, K)}, got {potentials.shape}")
    alpha = np.zeros((T, K), dtype=np.float64)
    logc = 0.0
    r = np.full(K, 1.0 / K) * potentials[0]
    c = r.sum()
    if c <= 0:
        raise ValueError("position 0 has zero total mass under its potentials")
    alpha[0] = r / c
    logc += math.log(c)
    for t in range(1, T):
        A = dense_transition(float(genetic_pos[t] - genetic_pos[t - 1]), K)
        r = (alpha[t - 1] @ A) * potentials[t]
        c = r.sum()
        if c <= 0:
            raise ValueError(f"position {t} has zero total mass; the chain dead-ends")
        alpha[t] = r / c
        logc += math.log(c)
    return alpha, logc


def enumerate_state_posterior(haplotypes: np.ndarray, genetic_pos: np.ndarray,
                              potentials: np.ndarray) -> np.ndarray:
    """The exact distribution over STATE paths under the given potentials.

    Returned as a (n_paths, T) index array's worth of probability, flattened to a
    dict-free dense vector in `itertools.product` order with zeros for excluded
    paths. Used to check that `ffbs`'s empirical state-path frequencies match the
    posterior it claims to sample -- a stronger check than comparing the output
    marginals, which can agree while the joint is wrong.
    """
    K, T = haplotypes.shape
    A = [dense_transition(float(genetic_pos[t + 1] - genetic_pos[t]), K)
         for t in range(T - 1)]
    out = np.zeros(K ** T, dtype=np.float64)
    for i, z in enumerate(itertools.product(range(K), repeat=T)):
        p = (1.0 / K) * potentials[0, z[0]]
        for t in range(1, T):
            p *= A[t - 1][z[t - 1], z[t]] * potentials[t, z[t]]
        out[i] = p
    s = out.sum()
    if s <= 0:
        raise ValueError("the potentials annihilate every state path")
    return out / s


# ---------------------------------------------------------------------------
# COMPOSITE RELEASES: diploid, multi-segment, and both at once
#
# ⚠ WHY THIS EXISTS. Until 2026-10-09 every falsification clause in this project
# was HAPLOID and single-chromosome, while the diploid path was the only one ever
# run -- all 18 releases of the first chr21 run carry `ploidy: 2`. There was no
# diploid primitive here at all, so the gate could not be applied to the code
# path that actually executes. That is missing REFERENCE CODE, not a missing
# assertion, and it is what this block supplies.
#
# THE ALGEBRA, because getting it wrong is the error the gate exists to catch.
# A release is made of INDEPENDENT chains -- one per (haplotype, chromosome)
# pair. The utility is the MEAN over haplotypes of a SUM over chromosomes:
#
#     u(p, y) = (1 / ploidy) * sum_h sum_c u_{h,c}
#
# where `u_{h,c}` is that factor's own weighted indicator and its `beta` sums to
# that chromosome's SHARE `w_c` of the one global budget (standing constraint 4,
# with sum_c w_c = 1). Then
#
#     Q(y) prop. prod_{h,c} R_D(y_{h,c}) * exp((eta_tau / ploidy) * u_{h,c})
#
# so each factor is tilted at `eta_tau / ploidy` with its own share-weighted
# beta, and `log Z` is the SUM of the factors' own normalisers. Two bounds follow
# and both are gates: `u` stays in [0,1], and `log Z` stays in [0, eta_tau].
#
# ⚠ ONE INVARIANT CATCHES BOTH COMPOSITION ERRORS, which is why it is asserted
# rather than documented: `sum over all factors of beta.sum() == ploidy`.
#   - Running each HAPLOTYPE at the full `eta_tau` instead of half doubles the
#     exponent, so `log Z` reaches 2 * eta_tau. This is the defect that was
#     ACCEPTED at every `tau` before the diploid bracket was added.
#   - Normalising `beta` to 1 within each CHROMOSOME instead of to `w_c` makes
#     the betas sum to `ploidy * n_chrom`, drives `u` to `n_chrom`, and turns one
#     release into `n_chrom` independent releases of the same person -- which
#     standing constraint 5 forbids. NO GATE IN THE REPO COULD DETECT THIS
#     before 2026-10-10.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Factor:
    """One independent chain of a composite release.

    `beta` sums to this factor's SHARE of the single global budget, NOT to 1.
    A factor whose beta sums to 1 is only correct for a haploid, single-chromosome
    release -- the degenerate case where share = 1 and ploidy = 1.
    """

    base: ExactBaseline
    target: np.ndarray       # the private input restricted to this factor
    beta: np.ndarray         # weights; sums to this factor's share
    label: str = ""

    @property
    def share(self) -> float:
        return float(self.beta.sum())


def composite_release(factors: "list[Factor]", eta: float, *, ploidy: int
                      ) -> tuple[np.ndarray, float, np.ndarray]:
    """`Q` over the PRODUCT support, `log Z`, and the per-factor `u` ranges.

    Returns `(Q_joint, log_z_total, shares)`. `Q_joint` is flattened in
    `itertools.product` order over the factors' supports, so `tv` and
    `ratio_spread` apply to it UNCHANGED -- the distance functions are not new
    code, which matters because new distance code checked by a new gate would be
    two unverified things agreeing.

    ⚠ The budget invariant is asserted here, at the only place that can see every
    factor at once. A per-factor check cannot catch it: each factor in isolation
    looks perfectly well formed.
    """
    if ploidy < 1:
        raise ValueError(f"ploidy must be >= 1, got {ploidy}")
    if not factors:
        raise ValueError("a release needs at least one factor")

    total_beta = sum(f.share for f in factors)
    if abs(total_beta - ploidy) > 1e-9:
        raise AssertionError(
            f"the factors' beta weights sum to {total_beta!r}, not ploidy = "
            f"{ploidy}. Each haplotype's weights must sum to 1 ACROSS its "
            f"chromosomes, so the total over {len(factors)} factors is the "
            f"ploidy. Shares seen: "
            f"{[(f.label or i, f.share) for i, f in enumerate(factors)]}.\n"
            f"  If each chromosome was normalised to 1 independently, this is "
            f"{len(factors)} independent releases of one genome, which compose "
            f"and degrade the guarantee (standing constraint 5)."
        )

    exponent = eta / ploidy
    per_factor: list[np.ndarray] = []
    log_z_total = 0.0
    for f in factors:
        u = utility_vector(f.base, f.target, f.beta)
        q, log_z = release_distribution(f.base, u, exponent)
        per_factor.append(q)
        log_z_total += log_z

    joint = per_factor[0]
    for q in per_factor[1:]:
        joint = np.outer(joint, q).ravel()
    return joint, float(log_z_total), np.array([f.share for f in factors])


def composite_utility(factors: "list[Factor]", *, ploidy: int) -> np.ndarray:
    """`u(p, y)` over the product support: the mean over haplotypes of the sum
    over chromosomes, assembled factor by factor.

    Kept separate from `composite_release` so the [0,1] bound can be checked on
    the composite utility itself, which is the quantity standing constraint 2 is
    about. `release_distribution` only ever sees one factor's `u`, which is in
    [0, share] and therefore passes its own bound even when the composite does
    not.
    """
    per_factor = [utility_vector(f.base, f.target, f.beta) for f in factors]
    acc = per_factor[0]
    for u in per_factor[1:]:
        acc = np.add.outer(acc, u).ravel()
    return acc / ploidy


def assert_log_z_bracket(log_z: float, eta: float, *, tol: float = 1e-9,
                         what: str = "log_z_p") -> None:
    """`log Z` in [0, eta_tau], for a composite release of any shape.

    ⚠ The bracket is on the TOTAL, never on a factor. A per-factor bracket is
    about `eta_tau / ploidy * w_c` and is roughly `ploidy * n_chrom` times
    tighter, so checking factors against the full `eta_tau` is a bracket that
    cannot fire -- the exact shape of the defect this module was written for.
    """
    if not (-tol <= log_z <= eta + tol):
        raise AssertionError(
            f"{what} = {log_z!r} is outside [0, eta_tau] = [0, {eta!r}].\n"
            f"  Z_p = sum_y R_D(y) exp(eta_tau * u(p,y)) with u in [0,1], so Z_p "
            f"is in [1, e^eta_tau] and its log is in [0, eta_tau] by definition. "
            f"A violation means the exponent applied is not the one calibrated "
            f"to tau -- most likely the full eta_tau where eta_tau/ploidy "
            f"belongs, or a beta normalised per chromosome rather than to its "
            f"share."
        )
