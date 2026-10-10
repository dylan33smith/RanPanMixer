"""B-THY-diploid-gate: the falsification gate for the path that actually runs.

⚠ WHY THIS FILE EXISTS. Every clause in `test_toy_enumeration.py` is HAPLOID and
single-chromosome. The diploid path is the only one that has ever been run -- all
18 releases of the first chr21 run carry `ploidy: 2` -- and decision E makes Arm 1
diploid across 22 autosomes under one genome-wide budget. So the project's only
falsifier did not cover the only code path in use, and `exact.py` had no diploid
primitive to apply it with. That was missing REFERENCE CODE, not a missing
assertion.

WHAT IS GATED HERE, and the last one is the point:

A. DIPLOID, one chromosome, on the 96-output toy -- product support 9,216.
B. TWO SEGMENTS, haploid, with declared unequal shares -- EXHAUSTIVE over all
   1,024 input pairs, not sampled.
C. DIPLOID x TWO SEGMENTS -- four independent chains, the Arm 1 shape.
D. THREE MUTATIONS. Each injects one real calibration error and asserts which
   clauses fire, AT EVERY `tau` on the grid. A gate that cannot fail is
   decoration, and clause (b) of the haploid gate is the project's own worked
   example: exact TV reaches only 27-43% of `tau`, so a doubled exponent passes
   it everywhere. The mutations are the evidence; the passing clauses are not.

⚠ THE INPUT SPACE IS THE HARD PART, NOT THE SUPPORT. The big toy has 2,673
fully-called haploid inputs, so 7.1 million diploid inputs and 5.1e13 pairs --
the exhaustive sweep the haploid gate runs is not available. So: CONSTRUCTED
extremal pairs, which saturate the bound and are therefore sharp, plus a seeded
random sample whose size and seed are recorded here so the coverage is stated
rather than implied. Where the space is small enough -- the two-segment haploid
case -- the sweep IS exhaustive and says so.
"""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mechanism import calibration, exact, toy  # noqa: E402

# The pre-registered diagnostic grid. `tau` = 0 is the floor and must reproduce
# the untilted product exactly; 0.9 is past any claim we would make and is here
# because a bound must be checked where it is loose.
TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90)

# Sampled coverage, stated rather than implied.
PAIR_SAMPLE = 400
PAIR_SEED = 20261010


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

def _beta(support: np.ndarray, share: float) -> np.ndarray:
    """Inverse-support weights scaled to sum to `share`, not to 1.

    Mirrors the real `utility.beta_v` shape (non-uniform, `w_v = 1/support(v)`)
    so the gate exercises uneven weights, then scales to the segment's SHARE of
    the one global budget. A `beta` summing to 1 per segment is mutation 2.
    """
    w = 1.0 / support.astype(np.float64)
    return w * (share / w.sum())


@pytest.fixture(scope="module")
def big():
    """The 96-output toy, its baseline, and share-1 beta. Haploid building block."""
    t = toy.toy_cohort()
    base = exact.enumerate_baseline(t.haplotypes, t.genetic_pos)
    return t, base, _beta(t.support, 1.0)


@pytest.fixture(scope="module")
def segs():
    """Two segments with declared unequal shares, baselines, and share-scaled beta."""
    a, b = toy.toy_segments()
    w = toy.segment_shares()
    ba = exact.enumerate_baseline(a.haplotypes, a.genetic_pos)
    bb = exact.enumerate_baseline(b.haplotypes, b.genetic_pos)
    return (a, b), (ba, bb), w, (_beta(a.support, w[0]), _beta(b.support, w[1]))


def _diploid_factors(base, beta, p1, p2):
    """One chromosome, two haplotypes: each at share 1, so the betas sum to 2."""
    return [exact.Factor(base, p1, beta, "hap0"),
            exact.Factor(base, p2, beta, "hap1")]


def _segmented_factors(bases, betas, targets, *, ploidy):
    """`ploidy` x `n_seg` factors. Each haplotype's betas sum to 1 across segments,
    so the total over all factors is the ploidy."""
    out = []
    for h in range(ploidy):
        for c, (bse, bet) in enumerate(zip(bases, betas)):
            out.append(exact.Factor(bse, targets[h][c], bet, f"hap{h}seg{c}"))
    return out


def _sharp_diploid(t):
    """The diploid pair that SATURATES the ratio bound.

    p = (D1, D1), q = (D2, D2). D1 and D2 differ at every position and both
    "copy D1 throughout" and "copy D2 throughout" are in the support, so
    `u_p - u_q` is +1 at one joint output and -1 at another. The range is exactly
    2 and the spread EQUALS `exp(epsilon_tau)`: attained, not approached.
    """
    d1, d2 = toy.sharp_pair(t)
    return (d1, d1), (d2, d2)


def _sample_pairs(inputs, n, seed):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(inputs), size=(n, 2))
    return [(inputs[i], inputs[j]) for i, j in idx if i != j]


# ---------------------------------------------------------------------------
# A. DIPLOID, ONE CHROMOSOME
# ---------------------------------------------------------------------------

def test_diploid_tau_zero_reproduces_the_product_baseline(big):
    """The floor. At `tau` = 0 the release is the untilted product of two
    independent draws from `R_D`, which is the only point where the answer is
    known in closed form without the mechanism."""
    t, base, beta = big
    p, _ = _sharp_diploid(t)
    q, log_z, _ = exact.composite_release(_diploid_factors(base, beta, *p),
                                          0.0, ploidy=2)
    expected = np.outer(base.prob, base.prob).ravel()
    assert np.allclose(q, expected, atol=1e-15), "tau = 0 is not the product baseline"
    assert abs(log_z) < 1e-15, f"log_z_p at tau = 0 is {log_z!r}, not 0"


def test_diploid_composite_utility_stays_in_the_unit_interval(big):
    """Standing constraint 2 on the COMPOSITE, which is the quantity the
    calibration is about. Each factor's own `u` is in [0, share] and passes its
    own bound even when the composite does not -- which is mutation 2."""
    t, base, beta = big
    inputs = toy.toy_inputs(t, kind="all")
    rng = np.random.default_rng(PAIR_SEED)
    for i, j in rng.integers(0, len(inputs), size=(60, 2)):
        u = exact.composite_utility(
            _diploid_factors(base, beta, inputs[i], inputs[j]), ploidy=2)
        assert u.min() >= -1e-12 and u.max() <= 1 + 1e-12, (
            f"composite u left [0,1]: [{u.min():.6g}, {u.max():.6g}]")


@pytest.mark.parametrize("tau", TAUS)
def test_diploid_log_z_is_bracketed(big, tau):
    """`log Z` in [0, eta_tau] for the PAIR, never per haplotype. The per-haplotype
    value is in [0, eta_tau/2], so checking a haplotype against the full eta_tau
    is a bracket that cannot fire."""
    t, base, beta = big
    eta = calibration.eta_tau(tau)
    inputs = toy.toy_inputs(t, kind="all")
    rng = np.random.default_rng(PAIR_SEED + 1)
    for i, j in rng.integers(0, len(inputs), size=(40, 2)):
        _, log_z, _ = exact.composite_release(
            _diploid_factors(base, beta, inputs[i], inputs[j]), eta, ploidy=2)
        exact.assert_log_z_bracket(log_z, eta, what=f"diploid log_z_p at tau={tau}")


@pytest.mark.parametrize("tau", TAUS)
def test_diploid_exact_tv_never_exceeds_tau(big, tau):
    """Theorem 1 on the diploid release object, exactly, over the product support."""
    t, base, beta = big
    eta = calibration.eta_tau(tau)
    p, q = _sharp_diploid(t)
    pairs = [(p, q)] + _sample_pairs(
        [(a, b) for a, b in zip(*[toy.toy_inputs(t)] * 2)][:0] or
        list(zip(toy.toy_inputs(t)[:60], toy.toy_inputs(t)[::-1][:60])),
        PAIR_SAMPLE // 8, PAIR_SEED + 2)
    worst = 0.0
    for pp, qq in pairs:
        q1, _, _ = exact.composite_release(_diploid_factors(base, beta, *pp), eta, ploidy=2)
        q2, _, _ = exact.composite_release(_diploid_factors(base, beta, *qq), eta, ploidy=2)
        worst = max(worst, exact.tv(q1, q2))
    assert worst <= tau + 1e-12, f"exact diploid TV {worst!r} exceeds tau = {tau}"


@pytest.mark.parametrize("tau", TAUS[1:])
def test_the_diploid_sharp_pair_saturates_the_ratio_bound(big, tau):
    """THE SHARP CLAUSE. The constructed diploid pair attains `exp(epsilon_tau)`
    exactly, so ANY inflation of the exponent fails it immediately. This is the
    only clause that catches a doubled tilt -- TV does not."""
    t, base, beta = big
    eta, eps = calibration.eta_tau(tau), calibration.epsilon_tau(tau)
    p, q = _sharp_diploid(t)
    q1, _, _ = exact.composite_release(_diploid_factors(base, beta, *p), eta, ploidy=2)
    q2, _, _ = exact.composite_release(_diploid_factors(base, beta, *q), eta, ploidy=2)
    spread = exact.ratio_spread(q1, q2)
    assert spread <= math.exp(eps) * (1 + 1e-9), (
        f"ratio_spread {spread!r} exceeds exp(epsilon_tau) = {math.exp(eps)!r}")
    assert spread >= math.exp(eps) * (1 - 1e-9), (
        f"ratio_spread {spread!r} does NOT saturate exp(epsilon_tau) = "
        f"{math.exp(eps)!r}. The pair has been blunted and the gate is no longer "
        f"sharp -- it would stop detecting a doubled exponent.")


# ---------------------------------------------------------------------------
# B. TWO SEGMENTS, HAPLOID -- EXHAUSTIVE
# ---------------------------------------------------------------------------

def _seg_inputs(tt):
    """Every composite haploid input: the product of the two segments' inputs."""
    ia, ib = toy.toy_inputs(tt[0]), toy.toy_inputs(tt[1])
    return [(a, b) for a in ia for b in ib]


@pytest.mark.parametrize("tau", TAUS)
def test_two_segment_budget_is_exhaustively_bounded(segs, tau):
    """EXHAUSTIVE over all 32 composite inputs and all 1,024 ordered pairs.

    Three clauses at once, because they share the expensive enumeration: `log Z`
    bracketed on the TOTAL, exact TV <= `tau`, and `ratio_spread` <=
    `exp(epsilon_tau)`. Also checks the tighter PER-SEGMENT statement
    `log_z_c <= eta_tau * w_c`, which is what makes the global bracket meaningful
    rather than merely true.
    """
    tt, bases, w, betas = segs
    eta, eps = calibration.eta_tau(tau), calibration.epsilon_tau(tau)
    inputs = _seg_inputs(tt)
    assert len(inputs) == 32, f"input space changed: {len(inputs)}"

    dists = []
    for tgt in inputs:
        factors = _segmented_factors(bases, betas, [tgt], ploidy=1)
        q, log_z, shares = exact.composite_release(factors, eta, ploidy=1)
        exact.assert_log_z_bracket(log_z, eta, what=f"two-segment log_z at tau={tau}")
        # the per-segment statement: each factor's own normaliser is bounded by
        # its OWN share of the budget, not by the whole of it
        for f, share in zip(factors, shares):
            u = exact.utility_vector(f.base, f.target, f.beta)
            _, lz_c = exact.release_distribution(f.base, u, eta)
            assert -1e-12 <= lz_c <= eta * share + 1e-12, (
                f"segment {f.label} log_z {lz_c!r} exceeds its share of the "
                f"budget, eta_tau * w_c = {eta * share!r}")
        dists.append(q)

    worst_tv, worst_spread = 0.0, 1.0
    for i in range(len(dists)):
        for j in range(len(dists)):
            if i == j:
                continue
            worst_tv = max(worst_tv, exact.tv(dists[i], dists[j]))
            worst_spread = max(worst_spread, exact.ratio_spread(dists[i], dists[j]))
    assert worst_tv <= tau + 1e-12, f"exact TV {worst_tv!r} exceeds tau = {tau}"
    assert worst_spread <= math.exp(eps) * (1 + 1e-9), (
        f"ratio_spread {worst_spread!r} exceeds exp(epsilon_tau) = {math.exp(eps)!r}")


@pytest.mark.parametrize("tau", TAUS[1:])
def test_the_two_segment_sharp_pair_saturates(segs, tau):
    """The saturating pair generalises across segments: match donor 0 on BOTH
    segments against donor 1 on both. `u_p - u_q` spans w_0 + w_1 = 1 in each
    direction, so the spread attains `exp(epsilon_tau)` whatever the shares are --
    which is also why unequal shares cannot hide a wrong normaliser here."""
    tt, bases, w, betas = segs
    eta, eps = calibration.eta_tau(tau), calibration.epsilon_tau(tau)
    p = [(tt[0].haplotypes[0].copy(), tt[1].haplotypes[0].copy())]
    q = [(tt[0].haplotypes[1].copy(), tt[1].haplotypes[1].copy())]
    q1, _, _ = exact.composite_release(_segmented_factors(bases, betas, p, ploidy=1),
                                       eta, ploidy=1)
    q2, _, _ = exact.composite_release(_segmented_factors(bases, betas, q, ploidy=1),
                                       eta, ploidy=1)
    spread = exact.ratio_spread(q1, q2)
    assert math.exp(eps) * (1 - 1e-9) <= spread <= math.exp(eps) * (1 + 1e-9), (
        f"the two-segment pair gives spread {spread!r}, not the saturating "
        f"exp(epsilon_tau) = {math.exp(eps)!r}")


# ---------------------------------------------------------------------------
# C. DIPLOID x TWO SEGMENTS -- the Arm 1 shape, four independent chains
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tau", TAUS)
def test_four_factor_release_is_bounded(segs, tau):
    """Diploid across two chromosomes: 2 x 2 = 4 chains, each at `eta_tau/2` with
    its segment's share. The total `log Z` must still land in [0, eta_tau] and the
    composite `u` in [0,1] -- the budget does not grow with the number of chains,
    which is the whole content of standing constraint 4."""
    tt, bases, w, betas = segs
    eta = calibration.eta_tau(tau)
    inputs = _seg_inputs(tt)
    rng = np.random.default_rng(PAIR_SEED + 3)
    for i, j in rng.integers(0, len(inputs), size=(40, 2)):
        factors = _segmented_factors(bases, betas, [inputs[i], inputs[j]], ploidy=2)
        assert len(factors) == 4
        u = exact.composite_utility(factors, ploidy=2)
        assert u.min() >= -1e-12 and u.max() <= 1 + 1e-12, (
            f"four-factor composite u left [0,1]: [{u.min():.6g}, {u.max():.6g}]")
        _, log_z, _ = exact.composite_release(factors, eta, ploidy=2)
        exact.assert_log_z_bracket(log_z, eta, what=f"4-factor log_z at tau={tau}")


@pytest.mark.parametrize("tau", TAUS)
def test_four_factor_exact_tv_never_exceeds_tau(segs, tau):
    """Theorem 1 on the Arm 1 release shape, exactly, over the 256-point product."""
    tt, bases, w, betas = segs
    eta = calibration.eta_tau(tau)
    inputs = _seg_inputs(tt)
    # constructed extremal pair first, then a seeded sample
    extremal = (([(tt[0].haplotypes[0].copy(), tt[1].haplotypes[0].copy())] * 2),
                ([(tt[0].haplotypes[1].copy(), tt[1].haplotypes[1].copy())] * 2))
    rng = np.random.default_rng(PAIR_SEED + 4)
    pairs = [extremal]
    for i, j, k, m in rng.integers(0, len(inputs), size=(PAIR_SAMPLE // 8, 4)):
        pairs.append(([inputs[i], inputs[j]], [inputs[k], inputs[m]]))
    worst = 0.0
    for pp, qq in pairs:
        q1, _, _ = exact.composite_release(
            _segmented_factors(bases, betas, pp, ploidy=2), eta, ploidy=2)
        q2, _, _ = exact.composite_release(
            _segmented_factors(bases, betas, qq, ploidy=2), eta, ploidy=2)
        worst = max(worst, exact.tv(q1, q2))
    assert worst <= tau + 1e-12, f"exact 4-factor TV {worst!r} exceeds tau = {tau}"


# ---------------------------------------------------------------------------
# D. THE MUTATIONS -- the evidence
#
# Each injects one real error and records WHICH clauses fire. A mutation that
# fires nothing means the gate above is decoration at that `tau`.
# ---------------------------------------------------------------------------

def _which_clauses_fire(factors, eta, eps, ploidy, *, declared_ploidy=None):
    """Run every clause and return the set that FAILED. Never raises."""
    fired = set()
    dp = ploidy if declared_ploidy is None else declared_ploidy
    try:
        exact.composite_release(factors, eta, ploidy=dp)
    except AssertionError:
        fired.add("budget_invariant")
    try:
        u = exact.composite_utility(factors, ploidy=dp)
        if u.min() < -1e-12 or u.max() > 1 + 1e-12:
            fired.add("u_in_unit_interval")
    except AssertionError:
        fired.add("u_in_unit_interval")
    if "budget_invariant" not in fired:
        _, log_z, _ = exact.composite_release(factors, eta, ploidy=dp)
        try:
            exact.assert_log_z_bracket(log_z, eta)
        except AssertionError:
            fired.add("log_z_bracket")
    return fired


@pytest.mark.parametrize("tau", TAUS[1:])
def test_mutation_full_eta_per_haplotype_is_caught(big, tau):
    """MUTATION 1: each haplotype tilted at the full `eta_tau` instead of half.

    This is the SUM aggregator in disguise and it is the defect that was ACCEPTED
    at every `tau` before a diploid bracket existed -- it returned a record whose
    published `log_z_p` violated its own stated invariant.
    """
    t, base, beta = big
    eta, eps = calibration.eta_tau(tau), calibration.epsilon_tau(tau)
    p, q = _sharp_diploid(t)
    factors = _diploid_factors(base, beta, *p)
    # the mutation: declare ploidy 1 while handing over two haplotypes' factors,
    # so the exponent applied is eta_tau rather than eta_tau/2
    fired = _which_clauses_fire(factors, eta, eps, 2, declared_ploidy=1)
    assert fired, (
        f"MUTATION 1 PASSED EVERY CLAUSE at tau = {tau}. A doubled diploid "
        f"exponent is undetectable and the gate is decoration.")
    # and the sharp clause must see it too, independently of the brackets
    q1, _, _ = exact.composite_release(factors, eta * 2, ploidy=2)
    q2, _, _ = exact.composite_release(_diploid_factors(base, beta, *q), eta * 2, ploidy=2)
    spread = exact.ratio_spread(q1, q2)
    assert spread > math.exp(eps) * (1 + 1e-9), (
        f"a doubled exponent gives spread {spread!r}, which still fits inside "
        f"exp(epsilon_tau) = {math.exp(eps)!r} -- the sharp clause missed it")


@pytest.mark.parametrize("tau", TAUS[1:])
def test_mutation_beta_normalised_per_segment_is_caught(segs, tau):
    """MUTATION 2: `beta` normalised to 1 within EACH chromosome instead of to its
    share `w_c`.

    ⚠ THIS IS THE ONE NOTHING IN THE REPO COULD DETECT before 2026-10-10. It is
    also the most plausible genome-wide error there is, because normalising to 1
    is exactly what the single-chromosome code does and the hard assertion in
    `utility.beta_v` demands. Its consequence is that one release becomes
    `n_chrom` independent releases of the same person, which compose.
    """
    tt, bases, w, betas = segs
    eta, eps = calibration.eta_tau(tau), calibration.epsilon_tau(tau)
    bad = [_beta(tt[0].support, 1.0), _beta(tt[1].support, 1.0)]  # each sums to 1
    tgt = [(tt[0].haplotypes[0].copy(), tt[1].haplotypes[0].copy())]
    factors = _segmented_factors(bases, bad, tgt, ploidy=1)
    fired = _which_clauses_fire(factors, eta, eps, 1)
    assert "budget_invariant" in fired, (
        f"MUTATION 2 did not trip the budget invariant at tau = {tau}. The "
        f"factors' betas sum to {sum(f.share for f in factors)!r} against a "
        f"ploidy of 1, so this is the clause that must fire.")
    assert "u_in_unit_interval" in fired, (
        f"MUTATION 2 left the composite utility inside [0,1] at tau = {tau}, so "
        f"standing constraint 2 would not have caught it either")


@pytest.mark.parametrize("tau", TAUS[1:])
def test_mutation_epsilon_in_place_of_eta_is_caught(segs, tau):
    """MUTATION 3: the tilt applied with `epsilon_tau` where `eta_tau` belongs.

    A one-character error (`epsilon_tau = 2 * eta_tau`). The haploid gate showed
    exact TV passes it at EVERY `tau` while `ratio_spread` catches it, and this
    asserts the same holds for a composite release -- i.e. that adding segments
    and haplotypes has not blunted the sharp clause.
    """
    tt, bases, w, betas = segs
    eta, eps = calibration.eta_tau(tau), calibration.epsilon_tau(tau)
    p = [(tt[0].haplotypes[0].copy(), tt[1].haplotypes[0].copy())]
    q = [(tt[0].haplotypes[1].copy(), tt[1].haplotypes[1].copy())]
    fp = _segmented_factors(bases, betas, p, ploidy=1)
    fq = _segmented_factors(bases, betas, q, ploidy=1)
    # the mutation: eps where eta belongs
    q1, log_z, _ = exact.composite_release(fp, eps, ploidy=1)
    q2, _, _ = exact.composite_release(fq, eps, ploidy=1)

    fired = set()
    try:
        exact.assert_log_z_bracket(log_z, eta)
    except AssertionError:
        fired.add("log_z_bracket")
    if exact.ratio_spread(q1, q2) > math.exp(eps) * (1 + 1e-9):
        fired.add("ratio_spread")
    if exact.tv(q1, q2) > tau + 1e-12:
        fired.add("tv")

    assert "ratio_spread" in fired, (
        f"MUTATION 3 was not caught by ratio_spread at tau = {tau}; it is the "
        f"clause that is supposed to be sharp")
    assert fired, f"MUTATION 3 passed every clause at tau = {tau}"


# ---------------------------------------------------------------------------
# E. THE CLAUSE THAT MAKES THIS A GATE ON THE MECHANISM
#
# ⚠ Everything above checks `exact.py` against itself, which is necessary and is
# NOT sufficient: a reference implementation agreeing with its own gate is two
# unverified things agreeing. These clauses tie the diploid reference to the
# diploid code path that actually runs.
#
# `log_z_p` is DETERMINISTIC -- no draw, no Monte-Carlo error -- so comparing it
# between `ffbs.release` and `exact.composite_release` pins the exponent, the
# factorisation across haplotypes and the normalisation, all to ~1e-12 and with
# zero sampling. That is a far stronger statement per unit of compute than any
# comparison of sampled frequencies, and it is exactly what the three mutations
# attack.
# ---------------------------------------------------------------------------

from mechanism import ffbs as F  # noqa: E402
from mechanism import hmm as H  # noqa: E402
from mechanism import utility as U  # noqa: E402


@pytest.fixture(scope="module")
def live():
    """The real mechanism on the toy: a CohortHMM, a Utility, and the baseline."""
    t = toy.toy_cohort()
    util = U.beta_v(t.support)
    m = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id=t.name)
    return t, util, m, F.baseline_log_norm(m)


@pytest.mark.parametrize("tau", TAUS)
def test_ffbs_diploid_log_z_matches_the_exact_reference(live, big, tau):
    """The mechanism's diploid `log_z_p` equals the exact composite reference.

    If `ffbs` applied the full `eta_tau` per haplotype rather than half, this
    would be off by a factor that grows with `eta_tau` -- mutation 1, detected
    without drawing anything. The two routes share no code: `exact` enumerates
    every output path with a Python loop and a dense transition matrix, `ffbs`
    does two O(TK) rank-one forward passes.
    """
    t, util, m, basenorm = live
    _, base, beta = big
    assert np.allclose(util.beta, beta, atol=1e-15), (
        "the live Utility's beta and the reference beta have diverged, so this "
        "comparison would no longer be like for like")

    inputs = toy.toy_inputs(t)
    rng = np.random.default_rng(PAIR_SEED + 5)
    for i, j in rng.integers(0, len(inputs), size=(12, 2)):
        tgt = np.stack([inputs[i], inputs[j]], axis=1)      # (T, 2), diploid
        assert tgt.shape == (t.T, 2)
        rel = F.release(m, t.haplotypes, tgt, util, tau, seed=7,
                        target_id="toy", baseline_norm=basenorm)
        _, log_z_ref, _ = exact.composite_release(
            _diploid_factors(base, beta, inputs[i], inputs[j]), rel_eta := 
            calibration.eta_tau(tau), ploidy=2)
        assert abs(rel.log_z_p - log_z_ref) < 1e-12, (
            f"ffbs diploid log_z_p = {rel.log_z_p!r} but the exact composite "
            f"reference gives {log_z_ref!r} at tau = {tau} (eta = {rel_eta!r}). "
            f"A mismatch here is the exponent, the per-haplotype factorisation, "
            f"or the baseline subtraction -- not sampling noise, because neither "
            f"side sampled.")


def test_ffbs_diploid_marginals_match_the_exact_single_chain(live, big):
    """Each haplotype's marginal output distribution is the haploid tilted
    distribution at `eta_tau/2`, and the two are independent.

    This is the sampled complement to the deterministic clause above: it checks
    that the DRAW follows the factorisation, not merely that the normaliser does.
    A shared or correlated RNG stream between the two haplotypes would pass the
    `log_z` clause and fail here -- and that defect is real history, the seed and
    seed + 1 bug of 2026-09-30 that made one release's haplotype 1 bit-identical
    to the next release's haplotype 0.
    """
    t, util, m, basenorm = live
    _, base, beta = big
    tau = 0.9                     # the loosest grid point: the largest tilt to see
    eta = calibration.eta_tau(tau)
    p1, p2 = toy.sharp_pair(t)
    tgt = np.stack([p1, p2], axis=1)

    n = 4000
    key = {row.tobytes(): i for i, row in enumerate(base.paths)}
    counts = np.zeros((2, base.M), dtype=np.int64)
    for s in range(n):
        rel = F.release(m, t.haplotypes, tgt, util, tau, seed=1000 + s,
                        target_id="toy", baseline_norm=basenorm)
        for h in range(2):
            counts[h, key[np.ascontiguousarray(
                rel.released_path[:, h]).tobytes()]] += 1

    for h, tgt_h in enumerate((p1, p2)):
        u = exact.utility_vector(base, tgt_h, beta)
        q_ref, _ = exact.release_distribution(base, u, eta / 2.0)
        emp = counts[h] / n
        # 3-sigma band on the multinomial, per cell with enough mass to judge
        sd = np.sqrt(np.maximum(q_ref, 1e-12) * (1 - q_ref) / n)
        off = np.abs(emp - q_ref) > 4.0 * sd + 1e-3
        assert not off.any(), (
            f"haplotype {h}'s sampled marginal departs from the exact tilted "
            f"distribution at eta_tau/2 in {off.sum()} of {base.M} cells "
            f"(worst |emp - ref| = {np.abs(emp - q_ref).max():.4g}, n = {n})")
