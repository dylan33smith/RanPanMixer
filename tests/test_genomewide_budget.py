"""B-IMP-genomewide: the budget is one per genome, not one per chromosome.

⚠ WHAT CHANGED AND WHY IT NEEDED A GATE OF ITS OWN. Until 2026-10-10
`utility.beta_v` hard-asserted `sum beta == 1` over whatever positions it was
handed, so "one `tau` per chromosome" was not merely the default -- it was the
only thing the code could express. Decision B makes `beta_v` sum to 1 across the
22 autosomes TOGETHER, each chromosome taking a SHARE. That turns a thing the
code could not say into a thing it can say wrongly, in one specific way:
normalising every chromosome to the whole budget, which is `n_chrom` independent
releases of one person and compounds the guarantee (standing constraint 5).

`tests/test_diploid_gate.py` already gates that error against the exact
reference. This file gates it against the MECHANISM, which only became possible
once the share existed -- the gate was deliberately built first.

THE THREE THINGS CHECKED HERE:
  A. the share arithmetic, and that it is exactly backward compatible at share 1
  B. the completeness checks, each shown to fire
  C. the composition error injected into the real sampler, end to end
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mechanism import calibration as C  # noqa: E402
from mechanism import ffbs as F  # noqa: E402
from mechanism import hmm as H  # noqa: E402
from mechanism import toy as TOY  # noqa: E402
from mechanism import utility  # noqa: E402
from mechanism import utility as U  # noqa: E402

TAUS = (0.05, 0.25, 0.50, 0.90)


@pytest.fixture(scope="module")
def two_chrom():
    """Two toy chromosomes with declared unequal shares, and their live models."""
    a, b = TOY.toy_segments()
    w = TOY.segment_shares()
    out = []
    for name, t, share in (("chrA", a, w[0]), ("chrB", b, w[1])):
        util = U.beta_v(t.support, share=share, chrom=name)
        m = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id=name)
        out.append((name, t, util, m, F.baseline_log_norm(m)))
    return out


# ---------------------------------------------------------------------------
# A. The share arithmetic
# ---------------------------------------------------------------------------

def test_share_one_is_byte_identical_to_the_old_behaviour():
    """BACKWARD COMPATIBILITY, and it is the reason `share` defaults to 1.0.

    Every number on record was produced by a `beta_v` that summed to 1. If adding
    the share changed that path by even a float, every recorded chr21 figure would
    silently stop being comparable to anything produced afterwards -- and nothing
    would raise, because both would be internally consistent.
    """
    t = TOY.toy_cohort()
    old_style = U.beta_v(t.support)
    explicit = U.beta_v(t.support, share=1.0)
    assert np.array_equal(old_style.beta, explicit.beta), (
        "beta_v(support) and beta_v(support, share=1.0) differ; the default is "
        "not the identity and every pre-2026-10-10 number is now incomparable")
    assert old_style.share == 1.0 and old_style.chrom == ""
    assert F.utility_id(old_style) == "indicator/w=inv_support/beta=chain", (
        "the utility_id for a whole-budget chain changed, so no new release can "
        "be compared to a recorded one")


def test_a_share_scales_the_weights_exactly_and_nothing_else():
    """`beta` at share w is w times `beta` at share 1, position for position.

    The shape of the weighting -- which positions matter relative to which -- is a
    property of the cohort (`w_v = 1/support(v)`) and must not move when the
    budget is split. If it did, two chromosomes would be weighting their positions
    on different principles while claiming one utility function.
    """
    t = TOY.toy_cohort()
    full = U.beta_v(t.support)
    for w in (0.7, 0.3, 1.0 / 22.0):
        part = U.beta_v(t.support, share=w, chrom="chrX")
        assert np.allclose(part.beta, full.beta * w, rtol=0, atol=1e-18)
        assert abs(part.beta.sum() - w) < 1e-12
        assert part.share == w


def test_the_share_is_stamped_on_the_utility_id():
    """Standing constraint 8: a chromosome at share 1 and the same chromosome at
    its genome-wide share are different utility functions and must not wear one
    name."""
    t = TOY.toy_cohort()
    a = F.utility_id(U.beta_v(t.support, share=0.7, chrom="chr1"))
    b = F.utility_id(U.beta_v(t.support, share=0.3, chrom="chr1"))
    c = F.utility_id(U.beta_v(t.support))
    assert a != b != c and a != c, f"utility ids collide: {a!r} {b!r} {c!r}"
    assert "share=0.7" in a and "chrom=chr1" in a


@pytest.mark.parametrize("share", [0.0, -0.1, 1.5, 2.0])
def test_an_impossible_share_is_refused(share):
    t = TOY.toy_cohort()
    with pytest.raises(ValueError, match="share must be in"):
        U.beta_v(t.support, share=share)


# ---------------------------------------------------------------------------
# B. The completeness checks, each shown to FIRE
# ---------------------------------------------------------------------------

def test_shares_summing_to_one_pass_and_summing_to_n_do_not():
    """THE CENTRAL CHECK. 22 chromosomes each normalised to the whole budget is
    22 releases of one person, and every one of those 22 `Utility` objects is
    individually valid -- which is exactly why this cannot live in `beta_v`."""
    t = TOY.toy_cohort()
    good = [U.beta_v(t.support, share=1.0 / 22.0, chrom=f"chr{i}")
            for i in range(1, 23)]
    assert abs(U.assert_shares_complete(good) - 1.0) < 1e-12

    bad = [U.beta_v(t.support, share=1.0, chrom=f"chr{i}") for i in range(1, 23)]
    with pytest.raises(AssertionError, match="sum to"):
        U.assert_shares_complete(bad)


def test_a_duplicated_chromosome_is_caught_even_when_the_total_is_right():
    """⚠ A total of 1 is NOT sufficient. Counting one chromosome twice while
    omitting another leaves the sum correct and the release wrong: the omitted
    chromosome is released at no privacy cost and the duplicated one at double."""
    t = TOY.toy_cohort()
    dup = [U.beta_v(t.support, share=0.5, chrom="chr1"),
           U.beta_v(t.support, share=0.5, chrom="chr1")]
    assert abs(sum(u.share for u in dup) - 1.0) < 1e-12, "the total is right"
    with pytest.raises(AssertionError, match="appears twice"):
        U.assert_shares_complete(dup)


def test_an_unnamed_chromosome_in_a_multi_chromosome_budget_is_refused():
    """Without a name the duplicate check cannot work, so the gap is refused
    rather than left to be discovered by its consequences."""
    t = TOY.toy_cohort()
    anon = [U.beta_v(t.support, share=0.5, chrom="chr1"),
            U.beta_v(t.support, share=0.5)]
    with pytest.raises(AssertionError, match="no `chrom`"):
        U.assert_shares_complete(anon)


# ---------------------------------------------------------------------------
# C. The real sampler, end to end -- and the composition error injected into it
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tau", TAUS)
def test_a_two_chromosome_genome_release_is_bounded(two_chrom, tau):
    """One release of one genome across two chromosomes, through `ffbs.release`.

    The whole point in one assertion: `log Z` summed over chromosomes lands in
    [0, eta_tau] and `u` in [0,1], even though each chromosome ran its own
    independent pair of chains.
    """
    eta = C.eta_tau(tau)
    parts = []
    for name, t, util, m, base in two_chrom:
        tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
        rel = F.release(m, t.haplotypes, tgt, util, tau, seed=11,
                        target_id="toytgt", graph_id="toy", baseline_norm=base)
        # the per-chromosome bracket is its OWN share of the budget
        assert -1e-9 <= rel.log_z_p <= eta * util.share + 1e-6, (
            f"{name}: log_z_p {rel.log_z_p!r} exceeds its share of the budget "
            f"{eta * util.share!r}")
        parts.append(rel)

    g = F.GenomeRelease(parts=tuple(parts))
    assert -1e-9 <= g.log_z_p <= eta + 1e-6
    assert 0.0 - 1e-9 <= g.u_achieved <= 1.0 + 1e-9
    assert g.u_achieved <= g.u_ceiling + 1e-12, (
        "the release achieved more utility than its own ceiling, which means the "
        "ceiling is computed on a different weighting than the draw")
    rec = g.record()
    assert rec["scope"] == "genome" and rec["n_chromosomes"] == 2
    assert abs(math.fsum(rec["shares"]) - 1.0) < 1e-12
    assert len(set(rec["parts"])) == 2, (
        "the two parts share a release_id, so one would overwrite the other")


@pytest.mark.parametrize("tau", TAUS)
def test_MUTATION_each_chromosome_at_the_whole_budget_is_caught(two_chrom, tau):
    """⚠⚠ THE MUTATION THIS CHUNK EXISTS FOR, now injectable into the mechanism.

    Before the share existed this error could not be written down: `beta_v`
    refused anything but `sum == 1`, so there was no way to express "this
    chromosome holds 70% of the budget" and therefore no way to express getting it
    wrong. Now the wrong version is one keyword away and is the most plausible
    genome-wide mistake there is, because normalising to 1 is what the
    single-chromosome code has always done.

    Each chromosome is released at `share = 1.0` -- individually valid, and
    `ffbs.release` accepts each one without complaint, which is the point.
    """
    eta = C.eta_tau(tau)
    parts, utils = [], []
    for name, t, _util, m, base in two_chrom:
        bad = U.beta_v(t.support, share=1.0, chrom=name)   # the mutation
        utils.append(bad)
        tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
        rel = F.release(m, t.haplotypes, tgt, bad, tau, seed=11,
                        target_id="toytgt", graph_id="toy", baseline_norm=base)
        parts.append(rel)

    # ⚠ Each PART passed its own bracket. That is the defect, not an oversight:
    # a per-chromosome check cannot see a budget spent twice.
    for rel, u in zip(parts, utils):
        assert -1e-9 <= rel.log_z_p <= eta * u.share + 1e-6

    fired = set()
    try:
        U.assert_shares_complete(utils)
    except AssertionError:
        fired.add("assert_shares_complete")
    try:
        F.GenomeRelease(parts=tuple(parts))
    except AssertionError:
        fired.add("GenomeRelease")

    assert fired == {"assert_shares_complete", "GenomeRelease"}, (
        f"at tau = {tau} the composition error was caught by {fired or 'NOTHING'}. "
        f"Both the utility-level and the release-level check must fire: the first "
        f"is what a driver can call before spending compute, the second is what "
        f"catches a driver that forgot to.")


def test_a_genome_release_refuses_parts_from_different_releases(two_chrom):
    """One release means one `(target, tau, seed, utility)`. Parts drawn at
    different `tau` or different seeds are separate releases of one genome, and
    separate releases compose."""
    parts = []
    for i, (name, t, util, m, base) in enumerate(two_chrom):
        tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
        parts.append(F.release(m, t.haplotypes, tgt, util,
                               0.25 if i == 0 else 0.5,     # different tau
                               seed=11, target_id="toytgt", graph_id="toy",
                               baseline_norm=base))
    with pytest.raises(AssertionError, match="disagree on"):
        F.GenomeRelease(parts=tuple(parts))


def test_a_genome_release_refuses_a_duplicated_chromosome(two_chrom):
    """Same trap as the utility-level check, one level up: the shares can total 1
    while one chromosome is counted twice and another is missing."""
    name, t, util, m, base = two_chrom[0]
    tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    half = U.beta_v(t.support, share=0.5, chrom="chrA")
    rel = F.release(m, t.haplotypes, tgt, half, 0.25, seed=11,
                    target_id="toytgt", graph_id="toy", baseline_norm=base)
    with pytest.raises(AssertionError, match="more than once"):
        F.GenomeRelease(parts=(rel, rel))


# ---------------------------------------------------------------------------
# D. THE ASSERTIONS THEMSELVES MUST FIRE
#
# ⚠ WHY THIS SECTION EXISTS, and it is the lesson of chunk 2. Three mutations of
# the new budget machinery -- reverting the diploid bracket to the full `eta`,
# dropping the budget from `_draw_one`, and disabling `GenomeRelease`'s share
# total -- ALL PASSED the tests above. The reason is general and worth stating:
#
#   LOOSENING AN ASSERTION CANNOT BE DETECTED BY A TEST THAT INDEPENDENTLY
#   CHECKS THE SAME BOUND.
#
# The tests above assert that the VALUES land in range. A mutation that only
# widens an assertion leaves every value untouched, so they still pass -- while
# the guard that was supposed to catch a real defect is now silent. Catching that
# needs tests that put the mechanism in a state where the bracket SHOULD fire and
# assert that it does.
#
# `GenomeRelease`'s share check had a second problem: it was hidden behind the
# genome-wide `log_z` check, which fires on the same input for a different reason.
# "GenomeRelease raised" is not evidence that the share check works. Each
# assertion is therefore isolated below so exactly one of them can fire.
# ---------------------------------------------------------------------------

def test_log_z_p_actually_enforces_the_budget_it_is_given():
    """The `budget` parameter is load-bearing, not decorative.

    Checked directly, because every other test observes `log_z_p` only through a
    value that happens to be in range.
    """
    t = TOY.toy_cohort()
    util = U.beta_v(t.support)
    m = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id=t.name)
    eta = C.eta_tau(0.9)
    psi = U.tilt_potentials(m.allowed, t.haplotypes, t.haplotypes[0], util.beta, eta)
    fwd = F.forward(m, psi)
    base = F.baseline_log_norm(m)

    lz = F.log_z_p(fwd, base, eta)                     # whole budget: fits
    assert 0.0 <= lz <= eta
    assert lz > 0.0, "the tilt did nothing, so this input cannot test a budget"

    # the same value against a budget it must exceed
    with pytest.raises(AssertionError, match=r"outside \[0, eta_tau \* share\]"):
        F.log_z_p(fwd, base, eta, budget=lz * 0.5)
    # and a budget it just fits inside
    F.log_z_p(fwd, base, eta, budget=lz * 1.000001)


def test_the_diploid_bracket_is_wired_to_the_share_not_to_eta():
    """A `Utility` whose declared `share` understates its own weights must be
    refused by the sampler.

    This is a real failure mode, not a contrivance: the share is a DECLARED
    quantity read from a manifest, so a release can be handed weights summing to
    one thing and a label saying another. If the bracket were against `eta`
    instead of `eta * share` it would accept this silently -- which is mutation M4
    of 2026-10-10, and nothing else in this file caught it.
    """
    import dataclasses
    t = TOY.toy_cohort()
    honest = U.beta_v(t.support)                       # beta sums to 1
    mislabelled = dataclasses.replace(honest, share=0.05, chrom="chr21")
    assert abs(mislabelled.beta.sum() - 1.0) < 1e-12, "the weights still sum to 1"
    assert mislabelled.share == 0.05, "but the label claims a twentieth of that"

    m = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id=t.name)
    base = F.baseline_log_norm(m)
    tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    with pytest.raises(AssertionError, match="eta_tau \\* share"):
        F.release(m, t.haplotypes, tgt, mislabelled, 0.9, seed=3,
                  target_id="toy", baseline_norm=base)


def test_the_genome_share_total_is_checked_in_isolation(two_chrom):
    """⚠ ISOLATED so that ONLY the share total can fire.

    The obvious way to test this -- two chromosomes each at share 1 -- also
    doubles the total `log_z`, so the genome-wide `log_z` bracket fires too and
    "GenomeRelease raised" proves nothing about the share check. That is how
    mutation M6 passed on 2026-10-10.

    Here both shares are UNDERSTATED (0.3 + 0.3 = 0.6), so every `log_z` is
    comfortably inside its bracket and the only thing wrong is the total. If this
    test stops failing when the share check is removed, the check is decoration.
    """
    parts = []
    for name, t, _util, m, base in two_chrom:
        util = U.beta_v(t.support, share=0.3, chrom=name)
        tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
        parts.append(F.release(m, t.haplotypes, tgt, util, 0.9, seed=11,
                               target_id="toytgt", graph_id="toy",
                               baseline_norm=base))
    eta = C.eta_tau(0.9)
    total_lz = math.fsum(p.log_z_p for p in parts)
    assert 0.0 <= total_lz <= eta, (
        f"the genome log_z {total_lz!r} is NOT inside [0, {eta!r}], so the "
        f"log_z bracket could fire and this test would stop isolating the "
        f"share check")
    assert abs(sum(p.share for p in parts) - 0.6) < 1e-12

    with pytest.raises(AssertionError, match="shares sum to"):
        F.GenomeRelease(parts=tuple(parts))


def test_the_genome_log_z_bracket_is_checked_in_isolation(two_chrom):
    """The mirror of the test above: shares that sum to 1 but a `log_z` total that
    does not fit, so only the `log_z` bracket can fire.

    Built by handing each part an `eta_tau` label from a smaller `tau` than the
    one it was drawn at -- which is what a release stamped with the wrong privacy
    parameter looks like, and is the one thing the share check cannot see.
    """
    import dataclasses
    parts = []
    for name, t, util, m, base in two_chrom:
        tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
        rel = F.release(m, t.haplotypes, tgt, util, 0.9, seed=11,
                        target_id="toytgt", graph_id="toy", baseline_norm=base)
        parts.append(dataclasses.replace(rel, eta_tau=C.eta_tau(0.05)))
    assert abs(sum(p.share for p in parts) - 1.0) < 1e-12, "the shares are correct"
    with pytest.raises(AssertionError, match="genome-wide log_z_p"):
        F.GenomeRelease(parts=tuple(parts))


def test_the_per_chain_bracket_is_derived_from_the_weights(monkeypatch):
    """`_draw_one` must bracket against `eta * sum(beta)`, not against `eta`.

    ⚠ THIS ONE NEEDS THE ERROR INJECTED, and the reason is worth recording. On
    correct inputs `log Z <= eta * sum(beta)` is a THEOREM, so no choice of target
    or `tau` can push a correct forward pass past it -- which means no black-box
    test can tell the tight bracket from the loose one. That is exactly why
    mutation M5 of 2026-10-10 (dropping the budget from `_draw_one`) survived every
    other test in this file.

    The bracket's real job is to catch a BROKEN FORWARD PASS, so the forward pass
    is what gets broken here: `log_norm` is inflated by an amount that sits between
    this chain's share of the budget and the whole budget. The tight bracket must
    refuse it; the loose one would not notice.
    """
    t = TOY.toy_cohort()
    share = 0.1
    util = U.beta_v(t.support, share=share, chrom="chrA")
    m = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id=t.name)
    base = F.baseline_log_norm(m)
    tau = 0.9
    eta = C.eta_tau(tau)
    exponent = eta / 2.0                       # what the diploid path passes down

    # an inflation the loose bracket accepts and the tight one must not
    bogus = exponent * share + 0.5 * (exponent - exponent * share)
    assert exponent * share < bogus < exponent, "the injected value proves nothing"

    real_forward = F.forward

    def inflated(hmm, potentials):
        fwd = real_forward(hmm, potentials)
        return F.Forward(alpha=fwd.alpha, log_c=fwd.log_c + (base + bogus) / len(fwd.log_c))

    monkeypatch.setattr(F, "forward", inflated)
    tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    with pytest.raises(AssertionError, match="eta_tau \\* share"):
        F.release(m, t.haplotypes, tgt, util, tau, seed=5,
                  target_id="toy", baseline_norm=base)


# ----------------------------------------------------------------------------------
# sum(beta) <= 1 on the TILT path, added 2026-10-12 after Dylan proposed accumulating
# the weights of positions where the sampler did not switch donors.
#
# ⚠⚠ WHY THIS GATE EXISTS AND WHY IT HAD TO BE HERE. `beta_v` asserts sum == share, but
# NOTHING REQUIRES A CALLER TO GO THROUGH `beta_v` -- `scripts/beta_concentration.py`
# legitimately hand-builds its allocations. And the one downstream check that looks like
# it would catch a mis-normalised beta cannot: `ffbs._draw_one` hands `log_z_p` a bracket
# of `eta * float(beta.sum())`, which WIDENS in lockstep with the very error it would
# need to catch. Mutation-verified 2026-10-12: with this gate removed, beta summing to
# 2.0 flows through `tilt_potentials` and the release carries a TRUE bound of
# tanh(eta * 2) = 0.80 while being STAMPED tau = 0.5.
#
# ⚠ The invariant is sum(beta) <= 1, not == 1, because a per-chromosome release carries
# that chromosome's SHARE and the 22 shares sum to 1 (decision J).
# ----------------------------------------------------------------------------------

def _tiny():
    hap = np.array([[0, 1, 0], [0, 0, 1]], dtype=np.int8)
    tgt = np.array([0, 1, 0], dtype=np.int8)
    allowed = np.ones((2, 3), dtype=bool)
    return allowed, hap, tgt


def test_tilt_refuses_beta_summing_above_one():
    """THE GATE. delta_u > 1 does not weaken the bound, it REMOVES it (constraint 2)."""
    allowed, hap, tgt = _tiny()
    beta = np.array([0.7, 0.7, 0.6])                      # sums to 2.0
    with pytest.raises(ValueError, match="above 1"):
        utility.tilt_potentials(allowed, hap, tgt, beta, math.atanh(0.5))


def test_tilt_accepts_a_full_budget_and_a_share():
    """CONTROL. Both legitimate cases must pass, or the gate is just blocking work."""
    allowed, hap, tgt = _tiny()
    eta = math.atanh(0.5)
    for beta in (np.array([0.2, 0.5, 0.3]),               # a whole budget, sums to 1
                 np.array([0.01, 0.02, 0.01])):           # a chromosome's share
        psi = utility.tilt_potentials(allowed, hap, tgt, beta, eta)
        assert np.isfinite(psi).all() and (psi >= 0).all()


def test_the_overrun_the_gate_prevents_is_quantified():
    """The gate's REASON, as arithmetic: an unnormalised beta silently raises tau.

    This is the number that makes the gate worth having rather than a style rule -- a
    release stamped tau = 0.5 while carrying sum(beta) = 1.979595 is really at TV 0.7959,
    and `p_succ_bound` cannot detect it because measured TV reaches only 27-43% of tau.
    """
    eta = math.atanh(0.5)
    assert math.isclose(math.tanh(eta * 1.0), 0.5, rel_tol=1e-12)        # the stamp
    assert math.isclose(math.tanh(eta * 1.979595), 0.795929, abs_tol=1e-6)
    assert math.isclose(math.tanh(eta * 2.0), 0.8, rel_tol=1e-12)
