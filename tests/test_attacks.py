"""The arm guard, and the optimal attacker's positive controls.

⚠ The guard tests matter more than they look. The combinations they refuse do not
error on their own — they return a plausible number that reads as a privacy result.
So the test is not "does the code raise", it is "is the number that would otherwise
be reported unobtainable".

⚠ And the attacker tests are POSITIVE CONTROLS, not performance measurements. The
repository forecasts that every attack will sit at chance at the operating
configuration (R5: signal-to-noise about 1:900). Standing constraint 7 says that
null is not evidence unless the attacker is shown to work where there IS signal, so
these establish the conditions under which it must succeed.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import attacks as A  # noqa: E402
from mechanism import calibration as C  # noqa: E402
from mechanism import exact, toy  # noqa: E402

TAUS = (0.05, 0.10, 0.25, 0.50, 0.90)


# ---------------------------------------------------------------------------
# The arm guard
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["db_search", "membership"])
def test_an_arm1_shaped_attack_is_refused_in_arm2(kind):
    """THE GUARD. These two are ill posed in arm 2 and would return a number anyway.

    Database search: the target is external by construction, so no database row
    produced the release and a 0% link rate is the experimental design.
    Membership: there is no positive class, and the release is a donor mosaic, so a
    high match fraction identifies DONORS rather than the target.
    """
    with pytest.raises(ValueError, match="ILL POSED in arm 2"):
        A.assert_well_posed(kind, "arm2")
    A.assert_well_posed(kind, "arm1")        # the same attack IS posed in arm 1


def test_the_refusal_explains_itself_rather_than_just_failing():
    """A guard whose message does not say WHY gets worked around.

    The refusal has to carry the reason the number would be meaningless, because
    the next person's instinct on seeing a raise is to pass the flag that silences
    it. There is no such flag, and the message says what to do instead.
    """
    with pytest.raises(ValueError) as e:
        A.assert_well_posed("db_search", "arm2")
    msg = str(e.value)
    assert "NO CORRECT ANSWER" in msg
    assert "experimental design, not privacy" in msg
    assert "lrt_pair" in msg, "the message must name the well-posed alternative"
    with pytest.raises(ValueError) as e2:
        A.assert_well_posed("membership", "arm2")
    assert "identifies DONORS" in str(e2.value)


def test_the_pair_test_is_well_posed_in_both_arms():
    """The theorem's own game needs no adaptation: p is the target's path, q is any
    other admissible path. No database, no panel, no membership label."""
    for arm in A.ARMS:
        A.assert_well_posed("lrt_pair", arm)


def test_an_attack_without_a_positive_control_is_refused():
    """Standing constraint 7, enforced at construction.

    The repository's own forecast is a chance-level null for every attack here, so
    an attack with no stated success condition can never have its failure reported
    as privacy. Making the field mandatory is cheaper than remembering the rule.
    """
    with pytest.raises(ValueError, match="no positive control"):
        A.AttackSpec(name="lazy", arm="arm2", kind="lrt_pair",
                     statistic="Lambda", chance_level=0.5, positive_control="   ")
    ok = A.AttackSpec(name="A2-ATK-LRT", arm="arm2", kind="lrt_pair",
                      statistic="eta*(u_p - u_q) - (log Z_p - log Z_q)",
                      chance_level=0.5,
                      positive_control="must exceed (1+tau)/2 against a mechanism "
                                       "miscalibrated with epsilon_tau for eta_tau")
    assert ok.arm == "arm2"


def test_an_arm1_spec_cannot_be_relabelled_into_arm2():
    """The arm is part of the spec's identity, checked in __post_init__, so an
    Arm 1 attack cannot be re-pointed at Arm 2 by editing a field."""
    A.AttackSpec(name="link", arm="arm1", kind="db_search", statistic="gap score",
                 chance_level=1.0 / 44.0, positive_control="links the unobfuscated "
                 "target at near-100%")
    with pytest.raises(ValueError, match="ILL POSED in arm 2"):
        A.AttackSpec(name="link", arm="arm2", kind="db_search", statistic="gap score",
                     chance_level=1.0 / 44.0, positive_control="x")


# ---------------------------------------------------------------------------
# The optimal attacker: it attains the bound, and it is bounded by it
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def toy_pair():
    t = toy.toy_cohort()
    base = exact.enumerate_baseline(t.haplotypes, t.genetic_pos)
    w = 1.0 / t.support.astype(np.float64)
    beta = w / w.sum()
    p, q = toy.sharp_pair(t)
    return t, base, beta, p, q


def _dists(base, beta, p, q, eta):
    u_p = exact.utility_vector(base, p, beta)
    u_q = exact.utility_vector(base, q, beta)
    qp, lzp = exact.release_distribution(base, u_p, eta)
    qq, lzq = exact.release_distribution(base, u_q, eta)
    return u_p, u_q, qp, qq, lzp, lzq


@pytest.mark.parametrize("tau", TAUS)
def test_the_optimal_attacker_never_beats_p_succ_bound(tau, toy_pair):
    """Theorem 1, read as an attack: the best possible accuracy is (1+tau)/2.

    Exact, by enumeration, over the saturating pair — so this is the worst case and
    not a sample from it. A value above the bound would falsify the implementation.
    """
    t, base, beta, p, q = toy_pair
    eta = C.eta_tau(tau)
    *_, qp, qq, _, _ = (None, None, *_dists(base, beta, p, q, eta)[2:])
    adv = A.exact_advantage(qp, qq)
    bound = C.p_succ_bound(tau)
    assert adv <= bound + 1e-12, (
        f"the optimal attacker reaches {adv!r} against p_succ_bound {bound!r} at "
        f"tau = {tau}. The implementation is wrong, not the theorem.")
    assert adv >= 0.5 - 1e-12


@pytest.mark.parametrize("tau", TAUS)
def test_the_implemented_statistic_attains_the_exact_optimum(tau, toy_pair):
    """POSITIVE CONTROL (b), and the strongest statement available in this project.

    Not "the attacker failed" but "the attacker provably ATTAINS the optimum, and
    the optimum is below the bound". The sign rule on `Lambda` is Neyman-Pearson
    optimal, so its accuracy computed against the exact distributions must equal
    `(1 + TV)/2` — if it does not, the statistic is mis-derived.
    """
    t, base, beta, p, q = toy_pair
    eta = C.eta_tau(tau)
    u_p, u_q, qp, qq, lzp, lzq = _dists(base, beta, p, q, eta)
    lam = A.lrt_statistic(u_p, u_q, eta, lzp, lzq)
    # accuracy of the sign rule, weighted by the true distributions (no sampling)
    acc = 0.5 * (qp[lam > 0].sum() + 0.5 * qp[lam == 0].sum()
                 + qq[lam < 0].sum() + 0.5 * qq[lam == 0].sum())
    exact_opt = A.exact_advantage(qp, qq)
    assert abs(acc - exact_opt) < 1e-12, (
        f"the sign rule on Lambda gives {acc!r} but the Neyman-Pearson optimum is "
        f"{exact_opt!r} at tau = {tau}. The statistic is not the likelihood ratio.")


def test_at_tau_zero_the_optimal_attacker_is_exactly_at_chance(toy_pair):
    """The floor. At `tau` = 0 the release carries no target information at all, so
    the two hypotheses induce the SAME distribution and no attacker can do better
    than a coin. Exactly 0.5, not approximately."""
    t, base, beta, p, q = toy_pair
    u_p, u_q, qp, qq, lzp, lzq = _dists(base, beta, p, q, 0.0)
    assert np.allclose(qp, qq, atol=1e-15)
    assert abs(A.exact_advantage(qp, qq) - 0.5) < 1e-15
    lam = A.lrt_statistic(u_p, u_q, 0.0, lzp, lzq)
    assert np.allclose(lam, 0.0, atol=1e-15), (
        "Lambda is not identically zero at tau = 0, so the statistic is reading "
        "target information the release does not carry")


@pytest.mark.parametrize("tau", TAUS)
def test_POSITIVE_CONTROL_a_miscalibrated_mechanism_is_detected(tau, toy_pair):
    """POSITIVE CONTROL (a): a mechanism broken on purpose must RAISE the attacker's
    advantage, measured against its own correct-mechanism value.

    ⚠ THIS IS WHAT LICENSES REPORTING A NULL. If the attacker cannot tell a
    mechanism miscalibrated with `epsilon_tau` in place of `eta_tau` — the
    highest-severity error class in `CLAUDE.md`, and a one-character mistake — from
    the correct one, then its failure against the correct mechanism says nothing.

    ⚠⚠ AND THE CONTROL IS A PAIRED COMPARISON, NOT A THRESHOLD TEST, BECAUSE THE
    FIRST VERSION OF IT WAS WRONG. It required the doubled exponent to push the
    attacker past `p_succ_bound`, and it FAILED at every `tau`: at `tau` = 0.9 the
    doubled mechanism reaches only 0.8412 against a bound of 0.95. The reason is
    already on record for exact TV — the worst TV on this toy reaches only 27.4% to
    43.1% of `tau`, so the bound is never approached, and doubling a quantity that
    sits at 40% of its ceiling leaves it inside. **`p_succ_bound` therefore cannot
    serve as a detection threshold for a miscalibration, for exactly the reason
    `tv_empirical` cannot.** What CAN detect it is the attacker's own advantage
    moving, which is what this asserts.
    """
    t, base, beta, p, q = toy_pair
    eta, eps = C.eta_tau(tau), C.epsilon_tau(tau)
    bound = C.p_succ_bound(tau)
    *_, qp_ok, qq_ok, _, _ = (None, None, *_dists(base, beta, p, q, eta)[2:])
    *_, qp_bad, qq_bad, _, _ = (None, None, *_dists(base, beta, p, q, eps)[2:])
    adv_ok = A.exact_advantage(qp_ok, qq_ok)
    adv_bad = A.exact_advantage(qp_bad, qq_bad)

    # the correct mechanism is held under its bound; the broken one need not be
    assert adv_ok <= bound + 1e-12
    # and the break must MOVE the attacker, by a margin no numerical noise explains
    assert adv_bad > adv_ok + 1e-6, (
        f"a doubled tilt exponent moves the optimal attacker only from {adv_ok!r} "
        f"to {adv_bad!r} at tau = {tau}. The attacker cannot see a mechanism broken "
        f"on purpose, so its null against the correct one is uninterpretable.")
    # ⚠ and record what the threshold test would have claimed, so the hole is
    # visible rather than merely avoided
    assert not (adv_bad > bound), (
        f"at tau = {tau} the doubled exponent DID push the attacker past "
        f"p_succ_bound ({adv_bad!r} > {bound!r}). That contradicts the recorded "
        f"finding that TV reaches only 27-43% of tau on this toy — if this fires, "
        f"re-derive, because a threshold control would now be valid and this "
        f"docstring is stale.")


def test_POSITIVE_CONTROL_the_attacker_wins_against_an_untilted_release(toy_pair):
    """POSITIVE CONTROL, the other bracket: the instrument must be able to win.

    Hand the attacker a degenerate 'release' that is the target's own path with
    certainty. It must be right every time. If it cannot identify the target from
    the target's own data, the attacker is broken rather than the mechanism private
    — and that is indistinguishable from a privacy result unless this is checked.
    """
    t, base, beta, p, q = toy_pair
    eta = C.eta_tau(0.5)
    u_p, u_q, _, _, _, _ = _dists(base, beta, p, q, eta)
    # the output that IS p's own alleles, and the one that is q's
    i_p = int(np.argmax(u_p))
    i_q = int(np.argmax(u_q))
    assert i_p != i_q, "the sharp pair no longer separates, so this control is void"
    lam = A.lrt_statistic(u_p, u_q, eta, 0.0, 0.0)   # equal normalisers: pure LR
    acc = A.empirical_accuracy(np.array([lam[i_p]]), np.array([lam[i_q]]))
    assert acc == 1.0, (
        f"given each input's own most-favoured output, the attacker scores {acc!r} "
        f"rather than 1.0 — it cannot win even when the answer is handed to it")
