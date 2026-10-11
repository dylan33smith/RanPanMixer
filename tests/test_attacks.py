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


# ---------------------------------------------------------------------------
# The pre-registration lock
# ---------------------------------------------------------------------------

def test_an_unsigned_prereg_blocks_an_attack(tmp_path):
    """⚠ Constraint 6, enforced in code rather than remembered.

    What the constraint guards against is the one failure a person cannot
    self-police: having seen a result, you cannot PROVE you did not pick the
    analysis to suit it. So the guard has to be something other than good
    intentions.
    """
    (tmp_path / "prereg_2026-10-10.md").write_text(
        "**Status: ⚠ ALL EIGHT FIELDS FILLED — AWAITING SIGNATURE.**\n")
    with pytest.raises(RuntimeError, match="NOT SIGNED"):
        A.require_signed_prereg(str(tmp_path))


def test_the_refusal_names_what_is_unfilled(tmp_path):
    """The useful response to the refusal is to go and fill the fields, so it has
    to say which ones are missing."""
    (tmp_path / "prereg_x.md").write_text(
        "**Status: DRAFT**\n\n[OPEN] the tau grid\n\n- [ ] 1. something\n- [ ] 2. else\n")
    with pytest.raises(RuntimeError) as e:
        A.require_signed_prereg(str(tmp_path))
    msg = str(e.value)
    assert "[OPEN] marker" in msg and "unchecked checklist item" in msg
    assert "deliberate act" in msg


def test_a_signed_prereg_allows_it_and_hands_back_its_identity(tmp_path):
    """A result must be stampable with the exact pre-registration it ran under.
    Two experiments under different pre-registrations are not comparable, and the
    stamp is how anyone finds that out."""
    (tmp_path / "prereg_2026-10-10.md").write_text(
        "**Status: SIGNED** git_rev `deadbee`\n\n- [x] 1. done\n")
    name, rev = A.require_signed_prereg(str(tmp_path))
    assert name == "prereg_2026-10-10.md" and rev == "deadbee"


def test_a_missing_prereg_is_refused_too(tmp_path):
    """No file at all must not read as 'nothing to check'."""
    with pytest.raises(RuntimeError, match="no pre-registration found"):
        A.require_signed_prereg(str(tmp_path))


def test_the_real_prereg_is_signed_and_opens_the_lock():
    """⚠ A LIVE CHECK ON THIS REPOSITORY, not a fixture. ✅ THE EXPECTATION FLIPPED ON
    2026-10-10, in the signing commit, which is what the previous version of this test
    demanded: it asserted the lock was HOLDING and said "when it is signed this test's
    expectation flips, and that flip should be a deliberate edit made in the signing
    commit -- not something discovered later." This is that edit.

    What it now guards is the opposite failure, and it is the one that matters from
    here: a signed pre-registration being silently UNSIGNED, reverted, renamed or
    emptied, which would let an attack run against a design nobody committed to. The
    `git_rev` is asserted non-empty because a signature with no code state attached
    cannot tell anyone WHAT was pre-registered.
    """
    name, rev = A.require_signed_prereg()          # must NOT raise
    assert name == "prereg_2026-10-10.md", name
    assert rev and len(rev) >= 7 and all(c in "0123456789abcdef" for c in rev), (
        f"the signature records no usable git_rev ({rev!r}); a result stamped with a "
        f"pre-registration that names no code state is not traceable")


def test_the_signed_prereg_has_no_open_fields_left():
    """⚠ Belt and braces on the lock's own crudeness. `require_signed_prereg` refuses
    on any `[OPEN]` marker or unticked box, so this asserts the file genuinely has
    none rather than trusting that the refusal logic is the only reader. All eight
    checklist items must be ticked -- seven filled plus field 6, reopened and re-closed
    on 2026-10-10 at the balanced 25-draw design.
    """
    text = (Path(__file__).resolve().parents[1] / "paper" / "prereg_2026-10-10.md").read_text()
    assert "[OPEN]" not in text
    assert "\n- [ ] " not in text
    assert text.count("\n- [x] ") == 8, (
        f"expected all 8 checklist items ticked, found {text.count(chr(10) + '- [x] ')}")
    assert "25 draws" in text and "7,000 releases" in text, (
        "the signed design must state the pinned draw count and release total")


# ---------------------------------------------------------------------------
# ⚠⚠ THE ONE-SIDED TRAP. `p_succ_bound` = (1 + tau)/2 bounds the EQUAL-PRIOR
# accuracy and NOTHING ELSE. On 2026-10-10 `scripts/arm2_attack.py` drew every
# release from `Q_p`, computed the one-sided rate P_p(Lambda > 0), and compared THAT
# to the bound -- a comparison that fires with the theorem fully intact. The suite
# had 265 passing tests and could not see it, because `empirical_accuracy` was
# correct and the script had reimplemented it inline.
# ---------------------------------------------------------------------------

def test_empirical_accuracy_is_balanced_not_one_sided() -> None:
    """A degenerate attacker that ALWAYS answers "p": one-sided 1.0, balanced 0.5."""
    always_p_under_p = np.ones(100)      # Lambda > 0 on every draw from p
    always_p_under_q = np.ones(100)      # ...and on every draw from q too
    one_sided = float((always_p_under_p > 0).mean())
    assert one_sided == 1.0
    bal = A.empirical_accuracy(always_p_under_p, always_p_under_q)
    assert bal == 0.5, (
        f"a constant 'always p' answer must score exactly chance when balanced, got "
        f"{bal}. If this returns the one-sided rate instead, every attack figure in "
        f"the paper is compared to a bound that does not apply to it.")


def test_one_sided_rate_can_exceed_the_bound_without_violating_it() -> None:
    """⚠ The derivation, as a test: balanced = 0.5*TPR + 0.5*(1-FPR) <= (1+TV)/2
    implies TPR <= TV + FPR, so a theorem-respecting attacker with FPR = 0.5 reaches
    TPR = tau + 0.5 > (1+tau)/2 for every tau > 0.

    This test exists so that anyone reinstating a one-sided comparison has to delete
    an assertion that spells out why it is wrong."""
    for tau in (0.05, 0.10, 0.25, 0.50, 0.90):
        bound = C.p_succ_bound(tau)
        n = 1000
        # ⚠ capped at n: TPR = tau + 0.5 exceeds 1 for tau > 0.5, where the
        # constraint TPR <= TV + FPR is no longer the binding one.
        k_tpr = min(n, int(round((tau + 0.5) * n)))     # TPR = min(1, tau + 0.5)
        under_p = np.concatenate([np.ones(k_tpr), -np.ones(n - k_tpr)])
        under_q = np.concatenate([np.ones(n // 2), -np.ones(n - n // 2)])  # FPR = 0.5
        one_sided = float((under_p > 0).mean())
        bal = A.empirical_accuracy(under_p, under_q)
        assert one_sided > bound + 1e-9, (
            f"tau={tau}: the one-sided rate {one_sided:.4f} should exceed the bound "
            f"{bound:.4f} -- that is the whole point of this test")
        assert bal <= bound + 1e-9, (
            f"tau={tau}: the BALANCED accuracy {bal:.4f} must respect the bound "
            f"{bound:.4f}; if it does not, the constructed example is wrong, not the "
            f"theorem")
