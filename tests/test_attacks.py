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

import ast
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


# ---------------------------------------------------------------------------
# ⚠ THE ADDENDUM MUST NOT BE ABLE TO IMPERSONATE THE PRE-REGISTRATION.
# `require_signed_prereg` globs `prereg_*.md` and returns the FIRST signed, clean
# match. A second file matching that pattern could therefore become the document an
# experiment stamps itself with -- if the original were renamed, unsigned or deleted,
# the lock would silently fall through to the addendum and report it as the
# pre-registration. So the addendum is deliberately NOT named `prereg_*`.
# ---------------------------------------------------------------------------

def test_the_decision_k_addendum_exists_and_is_signed() -> None:
    """Decision K post-dates the signature and changed the released object, so the
    record of that must exist and must itself be signed -- an unsigned addendum would
    leave the signed file's description of the chain standing uncorrected."""
    f = (Path(__file__).resolve().parents[1] / "paper"
         / "addendum_2026-10-11_decision_K.md")
    assert f.exists(), "the decision-K addendum is missing"
    text = f.read_text()
    assert "Status: SIGNED" in text
    assert "[OPEN]" not in text
    # the figures it exists to correct
    for n in ("22,134,562", "172,142", "4.5178e-08"):
        assert n in text, f"the addendum does not state {n}"
    assert "NOT NAMED" in text, (
        "the addendum must say WHY it is not named prereg_*, or the next person will "
        "name theirs that way")


# ---------------------------------------------------------------------------
# ⚠⚠ COMMON RANDOM NUMBERS, THE THIRD OCCURRENCE. `docs/bugs.md` carries "A headline
# finding measured under common random numbers" (2026-09-30) and `docs/memory.md`
# records seed reuse correlating two chromosomes of one release. The Arm 2 attack
# then derived its seed from the DRAW INDEX ALONE, so all 10 targets and all 7 `tau`
# reused 25 seeds. With a per-position tilt of 7.7e-07 the released path was
# bit-identical across `tau` at a fixed seed, the seven `tau` columns were one
# experiment reported seven times, and the naive standard error was ~2x too narrow.
#
# These tests exist so a fourth occurrence fails here rather than in a result table.
# ---------------------------------------------------------------------------

def _code_only(path: Path) -> str:
    """Source with COMMENTS AND DOCSTRINGS STRIPPED, for structural assertions.

    ⚠ Needed because the comments explaining a fix QUOTE the defective code they
    replaced -- which is exactly what you want in a comment and exactly what breaks a
    naive substring check. The first version of the two tests below failed on the
    very comments documenting the fix. `tokenize` drops comment tokens properly
    rather than guessing where a `#` is a comment and where it is inside a string.
    """
    import io
    import tokenize
    out = []
    prev_type = None
    for tok in tokenize.generate_tokens(io.StringIO(path.read_text()).readline):
        if tok.type == tokenize.COMMENT:
            continue
        # a bare STRING on its own logical line is a docstring
        if tok.type == tokenize.STRING and prev_type in (tokenize.INDENT,
                                                         tokenize.NEWLINE,
                                                         tokenize.NL, None):
            prev_type = tok.type
            continue
        out.append(tok.string)
        if tok.type not in (tokenize.NL, tokenize.NEWLINE):
            prev_type = tok.type
    # ⚠ WHITESPACE-FREE, and the needles below must be too. Joining tokens with
    # spaces turns `f(a, b)` into `f ( a , b )`, so a natural-looking substring
    # check fails on correct code -- which is what the first two versions of these
    # tests did. Normalising both sides removes the whole question.
    return "".join("".join(out).split())


def test_the_cell_seed_differs_across_target_tau_and_draw() -> None:
    """Every (target, tau, draw) cell must draw from its own stream."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "arm2", Path(__file__).resolve().parents[1] / "scripts" / "arm2_attack.py")
    # the derivation is a closure inside main(), so reproduce it exactly here and
    # pin the PROPERTY rather than the implementation
    import numpy as _np

    def cell(ti, taui, i, base=50_000):
        return int(_np.random.SeedSequence([base, ti, taui, i]).generate_state(1)[0])

    seeds = {cell(ti, taui, i)
             for ti in range(10) for taui in range(7) for i in range(25)}
    assert len(seeds) == 10 * 7 * 25, (
        f"{10*7*25 - len(seeds)} seed collision(s) across the grid. A collision means "
        f"two cells share a draw, which is the common-random-numbers failure this "
        f"test exists for.")
    # and it must be REPRODUCIBLE from the key
    assert cell(3, 4, 5) == cell(3, 4, 5)
    # changing ANY component must change the stream
    assert cell(3, 4, 5) != cell(4, 4, 5)
    assert cell(3, 4, 5) != cell(3, 5, 5)
    assert cell(3, 4, 5) != cell(3, 4, 6)


def test_the_attack_does_not_seed_from_the_draw_index_alone() -> None:
    """⚠ STRUCTURAL. The defect was literally `50_000 + i` in the job list. This
    fails if anyone reintroduces a seed that ignores the target or the tau.

    ⚠ IT ASSERTS THE THREE COMPONENTS, NOT ONE SPELLING. This pinned the literal
    text `_cell_seed(ti,taui,i)` until 2026-10-22, when the tau component became
    `SEED_GRID.index(tau)` so that a SUBSET of the grid keeps each point's own
    stream. That change preserved the property exactly and the test failed anyway,
    which is a test pinned to an incidental spelling rather than to the defect.
    """
    path = Path(__file__).resolve().parents[1] / "scripts" / "arm2_attack.py"
    assert "50_000+i" not in _code_only(path), (
        "the job list seeds from the draw index alone again -- all targets and all "
        "tau would share 25 streams (common random numbers, third occurrence)")
    # ⚠ The AST is built from the RAW file, not from `_code_only`, for two reasons:
    # `ast` drops comments and docstrings itself (so a comment quoting the old
    # defective call is invisible here, which is what `_code_only` exists for in a
    # substring check), and `_code_only` also strips the whitespace BETWEEN tokens,
    # so its output is not parseable Python at all.
    calls = [n for n in ast.walk(ast.parse(path.read_text()))
             if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "_cell_seed"]
    # ⚠ EVERY call site, not the first. The eta-driven mode added a second one on
    # 2026-10-23, and a gate that checked only `calls[0]` would have passed while the
    # new branch seeded however it liked.
    assert calls, "no _cell_seed call site at all; the per-cell seed derivation is gone"
    for k, call in enumerate(calls):
        assert len(call.args) == 3, (
            f"call site {k}: _cell_seed takes {len(call.args)} argument(s); the seed "
            f"must depend on the target index, the operating-point index AND the draw "
            f"index")
        names = {n.id for a in call.args for n in ast.walk(a) if isinstance(n, ast.Name)}
        for want, role in (("ti", "target index"), ("i", "draw index")):
            assert want in names, (
                f"call site {k}: the per-cell seed does not mention the {role} "
                f"({want!r}); every cell would share a stream over that axis "
                f"(common random numbers)")
        # ⚠⚠ ONE NAME PER MODE, and the list must grow when a mode does. The ceiling
        # mode added on 2026-10-28 seeds from `acc` (its operating point is the
        # `acc_ceiling` label, not a `tau` or an `eta`), and this whitelist was not
        # updated with it -- so the gate failed on a CORRECT call site and stayed red
        # from that commit until 2026-10-29. A gate that is red for a false reason is
        # as bad as one that cannot fire: it trains the next person to skip the suite.
        # ⚠ `acc` only, NOT `REPORTED_CEILINGS`. The real call site is
        # `_cell_seed(ti, CEIL_SEED_OFFSET + REPORTED_CEILINGS.index(acc), i)`, whose
        # varying part is `acc`. Admitting the container name would let
        # `REPORTED_CEILINGS.index(0.50)` -- a CONSTANT, so every ceiling column drawing
        # the same releases, which is the common-random-numbers defect itself -- sail
        # through the only gate that looks at these call sites.
        assert names & {"tau", "taui", "eta", "SEED_GRID", "ETA_GRID", "acc"}, (
            f"call site {k}: the seed does not mention the operating point; all "
            f"columns would reuse one set of streams, which is the defect measured "
            f"on 2026-10-10 -- the released path was bit-identical across tau at a "
            f"fixed seed")


def test_the_eta_mode_cannot_reuse_a_tau_cell_stream() -> None:
    """⚠⚠ The eta-driven points are a SECOND axis of operating points, and indexing
    them from 0 would have handed eta = 96 the stream of tau = 0. That is precisely
    occurrence four of this project's common-random-numbers defect
    (`docs/bugs.md`, 2026-10-22), so the offset is pinned here rather than trusted."""
    import importlib
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    M = importlib.import_module("arm2_attack")
    assert M.ETA_SEED_OFFSET >= len(M.SEED_GRID), (
        f"ETA_SEED_OFFSET {M.ETA_SEED_OFFSET} does not clear the tau grid "
        f"({len(M.SEED_GRID)} points), so an eta cell can collide with a tau cell")
    tau_idx = {M.SEED_GRID.index(t) for t in M.SEED_GRID}
    eta_idx = {M.ETA_SEED_OFFSET + M.ETA_GRID.index(e) for e in M.ETA_GRID}
    assert not (tau_idx & eta_idx), (
        f"the two operating-point axes share seed indices {sorted(tau_idx & eta_idx)}")
    # and the realised streams must differ, not merely the indices
    import numpy as np
    def cell(ti, pi, i):
        return int(np.random.SeedSequence([M.SEED_BASE, ti, pi, i]).generate_state(1)[0])
    seeds = {cell(3, p, 5) for p in sorted(tau_idx | eta_idx)}
    assert len(seeds) == len(tau_idx | eta_idx), "two operating points share a stream"


def test_decoy_selection_uses_a_stable_hash() -> None:
    """⚠ Python's str `hash` is randomised per interpreter by PYTHONHASHSEED, so a
    decoy chosen through it cannot be reproduced from the recorded seed -- measured,
    three invocations gave 1730677557, 783898805 and 1043054338 for 'HG00096'. A
    pre-registered experiment whose decoy set is irreproducible is not
    pre-registered in the sense standing constraint 6 means."""
    src = _code_only(Path(__file__).resolve().parents[1] / "scripts" / "arm2_attack.py")
    assert "hash(target_id)" not in src, (
        "decoy selection is back on Python's randomised str hash")
    assert "hashlib.sha256(target_id.encode())" in src, (
        "decoy selection must derive from a hash that is stable across interpreters")


def test_a_stable_hash_really_is_stable_across_processes() -> None:
    """Not a tautology: it runs a FRESH interpreter twice and compares, which is the
    only way to catch a per-process randomised hash."""
    import subprocess
    import sys as _sys
    prog = ("import hashlib;"
            "print(int.from_bytes(hashlib.sha256(b'HG00096').digest()[:4],'big'))")
    outs = {subprocess.run([_sys.executable, "-I", "-c", prog],
                           capture_output=True, text=True).stdout.strip()
            for _ in range(3)}
    assert len(outs) == 1, f"the 'stable' hash varied across processes: {outs}"
