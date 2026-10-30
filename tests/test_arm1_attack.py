"""GATES on `scripts/arm1_attack.py`, each injected with the error it exists to catch.

The spec for this driver (`docs/memory.md` 2026-10-27, B-ATK-arm1-spec) lists four
silent-wrong-answer paths. Three of them are things a reader cannot see in an output
table, so each gets a gate here that is shown to FIRE:

  1. `arm2_attack_a3a4._A3` is an UNKEYED module global. Arm 1's bi-allelic axis is per
     target and `t_chain` spans only 332,516..332,521, so many targets SHARE a T -- and
     where they do, a stale cache is served silently and every A3 number for that target
     is computed on another target's site axis.
  2. TWO POSINGS, ONLY ONE BOUNDED. The FIXED MECHANISM posing (one mechanism, two input
     paths) is what Theorem 1 bounds. Drawing `Q_q` from q's own leave-one-out panel is
     two mechanisms and is NOT bounded, so an accuracy from it compared to `acc_ceiling`
     would report a violation with the theorem intact. The gate is STRUCTURAL, because
     the defect would be an added call rather than a wrong value.
  3. The fitted `eta` moves on a refit, so a seed keyed on its VALUE would silently
     redraw every cell and two targets' "same" cell would not be the same grid point.

Plus the arithmetic that decides whether a null is a null at all: at `eta` = 0 the two
hypotheses are the SAME distribution, so EVERY statistic must score exactly 0.5 balanced
-- while TPR and TNR may sit at 1 and 0, which is the constant-answer artefact the
DEGENERATE flag exists to name.
"""

from __future__ import annotations

import ast
import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src"), str(REPO / "scripts")]

import arm1_attack as AA            # noqa: E402
import arm2_attack_a3a4 as A3M      # noqa: E402


# ---------------------------------------------------------------- 1. the _A3 cache
def _toy_panel(t: int, seed: int) -> np.ndarray:
    """A (K, T) cohort panel whose every column is bi-allelic, so `a3_sites` keeps it."""
    rng = np.random.default_rng(seed)
    hap = rng.integers(0, 2, size=(8, t)).astype(np.int16)
    hap[0, :] = 0          # both codes present at every column
    hap[1, :] = 1
    return hap


def test_a3_cache_RAISES_on_an_axis_this_driver_did_not_stamp():
    """An unstamped non-empty `_A3` could belong to any target, so it must RAISE.

    ⚠ THIS IS THE GATE THE FIRST VERSION DID NOT HAVE. `_a3_for` originally cleared the
    cache unconditionally, which made its length assertion unfalsifiable -- nothing can
    be stale if the cache is always empty. This test is what caught that, and it is the
    reason the cache is keyed instead of cleared.
    """
    A3M._A3.clear()
    A3M.a3_sites({"hap": _toy_panel(300, 1)})          # another caller, no stamp
    with pytest.raises(RuntimeError, match="did not stamp"):
        AA._a3_for({"target": "T", "panel": _toy_panel(300, 2), "t_chain": 300})
    A3M._A3.clear()


def test_a3_cache_recomputes_for_a_SAME_LENGTH_other_target():
    """⚠⚠ THE HAZARD THAT MATTERS: equal `t_chain`, different axis.

    Arm 1's targets mostly share a `t_chain`, so a length check CANNOT distinguish this
    target's axis from another's. The key can. Poison with a properly stamped axis from
    panel A, then ask for panel B of the SAME length, and the served frequencies must be
    B's. This test fails if the key comparison is ever dropped.
    """
    a, b = _toy_panel(300, 3), _toy_panel(300, 4)
    A3M._A3.clear()
    first = AA._a3_for({"target": "A", "panel": a, "t_chain": 300})
    a_af = first["alt_af"].copy()
    # ⚠ `_a3_for` returns the LIVE module global, not a copy -- clearing `_A3` empties
    # the dict the caller is holding. Copy what is needed before any clear. (The driver
    # is safe: `run_cell` holds `sites` for one target and nothing clears mid-target.)
    got_af = AA._a3_for({"target": "B", "panel": b, "t_chain": 300})["alt_af"].copy()
    assert got_af.shape == a_af.shape, "same-length axes, as designed"
    assert not np.allclose(got_af, a_af), (
        "_a3_for served A's axis for B. Both chains are 300 long, so only the key can "
        "catch this -- a length check passes on the wrong axis.")
    A3M._A3.clear()
    assert np.allclose(got_af, A3M.a3_sites({"hap": b})["alt_af"]), \
        "the served axis is not B's own either"
    A3M._A3.clear()


def test_a3_cache_is_served_when_the_panel_is_unchanged():
    """It is still a cache: the same panel twice must not recompute."""
    a = _toy_panel(300, 5)
    A3M._A3.clear()
    one = AA._a3_for({"target": "A", "panel": a, "t_chain": 300})
    two = AA._a3_for({"target": "A", "panel": a, "t_chain": 300})
    assert one is two, "the keyed cache stopped serving; every target would recompute"
    A3M._A3.clear()


# ------------------------------------------------- 2. the posing, structurally
def _calls_in(fn_name: str, names: tuple) -> list:
    """Every call to one of `names` inside `fn_name` of the driver."""
    tree = ast.parse((REPO / "scripts" / "arm1_attack.py").read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == fn_name)
    out = []
    for n in ast.walk(fn):
        if not isinstance(n, ast.Call):
            continue
        f = n.func
        nm = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
        if nm in names:
            out.append((nm, n.lineno, {k.arg for k in n.keywords}))
    return out


def test_run_cell_builds_no_second_mechanism():
    """FIXED MECHANISM posing: `run_cell` may not construct a model or a panel for q.

    ⚠ A value test cannot catch this. The PROCEDURE posing returns perfectly plausible
    accuracies -- it is simply not the quantity `acc_ceiling` bounds -- so the gate has
    to be on the SHAPE of the code.
    """
    bad = _calls_in("run_cell", ("from_arrays", "loo_chain", "masked_haplotypes",
                                 "baseline_log_norm", "beta_v", "target_state"))
    allowed = {"target_state"}        # p's own state, built once
    offenders = [c for c in bad if c[0] not in allowed]
    assert not offenders, (
        f"run_cell builds mechanism state at {offenders}. Theorem 1 bounds ONE "
        f"mechanism fed TWO input paths; a second panel or model is the PROCEDURE "
        f"posing, which is not bounded and whose accuracy would read as a violation.")


def test_both_releases_in_run_cell_use_the_same_model_and_panel():
    """Both `ffbs.release` calls must pass p's `st['model']` and `st['panel']`."""
    src = (REPO / "scripts" / "arm1_attack.py").read_text()
    tree = ast.parse(src)
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "run_cell")
    rels = [n for n in ast.walk(fn) if isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute) and n.func.attr == "release"]
    assert len(rels) == 2, f"expected exactly two releases in run_cell, found {len(rels)}"
    for r in rels:
        first_two = [ast.unparse(a) for a in r.args[:2]]
        assert first_two == ["st['model']", "st['panel']"], (
            f"a release in run_cell is driven by {first_two}, not p's own model and "
            f"panel -- that is the unbounded PROCEDURE posing.")


def test_driver_never_calls_arm2_load_state():
    """`load_state` caches the ARM 2 full-cohort panel (K = 88, share 1.0)."""
    src = (REPO / "scripts" / "arm1_attack.py").read_text()
    tree = ast.parse(src)
    hits = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Call)
            and getattr(n.func, "id", getattr(n.func, "attr", None)) == "load_state"]
    assert not hits, (
        f"arm1_attack.py calls load_state() at line(s) {hits}. That serves K = 88 at "
        f"share 1.0 to an Arm 1 cell and measures the wrong mechanism silently.")


# ------------------------------------------------------------- 3. the seeding
def test_cell_seed_is_keyed_on_the_label_index_not_the_eta():
    a = AA.cell_seed("HG00438", "0.55", 0, "D1_extremal", "p")
    b = AA.cell_seed("HG00438", "0.80", 0, "D1_extremal", "p")
    assert a != b, "two acc_ceiling labels share a seed -- common random numbers"
    assert a == AA.cell_seed("HG00438", "0.55", 0, "D1_extremal", "p"), "not reproducible"


def test_cell_seed_separates_hypotheses_conditions_draws_and_targets():
    base = ("HG00438", "0.80", 0, "D1_extremal", "p")
    assert AA.cell_seed(*base) != AA.cell_seed("HG00438", "0.80", 0, "D1_extremal", "q"), \
        "the two hypotheses share a seed, which correlates the arms of the test"
    assert AA.cell_seed(*base) != AA.cell_seed("HG00438", "0.80", 1, "D1_extremal", "p")
    assert AA.cell_seed(*base) != AA.cell_seed("HG00438", "0.80", 0, "D2_nearest", "p")
    assert AA.cell_seed(*base) != AA.cell_seed("HG00621", "0.80", 0, "D1_extremal", "p")


def test_cell_seed_refuses_a_label_with_no_fixed_place():
    with pytest.raises(ValueError, match="no fixed place"):
        AA.cell_seed("HG00438", "0.42", 0, "D1_extremal", "p")


# --------------------------------------------- 4. what a null is, arithmetically
def test_at_eta_zero_every_statistic_scores_exactly_one_half():
    """⚠ THE 0.50 COLUMN IS A HARNESS CHECK AND THIS IS WHY IT CANNOT FAIL HONESTLY.

    At `eta` = 0, `Q_p = Q_q = R_D`, so y has the same law under both hypotheses and for
    ANY statistic T: balanced = 0.5*(P(T>0) + 0.5P(T=0)) + 0.5*(P(T<0) + 0.5P(T=0))
    = 0.5 exactly. A 0.50 column away from 0.5 is a harness bug, never a result.
    """
    rng = np.random.default_rng(0)
    for _ in range(20):
        x = rng.normal(size=400)          # one sample, split into the two arms
        bal, se, tpr, tnr = AA._acc(x[:200].tolist(), x[200:].tolist())
        assert abs(0.5 * tpr + 0.5 * tnr - bal) < 1e-12
    # and the exact-tie case the mechanism actually produces at eta = 0
    bal, se, tpr, tnr = AA._acc([0.0] * 50, [0.0] * 50)
    assert bal == 0.5 and tpr == 0.5 and tnr == 0.5
    # ⚠ AND THE SE IS NOT ZERO HERE, WHICH IS WORTH KNOWING. Every Lambda is exactly 0,
    # so the accuracy is DETERMINISTICALLY 0.5 and its true variance is zero -- but the
    # binomial SE is computed from TPR = TNR = 0.5, where it is at its MAXIMUM:
    # 0.5*sqrt(2 * 0.25/50) = 0.05. So the 0.50 column's SE is an over-estimate, not an
    # uncertainty. It changes no verdict (the excess there is exactly 0, so z = 0 either
    # way) but the figure must not be read as sampling noise on a tied cell.
    assert math.isclose(se, 0.5 * math.sqrt(2 * 0.25 / 50))
    assert math.isclose(se, 0.05)


def test_degenerate_is_detected_where_balanced_accuracy_looks_clean():
    """A constant answer scores 0.5 balanced -- the clean-null face of a dead test."""
    lp = [1.0] * 100          # always "p"
    lq = [1.0] * 100          # ... under both hypotheses
    bal, se, tpr, tnr = AA._acc(lp, lq)
    assert bal == 0.5, "a constant answer scores exactly chance balanced"
    assert tpr == 1.0 and tnr == 0.0
    assert min(tpr, tnr) < AA.DEGENERATE_FLOOR, (
        "the DEGENERATE floor must catch this; otherwise it is reported as a null")


def test_acc_drift_assertion_fires():
    """`_acc` asserts its own two definitions agree; show the assertion is live."""
    real = AA.A.empirical_accuracy
    try:
        AA.A.empirical_accuracy = lambda a, b: 0.123
        with pytest.raises(AssertionError, match="definitions drifted"):
            AA._acc([1.0, -1.0], [-1.0, 1.0])
    finally:
        AA.A.empirical_accuracy = real


# ------------------------------------------------------------ 5. the ceiling
def test_ceiling_from_kl_is_pinsker_and_clamps():
    tv, acc = AA.ceiling_from_kl(0.72)
    assert math.isclose(tv, math.sqrt(0.72 / 2.0))
    assert math.isclose(acc, 0.5 * (1.0 + tv))
    assert AA.ceiling_from_kl(0.0) == (0.0, 0.5), "eta = 0 must give exactly 0.5"
    assert AA.ceiling_from_kl(1e9) == (1.0, 1.0), "TV must clamp at 1"
    # monotone, so a larger KL can never report a tighter bound
    xs = [AA.ceiling_from_kl(k)[1] for k in (0.02, 0.08, 0.18, 0.32, 0.5, 0.72)]
    assert xs == sorted(xs)


def test_lambda_sign_convention_positive_means_p():
    """Lambda = eta*(u_p - u_q) - (log Z_p - log Z_q), and `d_logz` IS that difference.

    The driver passes `log_z_p = 0.0` and `log_z_q = -d_logz`, which must reproduce the
    statistic exactly -- an off-by-sign here would invert TPR and TNR and the table would
    still look plausible.
    """
    eta, d_logz = 3.0, 1.25
    got = float(AA.A.lrt_statistic(np.array([0.9]), np.array([0.2]), eta,
                                   0.0, -d_logz)[0])
    assert math.isclose(got, eta * (0.9 - 0.2) - d_logz)
    assert got > 0, "p carries the higher utility here, so Lambda must favour p"


def test_extremal_decoy_is_the_shared_one_not_a_copy():
    """The fitted `eta` was solved against `exact_kl_stable.extremal`; a second copy here
    would let the decoy the ceiling was computed for drift from the one scored."""
    import exact_kl_stable as EKS
    assert AA.extremal is EKS.extremal
    p = np.array([[0, 1], [1, 0], [-1, 2]], dtype=np.int16)
    e = AA.extremal(p)
    called = p >= 0
    assert (e[called] != p[called]).all(), "every called allele must differ"
    assert (e[~called] == p[~called]).all(), "uncalled cells must stay uncalled"


# -------------------------------------------- 5. D3_random must never BE D2_nearest
# ⚠ WHY THIS GATE EXISTS. `pick_decoys` drew D3 uniformly from the SAME pool it took
# D2 from, so a draw could land on the nearest donor and did: Arm 2's HG01323 and Arm
# 1's HG00741 and HG01123 each had D3_random == D2_nearest in the shipped runs. Those
# targets contributed one genome to two conditions reported as independent -- 1/10 of
# Arm 2's D3 row and 2/44 of Arm 1's. Nothing raised, because a collision is a legal
# uniform draw. Both arms now REJECT the nearest, and these gates fire on a revert.
class _ScriptedRNG:
    """An rng whose draws are pinned, so a collision is forced rather than hoped for."""

    def __init__(self, seq):
        self.seq = list(seq)
        self.calls = 0

    def integers(self, n):
        self.calls += 1
        assert self.seq, "pick_decoys drew more times than the script allows"
        return self.seq.pop(0)


def _diploid(val: int, t: int = 6) -> np.ndarray:
    return np.full((t, 2), val, dtype=np.int16)


def _arm2_state():
    """Target plus three donors; the FIRST listed donor is the concordance maximum."""
    return {"paths": {"T0": _diploid(0), "NEAR": _diploid(0),
                      "FAR1": _diploid(1), "FAR2": _diploid(1)}}


def test_arm2_d3_rejects_the_nearest_when_the_draw_lands_on_it():
    import arm2_attack as A2
    st = _arm2_state()
    # others == [NEAR, FAR1, FAR2]; index 0 IS the nearest, so draw 0 then 1.
    rng = _ScriptedRNG([0, 1])
    d = A2.pick_decoys("T0", st, rng)
    assert d["D2_nearest"]["id"] == "NEAR"
    assert d["D3_random"]["id"] != d["D2_nearest"]["id"], (
        "D3_random came back AS D2_nearest -- one genome in two conditions")
    assert d["D3_random"]["id"] == "FAR1"


def test_arm2_d3_is_untouched_when_the_draw_does_not_collide():
    """The fix must be REJECTION, not a narrowed pool: a clean draw is left alone."""
    import arm2_attack as A2
    rng = _ScriptedRNG([2])
    d = A2.pick_decoys("T0", _arm2_state(), rng)
    assert d["D3_random"]["id"] == "FAR2"
    assert rng.calls == 1, (
        f"a non-colliding draw consumed {rng.calls} draws, so the pool was re-indexed "
        f"and every other target's D3 identity would move too")


def test_arm2_d3_raises_when_the_pool_cannot_satisfy_the_invariant():
    import arm2_attack as A2
    st = {"paths": {"T0": _diploid(0), "ONLY": _diploid(0)}}
    with pytest.raises(ValueError, match="cannot differ from D2_nearest"):
        A2.pick_decoys("T0", st, _ScriptedRNG([0]))


def test_arm1_d3_rejects_the_nearest_when_the_draw_lands_on_it():
    # sorted(others) == [A_NEAR, B_FAR, C_FAR]; index 0 IS the nearest.
    st = {"path": _diploid(0),
          "others": {"A_NEAR": _diploid(0), "B_FAR": _diploid(1), "C_FAR": _diploid(1)}}
    rng = _ScriptedRNG([0, 2])
    d = AA.pick_decoys(st, rng)
    assert d["D2_nearest"]["id"] == "A_NEAR"
    assert d["D3_random"]["id"] != d["D2_nearest"]["id"]
    assert d["D3_random"]["id"] == "C_FAR"


def test_arm1_d3_raises_when_the_pool_cannot_satisfy_the_invariant():
    st = {"path": _diploid(0), "others": {"ONLY": _diploid(0)}}
    with pytest.raises(ValueError, match="cannot differ"):
        AA.pick_decoys(st, _ScriptedRNG([0]))


def test_the_recorded_arm2_collision_is_the_one_this_gate_reproduces():
    """HG01323's REAL seed drew index 5 of 9, which was its own nearest donor.

    Not a synthetic case: the seed and the index are the shipped run's.
    """
    import hashlib
    import arm2_attack as A2
    stable = int.from_bytes(hashlib.sha256(b"HG01323").digest()[:4], "big")
    rng = np.random.default_rng([A2.DECOY_SEED, stable])
    assert int(rng.integers(9)) == 5, "the historical draw is not reproducing"
    # Rebuild a 9-donor pool whose index 5 is the concordance maximum, as it was.
    names = [f"D{i}" for i in range(9)]
    names[5] = "HG04211"
    paths = {"HG01323": _diploid(0)}
    for n in names:
        paths[n] = _diploid(0) if n == "HG04211" else _diploid(1)
    d = A2.pick_decoys("HG01323", {"paths": paths},
                       np.random.default_rng([A2.DECOY_SEED, stable]))
    assert d["D2_nearest"]["id"] == "HG04211"
    assert d["D3_random"]["id"] != "HG04211", (
        "the exact collision the shipped run had is still reachable")


# ------------------------- 6. the decoy is never the TARGET, in either arm
# ⚠ Verified empirically first (0 cells in either shipped run), then gated, because
# "it held last time" is not an invariant. Both arms exclude the target with an
# explicit guard in a comprehension -- `arm2_attack.py` inside `pick_decoys` itself,
# `arm1_attack.py` up in `load_state` -- and a dropped guard would make the target its
# OWN nearest decoy, since concordance(p, p) = 1.0 is the maximum possible. The
# measured accuracy would then collapse to chance and read as a flawless null.
def test_arm2_rejects_the_target_as_its_own_decoy():
    import arm2_attack as A2
    st = _arm2_state()
    d = A2.pick_decoys("T0", st, _ScriptedRNG([1, 2]))
    assert d["D2_nearest"]["id"] != "T0"
    assert d["D3_random"]["id"] != "T0"


def test_arm2_would_pick_the_target_without_its_guard():
    """The mutation this gate exists for, shown rather than asserted in a comment.

    With the target left in the pool, concordance(p, p) = 1.0 beats every real donor,
    so D2_nearest becomes the target itself.
    """
    import arm2_attack as A2
    st = _arm2_state()
    p = st["paths"]["T0"]
    pool_with_target = list(st["paths"])                      # the guard removed
    conc = {k: A2.concordance(p, st["paths"][k]) for k in pool_with_target}
    assert max(conc, key=conc.get) == "T0", (
        "if this is not the target then concordance is not maximal at identity and "
        "the guard this gate protects is not the thing that keeps the target out")


def test_both_arms_carry_the_target_exclusion_guard_in_source():
    """⚠ STRUCTURAL, and deliberately so. Arm 1's guard is in `load_state`, which needs
    the cohort arrays, so a behavioural gate on it would be skipif'd on a gitignored
    artifact -- the silent-coverage-hole defect `docs/bugs.md` 2026-10-29 records. The
    guard is one clause; assert the clause is there.
    """
    a2 = (REPO / "scripts" / "arm2_attack.py").read_text(encoding="utf-8")
    a1 = (REPO / "scripts" / "arm1_attack.py").read_text(encoding="utf-8")
    assert 'others = [k for k in st["paths"] if k != target_id]' in a2, (
        "arm2_attack.pick_decoys no longer excludes the target from the decoy pool")
    assert "for s, r in index.items() if s != target" in a1, (
        "arm1_attack.load_state no longer excludes the target from the decoy pool")
    # And the pool feeds BOTH conditions, so one guard covers D2 and D3 in each arm.
    assert "conc = {k: concordance(p, st[\"paths\"][k]) for k in others}" in a2
    assert 'conc = {s: _concordance(p, q) for s, q in st["others"].items()}' in a1
