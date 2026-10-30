#!/usr/bin/env python
"""A1 (the optimal distinguisher) and A3 on ARM 1's chr21 `acc_ceiling` grid.

Spec and every hazard below: `docs/memory.md` 2026-10-27 (B-ATK-arm1-spec). That entry
was written from an instrument check BEFORE any powered run, and this script implements
it. Read it before changing anything here.

⚠⚠ `acc_ceiling` IS A PER-PAIR QUANTITY AND THE DIRECTORY LABEL IS THE D1_extremal
CEILING ALONE. `kl_fit_verified.json` was solved with `q = extremal(p)`, so the label
0.80 is D1's ceiling; HG00438's D2_nearest ceiling at that label is 0.678293 and D3's is
0.696243. Comparing a D2 or D3 accuracy to the column header would test it against a
bound ~0.12 TOO LOOSE. So this script computes `acc_ceiling` PER (target, label,
condition) from `kl_stable` and prints it as a ROW. It also asserts per cell that D1 is
the max-TV decoy, which makes the label's conservatism measured rather than assumed.

⚠⚠ TWO POSINGS, ONLY ONE BOUNDED, AND THIS IMPLEMENTS (i).
  (i)  FIXED MECHANISM -- ONE mechanism (p's leave-one-out panel, chain, beta and eta)
       fed TWO input paths. Theorem 1 bounds exactly this, and `acc_ceiling` is computed
       for exactly this. Every `ffbs.release` below therefore passes p's `m`, `panel` and
       `util` whichever path is the input.
  (ii) PROCEDURE POSING -- `Q_q` drawn from q's OWN leave-one-out panel, which is what
       the real release procedure would do. That is TWO mechanisms and is NOT bounded by
       Theorem 1; an accuracy from (ii) compared to `acc_ceiling` would report a
       violation with the theorem fully intact. NOT implemented, deliberately.

⚠⚠ `calibration.p_succ_bound` RAISES ON EVERY CELL HERE, so neither existing `_report`
is reusable. Ceiling-mode releases record `tau` = 1.0 because `tanh(25494)` saturates in
float64, and `p_succ_bound` requires `tau` in [0, 1). Worse, `arm2_attack_a3a4._report`
does `if tau >= 1.0: continue`, which would skip EVERY cell and then print "no condition
exceeded p_succ_bound" having checked nothing. This script compares against the per-pair
Pinsker ceiling it computes itself and never imports either reporter.

⚠⚠ NEVER CALL `arm2_attack.load_state()`. It caches the ARM 2 full-cohort panel (K = 88,
share 1.0) in a module global, and serving that to an Arm 1 cell would silently measure
the wrong mechanism. Arm 1 state is per target -- K = 86 and the chain both move -- so it
is built here, per target, and nothing is cached across targets.

⚠⚠ `a3_sites` CACHES IN THE UNKEYED MODULE GLOBAL `_A3`. Arm 1's bi-allelic axis is per
target and measured `t_chain` spans only 332,516..332,521, so MANY TARGETS SHARE A T --
where they do, a stale cache is served silently and every A3 number for that target is
computed on another target's site axis. `_a3_for` clears the cache and then asserts the
returned mask length against this target's own chain, so a stale serve raises.

⚠ THE PANEL ASYMMETRY UNBALANCES THE TEST, AND IT IS NOT A BUG. The target is HELD OUT
while every D2/D3 decoy is IN the panel, so the decoy is exactly reproducible by the
mosaic and the target is not. That offset is a legitimate part of the Neyman-Pearson
statistic and sits inside the Theorem 1 bounds, but it drives TPR and TNR apart. So TPR
and TNR are printed as their OWN rows, and a row is auto-labelled DEGENERATE when
min(TPR, TNR) < 0.05 -- an unbalanced test landing near 0.5 is the constant-answer
artefact of `docs/bugs.md` wearing a clean-null face.

⚠ THE `acc_ceiling` 0.50 COLUMN IS A HARNESS CHECK, NOT A RESULT. The fit has six labels
(0.55..0.80); 0.50 is `eta` = 0, where `Q_p = Q_q = R_D` identically, so the true
balanced accuracy is EXACTLY 0.5 for every attacker and every condition. It cannot be a
violation. If that column is not 0.5 to within noise, the harness is wrong.

    python scripts/arm1_attack.py --draws 5 --targets 2          # smoke
    python scripts/arm1_attack.py --draws 100 --workers 22       # powered
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import statistics as stats
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src"), str(REPO / "scripts")]

import numpy as np                                             # noqa: E402

import attacks as A                                            # noqa: E402
from cohort import loo                                         # noqa: E402
from mechanism import ffbs, hmm, utility                       # noqa: E402

import arm2_attack_a3a4 as A3M                                 # noqa: E402
# ⚠ the acc-line arithmetic has ONE home, scripts/arm2_attack.py; a second copy
# is how the ceiling branch reached one driver and not the other.
from arm2_attack import _interval                              # noqa: E402
from exact_kl_stable import extremal, kl_stable                # noqa: E402

CHROM = "chr21"
CEILINGS = ("0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80")
CONDS = ("D1_extremal", "D2_nearest", "D3_random")
FIT = REPO / "logs" / "B_eta_arm1_chr21" / "kl_fit_verified.json"
SURVEY = REPO / "logs" / "B_loo_chain_survey.json"
WEIGHTS = "uniform"

# Disjoint from A1's 50_000, A3/A4's 70_000 and the power run's 90_000.
SEED_BASE = 80_000
DECOY_SEED = 20261028
# min(TPR, TNR) below this and the row is a constant answer, not a null.
DEGENERATE_FLOOR = 0.05

from arm1_chr21_kl import loo_share_for                        # noqa: E402


def _h(s: str) -> int:
    """sha256, never `hash()`: Python's str hash is PYTHONHASHSEED-salted and has
    already made a decoy choice irreproducible once in this project."""
    return int.from_bytes(hashlib.sha256(s.encode()).digest()[:4], "big")


def cell_seed(target: str, label: str, draw: int, cond: str | None, hyp: str) -> int:
    """One seed per (target, label, draw, condition, hypothesis).

    ⚠ KEYED ON THE LABEL'S INDEX IN `CEILINGS`, never on the fitted `eta`: Arm 1's
    `eta` is per target and moves on a refit (it already moved once, between
    `kl_fit.json` and `kl_fit_verified.json`), so keying on the value would redraw every
    cell silently and two targets' "same" cell would not be the same grid point.
    ⚠ And keyed on the HYPOTHESIS: reusing the Q_p seed for the Q_q draw would correlate
    the two arms of the test and bias the balance.
    """
    if label not in CEILINGS:
        raise ValueError(f"{label!r} has no fixed place in CEILINGS, so no stable seed")
    # ⚠ `cond is None` is the SHARED p-side draw and gets slot len(CONDS), its own
    # place -- NOT slot 0, which would silently be D1_extremal's and reintroduce common
    # random numbers between the shared draw and a per-condition one.
    if cond is None:
        ci = len(CONDS)
    elif cond in CONDS:
        ci = CONDS.index(cond)
    else:
        raise ValueError(f"{cond!r} is not a condition and is not the shared draw")
    return int(np.random.SeedSequence(
        [SEED_BASE, _h(target), CEILINGS.index(label), draw, ci,
         {"p": 0, "q": 1}[hyp]]).generate_state(1)[0])


def target_state(target: str) -> dict:
    """ARM 1 state for one held-out individual: K = 86 panel, its chain, its beta.

    Built exactly as `scripts/arm1_release.run_target` builds it, including the chr21
    GENOME-WIDE share. ⚠⚠ THE SHARE IS NOT 1.0: chr21 carries ~0.01502, which is the
    scope the `eta` was calibrated at. share=1.0 would make beta per position 66.56x
    larger and, since `c = kl/eta^2` goes as share^2, the cell labelled 0.55 would sit
    at `kl` ~ 88.6 instead of 0.0200 -- Pinsker vacuous -- with no crash.
    """
    d = REPO / "data" / "cohort" / CHROM
    lc = loo.loo_chain(d, target)
    ic = lc.in_chain
    hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
    ids = np.load(d / "haplotype_ids.npy", allow_pickle=True)
    index = loo.sample_index(ids)
    panel = loo.masked_haplotypes(d, target, haplotypes=hap_all)[:, ic]
    tgt = hap_all[list(index[target])][:, ic].T                 # (T, 2)
    m = hmm.from_arrays(panel, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic],
                        cohort_id=f"{CHROM}_loo_{target}")
    share = loo_share_for(target, json.loads(SURVEY.read_text()))
    util = utility.beta_v(lc.support, in_chain=ic, share=share, chrom=CHROM,
                          weights=WEIGHTS)
    # ⚠ The DECOY POOL is the other 43 COHORT MEMBERS, on THIS target's chain axis.
    # In arm 2 the pool is the nine other prepared targets; here every member is
    # available and in the panel, which is the asymmetry the module docstring states.
    others = {s: hap_all[list(r)][:, ic].T
              for s, r in index.items() if s != target}
    return {"target": target, "panel": panel, "model": m, "util": util,
            "base": ffbs.baseline_log_norm(m), "path": tgt, "others": others,
            "t_chain": int(lc.t_chain), "k_states": int(lc.k_states), "share": share}


def _concordance(a: np.ndarray, b: np.ndarray) -> float:
    """Share of positions where two diploid paths carry the same allele pair.
    Positions either one does not call are excluded."""
    known = (a >= 0).all(axis=1) & (b >= 0).all(axis=1)
    if not known.any():
        return 0.0
    return float((a[known] == b[known]).all(axis=1).mean())


def pick_decoys(st: dict, rng: np.random.Generator) -> dict:
    """The three pre-registered conditions, D1 built by the SHARED `extremal`.

    ⚠ `extremal` is imported from `scripts/exact_kl_stable.py` rather than reimplemented,
    because the fitted `eta` was solved against THAT construction. A second copy here
    would let the decoy the ceiling was computed for drift from the decoy the attack
    scores, and nothing would raise.
    """
    p = st["path"]
    conc = {s: _concordance(p, q) for s, q in st["others"].items()}
    nearest = max(conc, key=conc.get)
    keys = sorted(st["others"])
    # ⚠⚠ D3 MUST NOT BE D2 -- the same defect and the same fix as
    # `arm2_attack.pick_decoys`, measured HERE as 2 of 44 targets: HG00741 and
    # HG01123 both drew their own nearest donor (HG01928, HG01978), so the pooled
    # D3 row shared 2/44 of its trials IDENTICALLY with D2 while being reported as
    # an independent condition. Rejection keeps every non-colliding target's D3
    # unchanged, because the rng is per target (`default_rng([DECOY_SEED, _h(t)])`).
    if len(keys) < 2:
        raise ValueError(
            f"the decoy pool is {len(keys)} donor(s), so D3_random cannot differ "
            f"from D2_nearest ({nearest}).")
    rand = keys[int(rng.integers(len(keys)))]
    while rand == nearest:
        rand = keys[int(rng.integers(len(keys)))]
    return {
        "D1_extremal": {"path": extremal(p), "id": "constructed", "concordance": None},
        "D2_nearest": {"path": st["others"][nearest], "id": nearest,
                       "concordance": conc[nearest]},
        "D3_random": {"path": st["others"][rand], "id": rand,
                      "concordance": conc[rand]},
    }


def _a3_for(st: dict) -> dict:
    """`a3_sites` on THIS target's panel, with the unkeyed-cache hazard closed BY KEY.

    ⚠⚠ `arm2_attack_a3a4._A3` is a module global with NO KEY, and Arm 1's bi-allelic
    axis is per target. Measured `t_chain` spans only 332,516..332,521, so MANY TARGETS
    SHARE A T -- which is exactly why a length check cannot be the gate: it passes on
    another target's axis of the same length, and every A3 number for that target is
    then computed on the wrong sites, silently.

    ⚠⚠ AND CLEARING IT BLINDLY IS NOT A GATE EITHER -- that was the first version of
    this function, and it made the length assertion below UNFALSIFIABLE: nothing can be
    stale if the cache is always empty, so the assertion could never fire and was
    decoration. Measured by `tests/test_arm1_attack.py`, which could not make it raise.

    So the cache is KEYED on a digest of the panel it was derived from:
      * same key   -> serve it, which is the whole point of a cache;
      * other key  -> clear and recompute, the stale-serve path, caught by KEY and not
                      by length;
      * NO key     -> RAISE. A non-empty `_A3` this driver did not stamp was populated
                      by some other caller, so whose chain it belongs to is unknowable
                      and a length check cannot tell. Refusing beats guessing.
    """
    panel = np.ascontiguousarray(st["panel"])
    key = hashlib.sha256(panel.tobytes()).hexdigest()[:16]
    c = A3M._A3
    if c:
        if "_panel_key" not in c:
            raise RuntimeError(
                "arm2_attack_a3a4._A3 holds a bi-allelic axis this driver did not "
                "stamp, so which target's chain it belongs to is unknowable. Arm 1 "
                "targets share a t_chain, so a length check cannot tell them apart. "
                "Refusing to serve it -- clear it at the call site if it is yours.")
        if c["_panel_key"] != key:
            c.clear()
        else:
            return c
    s = A3M.a3_sites({"hap": panel})
    assert s is A3M._A3, "a3_sites stopped using the module cache; re-check this guard"
    assert len(s["mask"]) == st["t_chain"], (
        f"{st['target']}: a3_sites returned a {len(s['mask'])}-site axis for a "
        f"{st['t_chain']}-position chain -- its own derivation disagrees with the chain.")
    s["_panel_key"] = key
    return s


def _u_dip(y: np.ndarray, path: np.ndarray, beta) -> float:
    return utility.diploid_utility(
        *[float(utility.u_of_output(y[:, h], path[:, h], beta)[0]) for h in (0, 1)])


def ceiling_from_kl(kl: float) -> tuple:
    """(TV upper bound, acc_ceiling) by Pinsker. ⚠ Pinsker is an INEQUALITY: the true
    ceiling is closer to 0.5 than this, so an excess over it is conservative evidence."""
    tv = min(1.0, math.sqrt(max(kl, 0.0) / 2.0))
    return tv, 0.5 * (1.0 + tv)


def run_cell(job: tuple) -> dict:
    """One (target, label): n draws under each hypothesis, A1 and A3, all three decoys.

    State is built ONCE per target by the caller and the labels looped inside, because
    `hmm.from_arrays` over 332k positions at K = 86 dominates the cost.
    """
    target, labels, etas, draws = job
    t0 = time.time()
    st = target_state(target)
    sites = _a3_for(st)
    rng = np.random.default_rng([DECOY_SEED, _h(target)])
    decoys = pick_decoys(st, rng)
    p, beta = st["path"], st["util"].beta
    out = {"target": target, "t_chain": st["t_chain"], "k_states": st["k_states"],
           "share": st["share"], "cells": {}}

    for label in labels:
        eta = etas[label]
        # ⚠ The PAIR quantities first, once per condition: `kl_stable` gives the per-pair
        # ceiling AND `d_logz` (= log Z_p - log Z_q, baseline cancelled exactly), which
        # is the only normaliser term Lambda needs.
        pair = {}
        for cond in CONDS:
            q = decoys[cond]["path"]
            k = kl_stable(st["model"], st["panel"], beta, p, q, eta=eta)
            tv, ceil = ceiling_from_kl(k["kl"])
            pair[cond] = (q, k, tv, ceil)
        # ⚠⚠ ONE RELEASE UNDER p SERVES ALL THREE CONDITIONS, which is correct and not a
        # shortcut: a draw from `Q_p` does not depend on q. It is also what makes a
        # powered run affordable -- a release over 332k positions at K = 86 costs ~40 s,
        # so a p-draw per condition would triple the dominant cost for nothing.
        # ⚠ BUT IT MAKES THE THREE CONDITIONS CORRELATED through the shared draw, so
        # they must NEVER be pooled into a single accuracy. `_report` keeps them apart.
        lam = {c: {"p": [], "q": []} for c in CONDS}
        a3 = {c: {"p": [], "q": []} for c in CONDS}
        for i in range(draws):
            rp = ffbs.release(st["model"], st["panel"], p, st["util"], eta=eta,
                              seed=cell_seed(target, label, i, None, "p"),
                              target_id=target, loo_id=target,
                              baseline_norm=st["base"], stream_key=CHROM)
            y = rp.released_path
            u_p_y = _u_dip(y, p, beta)
            a3_p_y = A3M.a3_score(y, p, sites)
            for cond in CONDS:
                q, k, _, _ = pair[cond]
                lam[cond]["p"].append(float(A.lrt_statistic(
                    np.array([u_p_y]), np.array([_u_dip(y, q, beta)]),
                    eta, 0.0, -k["d_logz"])[0]))
                a3[cond]["p"].append(a3_p_y - A3M.a3_score(y, q, sites))
                # The q-side needs its OWN release per condition: the input path differs.
                # FIXED MECHANISM -- p's model, panel, util and eta, with q as the input.
                rq = ffbs.release(st["model"], st["panel"], q, st["util"], eta=eta,
                                  seed=cell_seed(target, label, i, cond, "q"),
                                  target_id=f"{target}~{cond}", loo_id=target,
                                  baseline_norm=st["base"], stream_key=CHROM)
                yq = rq.released_path
                lam[cond]["q"].append(float(A.lrt_statistic(
                    np.array([_u_dip(yq, p, beta)]), np.array([_u_dip(yq, q, beta)]),
                    eta, 0.0, -k["d_logz"])[0]))
                a3[cond]["q"].append(A3M.a3_score(yq, p, sites)
                                     - A3M.a3_score(yq, q, sites))
        cell = {"eta": eta, "conditions": {}}
        for cond in CONDS:
            q, k, tv, ceil = pair[cond]
            cell["conditions"][cond] = {
                "decoy": decoys[cond]["id"],
                "concordance": decoys[cond]["concordance"],
                "kl": k["kl"], "d_logz": k["d_logz"], "tv_pinsker": tv,
                "acc_ceiling": ceil,
                "A1": {"lam_p": lam[cond]["p"], "lam_q": lam[cond]["q"]},
                "A3": {"lam_p": a3[cond]["p"], "lam_q": a3[cond]["q"]},
            }
        # ⚠ D1 is asserted the max-TV decoy per cell, so the directory label's
        # conservatism is MEASURED rather than assumed -- the spec's one ✅ finding.
        if eta > 0.0:
            c1 = cell["conditions"]["D1_extremal"]["acc_ceiling"]
            for cond in ("D2_nearest", "D3_random"):
                cc = cell["conditions"][cond]["acc_ceiling"]
                assert cc <= c1 + 1e-9, (
                    f"{target}/{label}: {cond} ceiling {cc:.6f} EXCEEDS D1_extremal's "
                    f"{c1:.6f}, so the directory label is NOT conservative for this "
                    f"cell and every table keyed on it understates the bound.")
        out["cells"][label] = cell
    out["seconds"] = time.time() - t0
    return out


def _acc(lp: list, lq: list) -> tuple:
    """Balanced accuracy, SE, TPR, TNR. Ties split, A1's convention."""
    a, b = np.asarray(lp, float), np.asarray(lq, float)
    tpr = float(((a > 0).sum() + 0.5 * (a == 0).sum()) / len(a))
    tnr = float(((b < 0).sum() + 0.5 * (b == 0).sum()) / len(b))
    bal = A.empirical_accuracy(a, b)
    assert abs(bal - (0.5 * tpr + 0.5 * tnr)) < 1e-12, "accuracy definitions drifted"
    se = 0.5 * math.sqrt(tpr * (1 - tpr) / len(a) + tnr * (1 - tnr) / len(b))
    return bal, se, tpr, tnr


def _report(res: list, labels: tuple, draws: int, rev: str) -> None:
    n_t = len(res)
    print(f"# ARM 1 — A1 AND A3 ON THE `acc_ceiling` GRID\n")
    print(f"chr21 ONLY. K = 86 per-target leave-one-out. {n_t} targets x "
          f"{len(labels)} labels x {draws} draws per hypothesis per condition.")
    print(f"Weighting `{WEIGHTS}`; chr21 `beta_share` per target (~0.01502, NOT 1.0). "
          f"Code revision `{rev}`.")
    print("Methodology: `docs/plan.md` direction block, 2026-10-24; spec "
          "`docs/memory.md` 2026-10-27 (B-ATK-arm1-spec).\n")
    print("⚠⚠ POSING: FIXED MECHANISM — one mechanism (p's panel, chain, beta, eta), "
          "two input paths. That is what Theorem 1 bounds. The PROCEDURE posing "
          "(`Q_q` from q's own leave-one-out panel) is TWO mechanisms, is NOT bounded, "
          "and is not run.")
    print("⚠⚠ `acc_ceiling` IS A ROW, PER (target, label, condition) — the directory "
          "label is the D1_extremal ceiling ALONE. Reading a D2/D3 accuracy against the "
          "label would test it against a bound ~0.12 too loose.")
    print("⚠ The pooled accuracy is an AVERAGE over targets, so the bound it is tested "
          "against is the MEAN of the per-pair ceilings. The max is printed beside it; "
          "testing against the max would be a looser test that could hide an excess.")
    print("⚠ Pinsker is an INEQUALITY: every ceiling here is an OVER-estimate, so the "
          "true one is closer to 0.5 and an excess is conservative.")
    print("⚠ The 0.50 column is `eta` = 0, where `Q_p = Q_q = R_D` identically: the "
          "true accuracy is EXACTLY 0.5 for every attacker. It is a HARNESS CHECK.\n")

    viol, degen, near = [], [], []
    for atk in ("A1", "A3"):
        name = ("A1 — `lrt_pair`, the Neyman-Pearson OPTIMAL distinguisher"
                if atk == "A1" else
                "A3 — PanMixer's rarity-weighted gap score, as a 1-of-2 discriminator")
        print(f"\n## {name}\n")
        for cond in CONDS:
            print(f"### {cond}\n")
            print("| metric | " + " | ".join(labels) + " |")
            print("|---|" + "---|" * len(labels))
            rows: dict = {k: [] for k in
                          ("balanced_accuracy", "se", "ci95", "tpr", "tnr",
                           "ceil_mean", "ceil_max", "excess", "z", "verdict")}
            for lab in labels:
                cells = [r["cells"][lab]["conditions"][cond] for r in res]
                lp = [v for c in cells for v in c[atk]["lam_p"]]
                lq = [v for c in cells for v in c[atk]["lam_q"]]
                bal, se, tpr, tnr = _acc(lp, lq)
                # ⚠⚠ WHICH CEILING A POOLED ACCURACY IS BOUNDED BY, and it is the MEAN.
                # With equal draws per target the pooled balanced accuracy IS the average
                # of the per-target accuracies, so if each is below its own pair ceiling
                # then the average is below the AVERAGE of those ceilings. The mean is
                # therefore the correct bound and the test runs against it. The MAX is
                # printed beside it as the conservative reading -- it is the bound a
                # reader gets if they assume the worst pair dominates -- but testing
                # against the max would be a LOOSER test that could hide a real excess,
                # which is the opposite of what a violation check is for.
                cm = stats.fmean([c["acc_ceiling"] for c in cells])
                cx = max(c["acc_ceiling"] for c in cells)
                ex = bal - cm
                # ⚠⚠ se = 0 with a POSITIVE excess is perfect separation, not a
                # tie: forcing z to 0.0 there files a total break as "within".
                # A constant answer (tpr 1 / tnr 0) scores 0.5 and is caught by
                # DEGENERATE_FLOOR below, so it never reaches this.
                z = ex / se if se > 0 else (float("inf") if ex > 0 else 0.0)
                bad = min(tpr, tnr) < DEGENERATE_FLOOR
                rows["balanced_accuracy"].append(f"{bal:.6f}")
                rows["se"].append(f"{se:.6f}")
                # ⚠ `acc_ci95`, the ERROR BAR. The SE row alone leaves the reader to
                # do the arithmetic that decides whether an excess is inside noise;
                # where a cell sits above its pair ceiling, the claim being made is
                # that this interval reaches BELOW that ceiling.
                _lo, _hi = _interval(bal, se)
                rows["ci95"].append(
                    f"[{_lo:.6f}, {_hi:.6f}]" if se > 0 else "DEGENERATE (se = 0)")
                rows["tpr"].append(f"{tpr:.4f}")
                rows["tnr"].append(f"{tnr:.4f}")
                rows["ceil_mean"].append(f"{cm:.6f}")
                rows["ceil_max"].append(f"{cx:.6f}")
                rows["excess"].append(f"{ex:+.6f}")
                rows["z"].append(f"{z:+.2f}" if se > 0 else "n/a")
                v = "DEGENERATE" if bad else ("⚠ EXCEEDS" if z > 1.96 else "within")
                rows["verdict"].append(v)
                if bad:
                    degen.append((atk, cond, lab, tpr, tnr))
                elif z > 1.96:
                    viol.append((atk, cond, lab, bal, cm, z))
                elif ex > 0:
                    # ⚠⚠ ABOVE THE CEILING BUT NOT SIGNIFICANT. These reached the
                    # table's `excess` and `z` rows and stopped there: the SYNTHESIS
                    # listed only z > 1.96, so a reader skimming the summary saw "0
                    # cells exceed" while two did by point estimate. Shown now, with
                    # the interval that qualifies them, rather than absorbed.
                    near.append((atk, cond, lab, bal, cm, z, _lo))
            lbl = {"balanced_accuracy": "balanced accuracy*", "se": "SE",
                   "ci95": "`acc_ci95` [THE ERROR BAR]",
                   "tpr": "TPR (P_p[Lambda>0])", "tnr": "TNR (P_q[Lambda<0])",
                   "ceil_mean": "`acc_ceiling` THIS pair, mean — THE BOUND*",
                   "ceil_max": "`acc_ceiling` THIS pair, max (conservative)",
                   "excess": "excess over the mean", "z": "z", "verdict": "verdict"}
            for k, v in rows.items():
                print(f"| {lbl[k]} | " + " | ".join(v) + " |")
            d = [r["cells"][labels[-1]]["conditions"][cond] for r in res]
            print(f"\n- decoy at {labels[-1]}: "
                  f"`{d[0]['decoy']}`"
                  + (f", concordance {stats.median([x['concordance'] for x in d]):.4f} "
                     f"median" if d[0]["concordance"] is not None else "")
                  + f"; `d_logz` median {stats.median([x['d_logz'] for x in d]):+.3f}")
            print()

    print("\n## SYNTHESIS\n")
    print(f"- **{len(viol)} cell(s) exceed their own pair `acc_ceiling` at z > 1.96** "
          f"and **{len(degen)} are DEGENERATE** (min(TPR, TNR) < "
          f"{DEGENERATE_FLOOR:g}).")
    print(f"- **{len(near)} cell(s) exceed their pair `acc_ceiling` by POINT ESTIMATE "
          f"without reaching z > 1.96.** Reported, not claimed:")
    for atk, cond, lab, bal, ceil, z, lo in near:
        reaches = "reaches below" if lo <= ceil else "⚠⚠ does NOT reach below"
        print(f"  - {atk} {cond} @ {lab}: {bal:.6f} against a pair ceiling of "
              f"{ceil:.6f} (+{bal - ceil:.6f}, z = {z:+.2f}); its 95% lower limit "
              f"{lo:.6f} {reaches} the ceiling.")
    for atk, cond, lab, tpr, tnr in degen:
        print(f"  - DEGENERATE {atk} {cond} @ {lab}: TPR {tpr:.4f} / TNR {tnr:.4f} — a "
              f"CONSTANT ANSWER, not a null. It must not be read as the mechanism "
              f"holding.")
    for atk, cond, lab, bal, ceil, z in viol:
        print(f"  - ⚠⚠ {atk} {cond} @ {lab}: {bal:.6f} against a pair ceiling of "
              f"{ceil:.6f}, z = {z:+.2f}. Pinsker runs CONSERVATIVE, so a SIGNIFICANT "
              f"excess falsifies the IMPLEMENTATION, not the theorem.")
    print("- ⚠ An attacker that fails is not evidence the bound holds (standing "
          "constraint 4). A1 is the Neyman-Pearson optimum for its pair, so ITS null is "
          "the strongest available null — but A3's is not, and a DEGENERATE A3 row is no "
          "evidence at all.")
    print("- ⚠ D2/D3 decoys are cohort members IN the K = 86 panel while the target is "
          "HELD OUT, so the decoy is exactly reproducible by the mosaic and the target "
          "is not. That asymmetry is legitimate and inside the bounds, but it drives "
          "TPR and TNR apart; read those rows, not the balanced figure alone.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=5)
    ap.add_argument("--targets", type=int, default=0, help="first N targets (0 = all)")
    ap.add_argument("--labels", nargs="+", default=list(CEILINGS))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)
    for lab in a.labels:
        if lab not in CEILINGS:
            raise SystemExit(f"unknown label {lab!r}; CEILINGS is {CEILINGS}")
    labels = tuple(a.labels)

    fit = json.loads(FIT.read_text())
    assert fit["arm"] == "arm1" and fit["scope"] == CHROM, fit["arm"]
    assert fit["weights"] == WEIGHTS, (
        f"the fit was solved at weights={fit['weights']!r} but this script runs at "
        f"{WEIGHTS!r}; one eta would sit at two different acc_ceiling (constraint 5).")
    per = fit["per_target"]
    targets = sorted(per)
    if a.targets:
        targets = targets[: a.targets]

    # ⚠ THE FIT HAS NO '0.50' KEY (six labels, 0.55..0.80). `acc_ceiling` 0.50 is
    # `eta` = 0, where Q_p = Q_q = R_D identically -- so 0.0 is the RIGHT eta, not a
    # missing value to be filled in, and the column is a harness check.
    etas = {t: {lab: (0.0 if lab == "0.50" else per[t]["refined"][lab]["eta"])
                for lab in labels} for t in targets}

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    rev = os.popen("git -C %s rev-parse --short HEAD" % REPO).read().strip() or "unknown"
    print(f"# {len(targets)} targets x {len(labels)} labels x {a.draws} draws, "
          f"{a.workers} workers", file=sys.stderr)
    res, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_cell, (t, labels, etas[t], a.draws)): t for t in targets}
        for f in as_completed(futs):
            r = f.result()
            res.append(r)
            print(f"  {r['target']}  {r['seconds']:.1f}s  "
                  f"({len(res)}/{len(targets)})", file=sys.stderr)
    # ⚠ assert completeness from the RESULT SET, not from the loop: a filter cannot
    # tell "found nothing" from "ran nothing".
    got = {r["target"] for r in res}
    if got != set(targets):
        raise SystemExit(f"incomplete: missing {sorted(set(targets) - got)}")
    res.sort(key=lambda r: r["target"])
    print(f"# {len(res)}/{len(targets)} targets in {time.time() - t0:.1f}s",
          file=sys.stderr)
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm1", "chrom": CHROM, "weights": WEIGHTS, "draws": a.draws,
             "labels": list(labels), "revision": rev, "results": res}, indent=1) + "\n")
        print(f"# JSON -> {a.json}", file=sys.stderr)
    _report(res, labels, a.draws, rev)
    return 0


if __name__ == "__main__":
    sys.exit(main())
