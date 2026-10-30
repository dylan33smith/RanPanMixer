#!/usr/bin/env python
"""Arm 2's A3 and A4 — the two pre-registered attacks that had never been run.

⚠ IT REFUSES TO RUN WITHOUT A SIGNED PRE-REGISTRATION. Standing constraint 6 froze
the attacker set on 2026-10-10; A3 and A4 are rows 3 and 4 of that table and nothing
here chooses them after the fact.

    | A3 | rarity-weighted genotype linkage, re-posed as 1-of-2 | both |
    |    | PanMixer's score, with `ref_af = 1 - alt_af` fixed | 0.5 |
    |    | must link the UNOBFUSCATED target at near-100%, and must sit at chance
    |    | against a `tau` = 0 release |
    | A4 | segment-restricted likelihood ratio over one haplotype block | both |
    |    | A1 restricted to a block, multiplicity-corrected | 0.5 |
    |    | must find a planted block copied verbatim from the target |

WHY THESE TWO AT ALL, GIVEN A1 IS PROVABLY OPTIMAL. A1's statistic is the
Neyman-Pearson ratio for the WHOLE released path, so no attacker beats it on that
game. A3 and A4 change the game rather than trying to beat it:
  * **A3 is the attack a real adversary runs.** It is PanMixer's own linkage score,
    which needs no knowledge of `eta_tau`, of `beta_v`, of `Z_p` or that a mechanism
    was applied at all. If a weaker, realistic attacker beat the optimal one, the
    implementation would be wrong; the reason to run it is that it is the attack the
    threat model actually faces, and its weights come from PUBLIC allele frequencies.
  * **A4 asks whether the signal CONCENTRATES.** A1 sums one statistic over the whole
    chain, so a tilt that is loud in one block and silent elsewhere is averaged away.
    A4 restricts the same statistic to a block and maximises over blocks, which is
    strictly stronger per-block and strictly weaker after the multiplicity correction
    that maximising requires. Budget dilution is this project's central result, and
    "is the dilution uniform?" is the question A4 answers.

⚠⚠ THREE THINGS THE PRE-REGISTRATION LEFT OPEN, FIXED HERE BEFORE ANY OUTCOME WAS
INSPECTED, because choosing them afterwards is the garden of forking paths:
  1. **A3's allele frequencies come from the PUBLIC COHORT**, not from the target and
     not from an external panel. Constraint 1: the adversary may use `G` and `D`.
     Cohort-derived `alt_af` is exactly what a real adversary computes.
  2. **A4's block grid is {512, 4096, 32768} positions, all three reported**, never
     a best-of. 64 was considered and dropped with a reason that is not an outcome:
     a 64-position block carries ~2e-4 of the chain's weight, so its statistic is
     dominated by the baseline's own randomness and the multiplicity correction over
     ~5,200 blocks would swamp it. Smaller is not better here; it is noisier.
  3. **A4's planted block, for the positive control, is the block with the MOST
     p-vs-q differing called positions.** Deterministic, computable before any
     release is drawn, and the best case for the attacker -- which is what a positive
     control is for. Planting in a block where `p` and `q` agree would test nothing.

⚠ THE BOUND APPLIES TO THE BALANCED ACCURACY ONLY, so both attacks draw from BOTH
hypotheses. `(1 + tau)/2` bounds `0.5*P_p(L>0) + 0.5*P_q(L<0)` and nothing else; a
one-sided rate is not bounded by it and an attacker that always answers "p" scores
1.0 one-sided while sitting exactly at chance balanced.

    python scripts/arm2_attack_a3a4.py [--draws 250] [--targets 10] [--json out.json]
    python scripts/arm2_attack_a3a4.py --controls-only     # cheap; no releases
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import attacks as A
import store  # noqa: E402
import runner as R  # noqa: E402
from mechanism import calibration as C  # noqa: E402
from mechanism import ffbs, utility  # noqa: E402

# ⚠ The shared machinery is IMPORTED from A1's runner, not copied. `load_state`,
# `pick_decoys` and `concordance` define the cohort, the model, the decoy conditions
# and the seeding discipline; a second copy would drift and then the two attacks
# would no longer be running on the same experiment.
import arm2_attack as A1  # noqa: E402  -- for the WEIGHTS global, which must be SET
from arm2_attack import (CEIL_SEED_OFFSET, CHROM, DECOY_SEED, ETA_GRID,  # noqa: E402
                         ETA_SEED_OFFSET, REPORTED_CEILINGS, ceiling_grid,
                         EXPLORATORY_TAUS, SEED_GRID, TAUS, load_state,
                         pick_decoys)

SEED_BASE = 70_000           # disjoint from A1's 50_000
BLOCK_SIZES = (512, 4096, 32768)

# ⚠⚠⚠ A3's PRE-REGISTERED POSITIVE CONTROL CANNOT FAIL, AND A GATE THAT CANNOT FAIL IS
# DECORATION (CLAUDE.md, Agent Behavior). Proved algebraically and checked on the real
# data 2026-10-13: their score sums `-log P_HWE` ONLY over sites where the genotypes
# AGREE, mismatches contribute zero, and the weights are non-negative (min 0.022857).
# At a site where `g_p == g_q` the two candidates' weights are identical. So
#     score(p, q) = sum over {v : g_p = g_q}  w(g_p)   <=   sum over {v : p called} w(g_p) = score(p, p)
# -- a subset sum over non-negative terms -- and `Lambda > 0` is FORCED the moment one
# site has `g_p != g_q` (measured: 45,411 such sites for the first target against its
# nearest decoy, score 96,316.8 against 43,294.9). The control therefore verifies the
# PLUMBING (site axis, genotype encoding, finite non-negative weights) and says nothing
# whatever about POWER. It is reported as passing AND as unable to fail.
#
# So a real power measurement is added, LABELLED A DIAGNOSTIC rather than an attack
# outcome -- which is what constraint 6 permits without a new pre-registration, exactly
# as PanMixer's evaluators were run as a labelled descriptive column. It asks the
# question the broken control was meant to answer: HOW MUCH verbatim target copy must a
# release contain before A3 notices? A tau = 0 release (zero target information) has a
# declared fraction `f` of its positions overwritten with the target verbatim, and the
# balanced accuracy is measured against `f`. The subsets are PREFIXES OF ONE RANDOM
# PERMUTATION, so the curve is monotone by construction and no subset-choice noise is
# mistaken for a threshold.
LEAK_F = (0.0, 0.0001, 0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 1.0)
POWER_SEEDS = 5
POWER_SEED_BASE = 90_000

_A3: dict = {}


# ---------------------------------------------------------------------------
# A3 — PanMixer's rarity-weighted genotype linkage, re-posed as 1-of-2
# ---------------------------------------------------------------------------

def _cell_seed(ti: int, pointi: int, k: int) -> int:
    """This driver's per-cell release stream, keyed on the WHOLE cell.

    ⚠ The `34` tag keeps these streams disjoint from `arm2_attack.py`'s, which uses the
    same SEED_BASE and the same (target, point, draw) key; without it the two attacks
    would silently score the SAME releases.
    ⚠ Hoisted to module level so the seed-uniqueness gate tests the REAL function.
    """
    return int(np.random.SeedSequence(
        [SEED_BASE, 34, ti, pointi, k]).generate_state(1)[0])


def a3_sites(st: dict) -> dict:
    """The bi-allelic site axis, the cohort `alt_af`, and the HWE rarity weights.

    THEIR SCORE, STATED EXACTLY (`gap_score.py:30-83`). For a candidate `c` and a
    released path `y`, summing over sites where BOTH are fully called:

        score(y, c) = sum_v  1[g_v(y) == g_v(c)] * (-log(P_HWE(g_v(c)) + 1e-10))

    with `g` the count of ALT alleles in {0,1,2} and
    `P_HWE(0) = q^2`, `P_HWE(1) = 2pq`, `P_HWE(2) = p^2`, `p = alt_af`, `q = 1 - p`.

    ⚠ FOUR PROPERTIES OF THEIR SCORE THAT DECIDE WHAT A NULL MEANS:
      * a MISMATCH contributes ZERO, not a penalty. The score rewards agreement and
        never punishes disagreement, so it rises with the number of sites a candidate
        happens to share with the release for any reason.
      * the weight is a property of the CANDIDATE, not of the release, so
        `score(y, p)` and `score(y, q)` are weighted differently. That is their
        design and is preserved.
      * RARE matching genotypes dominate: `-log` of a small HWE probability is large.
      * `ref_af = 1 - alt_af` is FIXED, per the pre-registration. Their own
        `build_cache` keeps only sites whose stored ref and alt frequencies already
        sum to 1 (`gap_score.py:101-102`), so this matches their effective axis while
        removing the dependence on a stored column that need not be consistent.

    ⚠ BI-ALLELIC IS DERIVED FROM THE COHORT, which is public. A position qualifies
    when its called donor alleles take EXACTLY TWO distinct codes; the lower code is
    REF and the higher is ALT. Our allele codes run to 88 (multi-allelic sites exist),
    and their genotype-as-a-sum encoding is only defined on bi-allelic sites.
    """
    if _A3:
        return _A3
    hap = st["hap"]                                  # (K, T) allele codes, -1 = no call
    called = hap >= 0
    n_called = called.sum(axis=0)
    lo = np.where(called, hap, np.int16(1 << 14)).min(axis=0)
    hi = np.where(called, hap, np.int16(-1)).max(axis=0)
    # exactly two distinct called codes: two values present, and nothing outside them
    only_two = ((~called) | (hap == lo[None, :]) | (hap == hi[None, :])).all(axis=0)
    bi = (hi > lo) & only_two & (n_called > 0)

    ref, alt = lo[bi], hi[bi]
    alt_count = ((hap[:, bi] == alt[None, :]) & called[:, bi]).sum(axis=0)
    p = alt_count / n_called[bi]
    # Their line 42 removes p == 0 or p == 1. On a cohort-derived axis both codes are
    # present by construction, so this cannot bind -- ASSERT it rather than assume it,
    # because if it ever binds the weights go to -log(1e-10) and the score explodes.
    assert ((p > 0) & (p < 1)).all(), "a bi-allelic site has a degenerate alt_af"
    q = 1.0 - p
    # -log of the HWE probability of each genotype, per site. Their +1e-10 preserved.
    w = np.stack([-np.log(q * q + 1e-10),
                  -np.log(2.0 * p * q + 1e-10),
                  -np.log(p * p + 1e-10)])        # (3, n_bi)
    _A3.update(mask=bi, ref=ref, alt=alt, alt_af=p, logw=w, n_bi=int(bi.sum()))
    return _A3


def a3_genotype(path: np.ndarray, s: dict) -> tuple:
    """(g, valid) on the bi-allelic axis. `g` is the ALT count; valid = both called."""
    sub = path[s["mask"]]                                   # (n_bi, 2)
    ok = (sub >= 0).all(axis=1)
    g = ((sub[:, 0] == s["alt"]).astype(np.int8)
         + (sub[:, 1] == s["alt"]).astype(np.int8))
    return g, ok


def a3_score(y: np.ndarray, cand: np.ndarray, s: dict) -> float:
    """Their score for one (release, candidate) pair."""
    gy, vy = a3_genotype(y, s)
    gc, vc = a3_genotype(cand, s)
    both = vy & vc
    w = s["logw"][gc, np.arange(len(gc))]
    return float((w * ((gy == gc) & both)).sum())


# ---------------------------------------------------------------------------
# A4 — the same statistic as A1, restricted to a block and maximised over blocks
# ---------------------------------------------------------------------------

def _edges(t: int, size: int) -> np.ndarray:
    return np.arange(0, t, size)


def a4_block_stats(y: np.ndarray, p: np.ndarray, q: np.ndarray,
                   beta: np.ndarray, eta: float) -> dict:
    """`Lambda` restricted to each block, for every declared block size.

    Per position, summed over both haplotypes, the A1 statistic's data-dependent part:

        contrib_v = (eta/2) * beta_v * ( 1[y_v = p_v] - 1[y_v = q_v] )

    ⚠ `eta/2` IS THE REALISED EXPONENT, not `eta`. The diploid release factorises and
    each haplotype runs at half the budget; reporting the full `eta` here would be the
    factor-of-two error this project has already made twice.

    ⚠ THE `log Z` TERM IS ABSENT, AND THAT IS CORRECT FOR A BLOCK. A1 subtracts
    `log Z_p - log Z_q`, which does NOT decompose over blocks. But that term is a
    constant for the pair -- it does not depend on `y` -- so within one block it is a
    fixed offset, and a sign rule at threshold 0 is the restricted likelihood ratio up
    to that offset. The offset is absorbed where a threshold is chosen; it cannot
    create block-to-block structure, which is what A4 is looking for.

    ⚠ `1[y = p]` uses the SAME match rule as `utility.u_of_output`: a -1 never
    matches, so a position neither candidate calls contributes exactly 0 to both.
    """
    mp = np.zeros(len(beta))
    mq = np.zeros(len(beta))
    for h in (0, 1):
        yh, ph, qh = y[:, h], p[:, h], q[:, h]
        mp += ((yh == ph) & (ph >= 0) & (yh >= 0))
        mq += ((yh == qh) & (qh >= 0) & (yh >= 0))
    contrib = (eta / 2.0) * beta * (mp - mq)
    return {s: np.add.reduceat(contrib, _edges(len(beta), s)).astype(np.float32)
            for s in BLOCK_SIZES}


# ---------------------------------------------------------------------------
# the trial
# ---------------------------------------------------------------------------

def _one_trial(job: tuple) -> dict:
    """One release under `Q_p` and one under `Q_q` per condition, scored by A3 and A4.

    ⚠ The two hypotheses get INDEPENDENT seeds, as in A1: reusing the `Q_p` seed would
    correlate the two arms of the test and bias the balance.
    """
    # ⚠ `point` is the COLUMN LABEL, `eta` is the tilt applied. In tau mode they are
    # tau and arctanh(tau); in eta mode point IS the eta and there is no
    # representable tau, so eta is PASSED rather than re-derived.
    target_id, point, eta, seed = job
    tau = point
    st = load_state()
    s3 = a3_sites(st)
    beta = st["util"].beta
    p = st["paths"][target_id]
    stable = int.from_bytes(hashlib.sha256(target_id.encode()).digest()[:4], "big")
    decoys = pick_decoys(target_id, st, np.random.default_rng([DECOY_SEED, stable]))

    rel = ffbs.release(st["model"], st["hap"], p, st["util"], eta=eta, seed=seed,
                       target_id=target_id, graph_id=CHROM,
                       baseline_norm=st["base"])
    y = rel.released_path

    out = {"target": target_id, "tau": tau, "seed": seed, "conditions": {}}
    for name, d in sorted(decoys.items()):
        qd = d["path"]
        cell: dict = {"decoy": d["id"], "concordance": d["concordance"]}

        # ---- A3 under Q_p
        cell["a3_under_p"] = a3_score(y, p, s3) - a3_score(y, qd, s3)
        # ---- A4 under Q_p
        bp = a4_block_stats(y, p, qd, beta, eta)

        # ---- the second hypothesis: a release from Q_q, independent seed
        q_seed = int(seed) ^ (0x5EED * (1 + sorted(decoys).index(name)))
        rel_q = ffbs.release(st["model"], st["hap"], qd, st["util"], eta=eta,
                             seed=q_seed, target_id=f"{target_id}~{name}",
                             graph_id=CHROM, baseline_norm=st["base"])
        yq = rel_q.released_path
        cell["a3_under_q"] = a3_score(yq, p, s3) - a3_score(yq, qd, s3)
        bq = a4_block_stats(yq, p, qd, beta, eta)

        cell["a4_under_p"] = {str(k): v.tolist() for k, v in bp.items()}
        cell["a4_under_q"] = {str(k): v.tolist() for k, v in bq.items()}
        out["conditions"][name] = cell
    return out


def _power_trial(job: tuple) -> dict:
    """A3's power against a KNOWN leak of declared size. A labelled diagnostic.

    Under each hypothesis the release is an untilted (`tau` = 0) draw -- carrying zero
    information about either candidate -- with a prefix of one random permutation of
    the chain overwritten by that hypothesis's own path, verbatim. `Lambda` is A3's
    statistic on the result, so the balanced accuracy at each `f` is A3's power against
    a leak of exactly that size.
    """
    target_id, seed = job
    st = load_state()
    s3 = a3_sites(st)
    p = st["paths"][target_id]
    t = len(st["util"].beta)
    stable = int.from_bytes(hashlib.sha256(target_id.encode()).digest()[:4], "big")
    decoys = pick_decoys(target_id, st, np.random.default_rng([DECOY_SEED, stable]))
    perm = np.random.default_rng([POWER_SEED_BASE, stable, seed]).permutation(t)

    rel_p = ffbs.release(st["model"], st["hap"], p, st["util"], 0.0, seed=seed,
                         target_id=f"{target_id}~pw", graph_id=CHROM,
                         baseline_norm=st["base"])
    out = {"target": target_id, "seed": seed, "conditions": {}}
    for name, d in sorted(decoys.items()):
        qd = d["path"]
        rel_q = ffbs.release(st["model"], st["hap"], qd, st["util"], 0.0,
                             seed=seed ^ (0x9000 + sorted(decoys).index(name)),
                             target_id=f"{target_id}~pw~{name}", graph_id=CHROM,
                             baseline_norm=st["base"])
        cell = {"decoy": d["id"], "lambda_p": [], "lambda_q": []}
        for fr in LEAK_F:
            k = int(round(fr * t))
            idx = perm[:k]
            yp = rel_p.released_path.copy(); yp[idx] = p[idx]
            yq = rel_q.released_path.copy(); yq[idx] = qd[idx]
            cell["lambda_p"].append(a3_score(yp, p, s3) - a3_score(yp, qd, s3))
            cell["lambda_q"].append(a3_score(yq, p, s3) - a3_score(yq, qd, s3))
        out["conditions"][name] = cell
    return out


# ---------------------------------------------------------------------------
# the positive controls — run FIRST, because a failed control voids the null
# ---------------------------------------------------------------------------

def controls(st: dict, tids: list) -> dict:
    """A3's and A4's pre-registered positive controls. No releases for A3.

    ⚠ THESE DECIDE WHETHER THE MAIN RESULT MEANS ANYTHING. Standing constraint 7: an
    attacker that fails is a lower bound on distinguishability. A null from an
    attacker that cannot succeed even when handed the answer is not a lower bound on
    anything -- it is a broken instrument, and this project has already shipped four
    false findings from instruments nobody checked.
    """
    s3 = a3_sites(st)
    beta = st["util"].beta
    out: dict = {"a3": {"trials": [], "n_correct": 0, "n": 0},
                 "a4": {"trials": [], "n_correct": 0, "n": 0}}

    for tid in tids:
        p = st["paths"][tid]
        stable = int.from_bytes(hashlib.sha256(tid.encode()).digest()[:4], "big")
        decoys = pick_decoys(tid, st, np.random.default_rng([DECOY_SEED, stable]))
        for name, d in sorted(decoys.items()):
            qd = d["path"]

            # ---- A3 (a): the UNOBFUSCATED target must be linked. y IS p.
            lam = a3_score(p, p, s3) - a3_score(p, qd, s3)
            out["a3"]["trials"].append({"target": tid, "condition": name,
                                        "lambda": lam, "correct": bool(lam > 0)})
            out["a3"]["n"] += 1
            out["a3"]["n_correct"] += int(lam > 0)

            # ---- A4: a block copied VERBATIM from the target must be found.
            # The planted block is the one with the most p-vs-q differing called
            # positions -- declared in the docstring, fixed before any outcome.
            diff = np.zeros(len(beta))
            for h in (0, 1):
                diff += ((p[:, h] != qd[:, h]) & (p[:, h] >= 0) & (qd[:, h] >= 0))
            size = 4096
            ed = _edges(len(beta), size)
            per_block = np.add.reduceat(diff, ed)
            planted = int(per_block.argmax())
            assert per_block[planted] > 0, (
                f"{tid}/{name}: no block has a differing called position, so a "
                f"planted block cannot be detected and this control is void")
            lo = ed[planted]
            hi = min(lo + size, len(beta))
            # a tau = 0 release -- NO tilt at all -- then the block overwritten with
            # the target verbatim. If A4 cannot find a verbatim copy against an
            # untilted background, it cannot find anything.
            rel = ffbs.release(st["model"], st["hap"], p, st["util"], 0.0,
                               seed=SEED_BASE + planted, target_id=f"{tid}~ctl",
                               graph_id=CHROM, baseline_norm=st["base"])
            y = rel.released_path.copy()
            y[lo:hi] = p[lo:hi]
            # scored at a positive eta so the statistic is not identically zero
            bs = a4_block_stats(y, p, qd, beta, C.eta_tau(0.99))
            found = int(np.asarray(bs[size]).argmax())
            ok = found == planted
            out["a4"]["trials"].append({"target": tid, "condition": name,
                                        "planted_block": planted, "found_block": found,
                                        "differing_positions": float(per_block[planted]),
                                        "correct": bool(ok)})
            out["a4"]["n"] += 1
            out["a4"]["n_correct"] += int(ok)

    for k in ("a3", "a4"):
        out[k]["rate"] = out[k]["n_correct"] / out[k]["n"]
    return out


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

def _z_for_one_sided(alpha: float) -> float:
    """`z` with `P(Z > z) = alpha`. Reused from A1's runner by import."""
    from arm2_attack import _z_for_one_sided as f
    return f(alpha)


# ⚠ Z95, _interval, _fmt and classify_excess are IMPORTED from `arm2_attack`, never
# redefined here. Two copies of the acc-line arithmetic is how the ceiling branch
# reached one driver and not the other; one definition means one mutation point and
# one set of gates for both.
from arm2_attack import Z95, _fmt, _interval, classify_excess  # noqa: E402,F401


def _acc(lp: np.ndarray, lq: np.ndarray) -> tuple:
    """Balanced accuracy and its SE, by A1's convention. Ties split."""
    tpr = ((lp > 0).sum() + 0.5 * (lp == 0).sum()) / len(lp)
    tnr = ((lq < 0).sum() + 0.5 * (lq == 0).sum()) / len(lq)
    bal = A.empirical_accuracy(lp, lq)
    assert abs(bal - (0.5 * tpr + 0.5 * tnr)) < 1e-12, "accuracy definitions drifted"
    se = 0.5 * math.sqrt(tpr * (1 - tpr) / len(lp) + tnr * (1 - tnr) / len(lq))
    return bal, se, tpr, tnr


CONDS = ("D1_extremal", "D2_nearest", "D3_random")


def _report(res: list, ctl: dict, taus: tuple, prereg: str, rev: str,
            tids: list, n_per: int, pw: list | None = None,
            ceilings: dict | None = None,
            prior_table: dict | None = None) -> tuple[dict, int]:
    table: dict = {"a3": {}, "a4": {}}
    is_ceil = ceilings is not None
    knob = "acc_ceiling" if is_ceil else "tau"

    def bound_at(t: float) -> float:
        """⚠⚠ THE LINE A MEASURED ACCURACY IS TESTED AGAINST -- the whole reason this
        function exists rather than a bare `C.p_succ_bound(t)` call.

        In CEILING mode the column label IS the measured `acc_ceiling`, so the bound is
        the label itself. Handing that label to `p_succ_bound` instead returns
        (1 + acc)/2 -- 0.7500 at the 0.50 column, 0.9000 at 0.80 -- a line 0.10 to 0.25
        ABOVE the real one, under which a cell can exceed the true ceiling by 0.2 and
        still print as comfortably within. THAT IS WHAT THIS DRIVER DID until
        2026-10-29, for every A3 and A4 cell in `logs/B_acc_arm2/a34_ceil.md`;
        `scripts/arm2_attack.py` never did, because it got the ceiling branch and this
        file did not. Gated by tests/test_acc_line_interval.py, which injects a value
        BETWEEN the two lines -- the 0.96 injection on record could not discriminate
        them, since it clears both.
        """
        return float(t) if is_ceil else C.p_succ_bound(t)

    def bound_row(label_gate: bool = True) -> None:
        """`p_succ_bound` is VACUOUS in ceiling mode -- the realised `tau` is 1.0 to
        float64 at these `eta_tau` -- so it is printed as such and the real line gets
        its own row, matching `arm2_attack.py`'s layout exactly."""
        # ⚠⚠ EVERY NUMBER IN THIS ROW COMES FROM `bound_at`, deliberately. If the row
        # re-derived the line itself, the PRINTED bound and the CHECKED bound could
        # drift apart -- a table showing 0.80 while testing against 0.90 -- and that
        # split is strictly worse than either error alone. One source, so one mutation
        # breaks both the row and the check.
        star = " *" if label_gate else ""
        if is_ceil:
            # ⚠ NOT vacuous at every column. `ceiling_grid` sets the 0.50 column's eta
            #   to 0 BY CONSTRUCTION (the untilted control), where tau = 0 exactly and
            #   `p_succ_bound` is exactly 0.5000 -- a real bound, and the tightest one
            #   in the table. Printing "~1.0 (vacuous)" there tells a reader no
            #   proof-derived bound applies to the control column, which is backwards.
            def _psb(t) -> str:
                etas = [ceilings[_t][f"{float(t):.2f}"] for _t in ceilings]
                return "0.5000 (exact, eta = 0)" if max(etas) == 0 else "~1.0 (vacuous)"
            print(f"| `p_succ_bound`{star} | "
                  + " | ".join(_psb(t) for t in taus) + " |")
            print(f"| `acc_ceiling` [THE BOUND, measured per target]{star} | "
                  + " | ".join(f"{bound_at(t):.4f}" for t in taus) + " |")
        else:
            print(f"| `p_succ_bound`{star} | " + " | ".join(
                ("~1.0 (vacuous)" if t >= 1.0 else f"{bound_at(t):.4f}")
                for t in taus) + " |")
    print(f"methodology: {prereg} (git_rev {rev})\n")
    print(f"ARM 2. {len(tids)} targets x {len(taus)} {knob} x {n_per} draws = "
          f"{n_per * len(tids)} trials per ({knob}, condition) PER HYPOTHESIS -- "
          f"matching A1's realised design of 250. chr21, scope chain, beta_share = 1.\n")
    if is_ceil:
        print(f"⚠⚠ CEILING MODE: each of the {len(tids)} targets runs at ITS OWN "
              f"`eta_tau`, solved from its own `c_quadratic`. The column is "
              f"`acc_ceiling` and IS the bound -- measured at the solved eta and "
              f"refined once, not fitted. ⚠ Pinsker is an inequality, so the TRUE "
              f"ceiling is closer to 0.5 than every column label, and the looseness "
              f"runs in the direction that flatters us.")
        for _lab in taus:
            _v = sorted(ceilings[_t][f"{float(_lab):.2f}"] for _t in tids)
            print(f"      {float(_lab):.2f}: eta {_v[0]:,.2f} .. {_v[-1]:,.2f} "
                  f"(median {_v[len(_v) // 2]:,.2f})")
        print()

    # ---------------- the controls first
    print("### The pre-registered positive controls\n")
    print("| control | requirement | result | passed? |")
    print("|---|---|---|---|")
    a3r, a4r = ctl["a3"]["rate"], ctl["a4"]["rate"]
    a3_ok, a4_ok = a3r >= 0.95, a4r >= 0.95
    print(f"| A3 (a) links the UNOBFUSCATED target | near-100% over "
          f"{ctl['a3']['n']} (target, decoy) pairs | **{a3r:.4f}** "
          f"({ctl['a3']['n_correct']}/{ctl['a3']['n']}) | "
          f"{'**yes, but see below**' if a3_ok else '**NO**'} |")
    print(f"| A4 finds a block copied VERBATIM from the target | the argmax block is "
          f"the planted one, over {ctl['a4']['n']} pairs | **{a4r:.4f}** "
          f"({ctl['a4']['n_correct']}/{ctl['a4']['n']}) | "
          f"{'**yes**' if a4_ok else '**NO**'} |")
    print()
    print("⚠⚠⚠ **A3's PRE-REGISTERED CONTROL CANNOT FAIL, SO ITS 1.0000 IS NOT A POWER "
          "STATEMENT.** Their score sums `-log P_HWE` only where genotypes AGREE, "
          "mismatches contribute zero, and the weights are non-negative, so "
          "`score(p, q)` is a SUBSET SUM of `score(p, p)` and `Lambda > 0` is "
          "algebraically forced. The control verifies the site axis, the genotype "
          "encoding and the weights -- the plumbing -- and nothing about power. "
          "A gate that cannot fail is decoration, so the power curve below replaces "
          "it as the actual power statement.")
    print("✅ **A4's control is NOT vacuous and it passed:** a 4,096-position block "
          "copied verbatim from the target was the argmax block in "
          f"{ctl['a4']['n_correct']} of {ctl['a4']['n']} pairs, against an untilted "
          "background. A verbatim leak of that size IS detectable by A4.")
    print()

    # ---------------- A3's real power statement
    if pw:
        print("### A3's power against a KNOWN leak — a labelled DIAGNOSTIC, not an "
              "attack outcome\n")
        print("An untilted (`tau` = 0) release with a declared fraction `f` of its "
              "positions overwritten by the target VERBATIM. Balanced accuracy against "
              "`f`, so each column is A3's power against a leak of exactly that size. "
              "Subsets are prefixes of one random permutation, so the curve is monotone "
              "by construction.\n")
        fh = "| condition | " + " | ".join(f"f={x:g}" for x in LEAK_F) + " |"
        print(fh); print("|" + "---|" * (len(LEAK_F) + 1))
        for cond in CONDS:
            cells = []
            for i, _fr in enumerate(LEAK_F):
                lp = np.array([r["conditions"][cond]["lambda_p"][i] for r in pw])
                lq = np.array([r["conditions"][cond]["lambda_q"][i] for r in pw])
                bal, se, _tpr, _tnr = _acc(lp, lq)
                # ⚠ the SE was DISCARDED here until 2026-10-29 (`_acc(lp, lq)[0]`).
                # n is only 50 per arm on this curve against 250 in the main tables,
                # so its SE is ~0.048 -- more than double -- and the f = 0 column is
                # a genuine `eta_tau` = 0 release whose acc line is exactly 0.50.
                # An f = 0 cell reading 0.59 is therefore WITHIN noise of chance, and
                # without the interval it reads as signal.
                table.setdefault("power", {}).setdefault(cond, {})[str(_fr)] = {
                    "balanced": bal, "se": se, "n_per_arm": len(lp)}
                cells.append(_fmt(bal, se))
            print(f"| {cond} | " + " | ".join(cells) + " |")
        # the smallest declared leak A3 resolves at all
        for cond in CONDS:
            curve = [(x, table["power"][cond][str(x)]["balanced"]) for x in LEAK_F]
            hit = next((x for x, b in curve if b >= 0.75), None)
            z0 = table["power"][cond]["0.0"]
            # ⚠ through `_fmt`, never a hand-rolled interval: this bullet printed
            #   "[0.5000, 0.5000] ... whose 95% bar reaches below it" for a zero-SE
            #   cell, which is the exact false-precision claim `_fmt` exists to stop.
            verdict, _z, lo0, _hi0 = classify_excess(
                z0["balanced"], z0["se"], 0.5, Z95)
            if not z0["se"] > 0:
                qual = "and carries NO bar, because it carries no information"
            elif z0["balanced"] <= 0.5:
                qual = "sits at or below the chance line"
            else:
                qual = (f"whose 95% bar {'reaches' if lo0 <= 0.5 else '⚠⚠ does NOT reach'}"
                        f" below it")
            print(f"\n- **{cond}:** A3 first reaches 0.75 accuracy at "
                  + (f"**f = {hit:g}**" if hit is not None
                     else "**no declared f** (never reaches 0.75)")
                  + f"; at f = 0 it sits at {_fmt(z0['balanced'], z0['se'])} "
                    f"against a chance line of 0.5000, {qual} "
                    f"(SE {z0['se']:.4f}, n = {z0['n_per_arm']} per arm).")
        print()

    if not res:
        # ⚠ --controls-only: no trials, so no bound check ran and there is
        #   nothing to fail on. `viol` is not bound this far up the function.
        return table, 0

    # ---------------- A3
    print("### A3 — rarity-weighted genotype linkage, re-posed as 1-of-2\n")
    print(f"Metrics are rows; `{knob}` is the column. Balanced accuracy, ties split.\n")
    print(f"⚠ Each cell is `attacker_accuracy` [`acc_ci95`] -- the point estimate and "
          f"its ERROR BAR, {Z95:.2f} x SE either side, clipped to [0, 1]. The bar is "
          f"what makes an excess over the line checkable rather than asserted; where a "
          f"cell sits above the line, read whether its lower limit falls below it. "
          f"⚠ Two cells carry no bar, because SE is exactly 0 and they are opposites: "
          f"`DEGENERATE` is a CONSTANT ANSWER at 0.5000 carrying no information, and "
          f"`EXACT` is PERFECT SEPARATION at 1.0000 carrying all of it.\n")
    hdr = "| metric | " + " | ".join(str(t) for t in taus) + " |"
    print(hdr); print("|" + "---|" * (len(taus) + 1))
    for cond in CONDS:
        cells = []
        for tau in taus:
            rows = [r["conditions"][cond] for r in res if r["tau"] == tau]
            lp = np.array([c["a3_under_p"] for c in rows])
            lq = np.array([c["a3_under_q"] for c in rows])
            bal, se, tpr, tnr = _acc(lp, lq)
            table["a3"].setdefault(cond, {})[tau] = {
                "balanced": bal, "se": se, "tpr": tpr, "tnr": tnr,
                "n_per_arm": len(lp)}
            cells.append(_fmt(bal, se))
        print(f"| {cond} balanced accuracy * | " + " | ".join(cells) + " |")
    bound_row()

    # ---------------- A4
    print("\n### A4 — segment-restricted likelihood ratio, maximised over blocks\n")
    print("⚠ EVERY declared block size is reported. The MAX over blocks is an "
          "optimistically biased estimate by construction -- maximising over `B` "
          "noisy accuracies -- so the Bonferroni-corrected comparison below is the "
          "only claimable one, and the raw max is shown beside it to make the size "
          "of that bias visible.\n")
    # ⚠⚠ TWO PASSES, and the split is load-bearing. The SIMULTANEOUS interval needs
    # `z_a4`, which needs `n_cells`, which is only known once every block-cell has
    # been counted -- so nothing can be printed until the whole table is computed.
    # ⚠⚠ AND THE SOURCE MAY BE THE STORED TABLE, NOT `res`. The persisted trials have
    # every `a4_*` key STRIPPED (~8M numbers), so a `--from-json` re-table CANNOT
    # recompute these from raw statistics -- it reads the per-cell summary the first
    # run stored. Before 2026-10-29 this branch did not exist and `--from-json` raised
    # KeyError on `a4_under_p`, which is why the A4 half of a table could never be
    # re-emitted without redrawing 7,000 releases.
    have_raw = bool(res) and "a4_under_p" in next(iter(res[0]["conditions"].values()))
    if not have_raw and not prior_table:
        raise SystemExit(
            "A4 needs either raw per-block statistics in `res` or a stored table to "
            "re-read. The trials in a persisted JSON carry neither, so pass the whole "
            "JSON's `table` as prior_table.")

    def _cell(d: dict, tau) -> dict:
        """A JSON round-trip turns the float column key into a string."""
        return d[tau] if tau in d else d[str(tau)]

    n_cells = 0
    for size in BLOCK_SIZES:
        for cond in CONDS:
            for tau in taus:
                if have_raw:
                    rows = [r["conditions"][cond] for r in res if r["tau"] == tau]
                    lp = np.array([c["a4_under_p"][str(size)] for c in rows])  # (n, B)
                    lq = np.array([c["a4_under_q"][str(size)] for c in rows])
                    nb = lp.shape[1]
                    accs = np.array([_acc(lp[:, b], lq[:, b])[0] for b in range(nb)])
                    ses = np.array([_acc(lp[:, b], lq[:, b])[1] for b in range(nb)])
                    b_max = int(accs.argmax())
                    cell = {
                        "n_blocks": nb, "max_balanced": float(accs[b_max]),
                        "argmax_block": b_max, "se_at_max": float(ses[b_max]),
                        "mean_balanced": float(accs.mean()),
                        "n_per_arm": int(lp.shape[0]),
                        # ⚠ the per-block accuracies are KEPT so the table can be
                        # rebuilt without redrawing 7,000 releases.
                        "block_balanced": [round(float(x), 6) for x in accs]}
                else:
                    cell = dict(_cell(prior_table["a4"][str(size)][cond], tau))
                table["a4"].setdefault(str(size), {}).setdefault(cond, {})[tau] = cell
                n_cells += cell["n_blocks"]

    z_a4_sel = _z_for_one_sided(0.05 / n_cells)
    print(f"⚠ Two intervals per cell, because the statistic is a MAXIMUM. "
          f"`acc_ci95` (square brackets) is {Z95:.2f} x SE at the argmax block and is "
          f"ANTI-CONSERVATIVE for a max -- it ignores the selection over B blocks. "
          f"`(simultaneous)` is {z_a4_sel:.2f} x SE, Bonferroni over all "
          f"{n_cells:,} block-cells, and is the interval the bound check below uses. "
          f"The gap between them IS the selection penalty.\n")
    print(hdr); print("|" + "---|" * (len(taus) + 1))
    for size in BLOCK_SIZES:
        for cond in CONDS:
            cells = []
            for tau in taus:
                c = table["a4"][str(size)][cond][tau]
                if not c["se_at_max"] > 0:
                    cells.append(_fmt(c["max_balanced"], 0.0))
                else:
                    slo, shi = _interval(c["max_balanced"], c["se_at_max"], z_a4_sel)
                    cells.append(f"{_fmt(c['max_balanced'], c['se_at_max'])} "
                                 f"({slo:.4f}, {shi:.4f})")
            print(f"| block {size}: {cond} MAX over blocks | " +
                  " | ".join(cells) + " |")
    bound_row()

    # ---------------- the bound check, Bonferroni over every cell examined
    a3_cells = len(CONDS) * len(taus)
    z_a3 = _z_for_one_sided(0.05 / a3_cells)
    z_a4 = _z_for_one_sided(0.05 / n_cells)
    _bname = "`acc_ceiling`" if is_ceil else "`p_succ_bound`"
    _skipped = "" if is_ceil else (
        " ⚠ SKIPPED where the point is an `eta_tau` (`p_succ_bound` is 1.0 to float64 "
        "there and cannot fail; compare against the exact `acc_ceiling` in "
        "logs/B_eta_sweep/).")
    print(f"\n**THE BOUND CHECK.**{_skipped} One-sided on the excess over {_bname}, "
          f"Bonferroni-corrected across the cells actually examined: A3 over "
          f"{a3_cells} cells (z_crit = {z_a3:.2f}), A4 over {n_cells:,} "
          f"block-cells (z_crit = {z_a4:.2f}).\n")
    if is_ceil:
        print(f"⚠⚠ THE LINE IS THE COLUMN LABEL, not `p_succ_bound`. In this mode each "
              f"target ran at its own `eta_tau` solved to land on that label, so the "
              f"label IS the measured `acc_ceiling` and IS the bound. Handing the label "
              f"to `p_succ_bound` instead would test against (1 + label)/2 -- 0.7500 at "
              f"the 0.50 column up to 0.9000 at 0.80 -- which is the defect fixed on "
              f"2026-10-29.\n")

    # ⚠⚠ THREE OUTCOMES, not two, and the middle one is new (Dylan, 2026-10-29).
    #   viol        -- significant at the Bonferroni-corrected z. THE CLAIMABLE test.
    #   bar_clears  -- the 95% ERROR BAR does not reach below the line. Not claimable
    #                  on its own (see the multiplicity note printed below), but it is
    #                  the case the reader must be shown rather than have absorbed by
    #                  a corrected threshold.
    #   reported    -- above the line by point estimate, bar reaches below. Noise.
    viol, bar_clears, reported, degen = [], [], [], []
    n_examined = 0

    def _judge(what: str, bal: float, se: float, t, z_crit: float, extra: str = "") -> None:
        nonlocal n_examined
        bound = bound_at(t)
        n_examined += 1
        verdict, z, lo, hi = classify_excess(bal, se, bound, z_crit)
        if verdict == "within":
            return
        if verdict == "degenerate":
            degen.append(f"{what} at {knob}={t}: {bal:.4f} vs bound {bound:.4f}, "
                         f"se = 0 -- DEGENERATE, no evidence either way" + extra)
            return
        msg = (f"{what} at {knob}={t}: {bal:.4f} [{lo:.4f}, {hi:.4f}] vs bound "
               f"{bound:.4f} (+{bal - bound:.4f}, z = {z:+.2f})" + extra)
        {"viol": viol, "bar_clears": bar_clears, "reported": reported}[verdict].append(msg)

    for cond in CONDS:
        for tau in taus:
            if not is_ceil and tau >= 1.0:
                continue
            c = table["a3"][cond][tau]
            _judge(f"A3 {cond}", c["balanced"], c["se"], tau, z_a3)
    for size in BLOCK_SIZES:
        for cond in CONDS:
            for tau in taus:
                if not is_ceil and tau >= 1.0:
                    continue
                c = table["a4"][str(size)][cond][tau]
                # ⚠ the SE is the argmax block's, UNADJUSTED for the selection; the
                # correction lives in z_a4, so the 95% bar here is anti-conservative.
                _judge(f"A4 block {size} {cond}", c["max_balanced"], c["se_at_max"],
                       tau, z_a4,
                       extra=f", block {c['argmax_block']} of {c['n_blocks']}")

    exp_false = 0.05 * n_examined
    if reported:
        print(f"⚠ {len(reported)} cell(s) exceed {_bname} by POINT ESTIMATE, each with "
              f"a 95% bar REACHING BELOW the line. Reported, not claimed:")
        for m in reported[:12]:
            print(f"    {m}")
        if len(reported) > 12:
            print(f"    ... and {len(reported) - 12} more (see the JSON)")
    if degen:
        print(f"\n⚠⚠ {len(degen)} cell(s) sit above {_bname} with se = 0 -- a CONSTANT "
              f"ANSWER, which is neither a violation nor a null:")
        for m in degen:
            print(f"    {m}")
    if bar_clears:
        print(f"\n⚠⚠ {len(bar_clears)} cell(s) exceed {_bname} and their 95% bar does "
              f"NOT reach below the line. Not significant after correction, so not "
              f"claimable -- but stated here rather than absorbed, because this is the "
              f"case an uncorrected reader would call a violation:")
        for m in bar_clears:
            print(f"    {m}")
    # ⚠⚠ THE VERDICT GOES IN THE JSON, not only in the prose. A re-table, a figure or
    # a later audit should not have to re-derive which cells were above the line and
    # under which bound -- that re-derivation is how the wrong line went unnoticed.
    table["bound_check"] = {
        "knob": knob, "bound_is_column_label": is_ceil, "n_examined": n_examined,
        "z_crit_a3": z_a3, "z_crit_a4": z_a4,
        "bar_half_width_z": Z95,
        "n_viol": len(viol), "n_bar_clears": len(bar_clears),
        "n_reported": len(reported), "n_degenerate": len(degen),
        "viol": list(viol), "bar_clears": list(bar_clears),
        "reported": list(reported), "degenerate": list(degen)}
    # ⚠⚠ THE 5% RATE APPLIES TO A3's CELLS ONLY. A4's are MAXIMA over hundreds of
    # blocks, and the chance that a max over B noisy blocks lands above the line is far
    # higher than 5% -- it is what `z_a4` exists to correct for. Quoting 0.05 x 84 would
    # make a statement about A4 that cannot be falsified by any A4 result.
    print(f"\n⚠ MULTIPLICITY. {n_examined:,} cells were examined. Of these, "
          f"{a3_cells} are A3 SINGLE estimates, where a one-sided 5% level puts about "
          f"{0.05 * a3_cells:,.1f} above the line BY CHANCE even if the mechanism is "
          f"exact. The other {n_examined - a3_cells} are A4 MAXIMA over hundreds of "
          f"blocks each ({n_cells:,} block-cells in total), whose chance of landing "
          f"above the line is far higher than 5% and is NOT quoted here -- that is "
          f"precisely what the z_crit of {z_a4:.2f} corrects for. A bare count of 95% "
          f"exceedances is not evidence in either attack; the corrected z is.")
    if viol:
        print(f"\n⚠⚠⚠ {len(viol)} cell(s) exceed {_bname} SIGNIFICANTLY. A "
              f"measured accuracy above the bound falsifies the IMPLEMENTATION, not "
              f"the theorem:")
        for m in viol:
            print(f"    {m}")
    else:
        print(f"\nno condition exceeded {_bname} significantly, in either attack")

    # ---------------- what the nulls do and do not mean
    print()
    if not (a3_ok and a4_ok):
        print("⚠⚠⚠ **A POSITIVE CONTROL FAILED, SO THE NULLS ABOVE ARE NOT PRIVACY "
              "EVIDENCE.** Standing constraint 7: an attacker that fails is a lower "
              "bound on distinguishability, and an attacker that fails WHEN HANDED "
              "THE ANSWER is a broken instrument that bounds nothing. Fix the "
              "instrument before quoting a single cell above.")
    else:
        print("✅ Both positive controls pass, so the nulls are lower bounds on "
              "distinguishability rather than broken instruments. ⚠ They remain LOWER "
              "bounds: the upper bound is `acc_ceiling`, which is proof-derived.")
    n_cell = table["a3"][CONDS[0]][taus[0]]["n_per_arm"]
    assert n_cell == n_per * len(tids), (
        f"per-cell n is {n_cell} but {n_per} draws x {len(tids)} targets is "
        f"{n_per * len(tids)}; the printed delta would describe a different experiment")
    se0 = 0.879 / math.sqrt(n_cell)
    print(f"⚠ Minimum detectable advantage at {n_cell} trials per arm per cell "
          f"({n_per} draws x {len(tids)} targets) is {se0:.4f} "
          f"UNCORRECTED. " + (
              (lambda unres: (
                  f"A cell resolves only where the bound exceeds {0.5 + se0:.4f}, so "
                  f"the {', '.join(f'{u:.2f}' for u in unres)} column(s) do NOT "
                  f"resolve below it in either attack"
              ) if unres else
                  f"Every column's bound exceeds {0.5 + se0:.4f}, so every one resolves"
              )([float(t) for t in taus if float(t) <= 0.5 + se0])
              if is_ceil else
              "so tau = 0.05 and 0.10 do not resolve below the bound in either attack")
          + f"; those rows are reported, not claimed, and the PROOF carries them. "
            f"A4's corrected threshold is {z_a4/1.645:.2f}x wider still.")
    print("⚠ D2's decoy pool is the other prepared targets, not the 3,202-sample "
          "panel, so it is a LOWER BOUND on decoy difficulty and flatters us.")
    print("⚠ The three conditions share each release and must not be pooled.")
    print(f"⚠ A3 scores on {a3_sites(load_state())['n_bi']:,} bi-allelic cohort sites "
          f"of {len(load_state()['util'].beta):,} chain positions; the rest carry "
          f"weight the attack cannot see, which flatters us.")
    return table, len(viol)


def _cap(plan, n):
    """Lower a worker plan to `n`. ⚠ May only LOWER: the plan's bound is memory or
    cores, not a preference, so raising it is how a run gets OOM-killed at hour three.
    See the note beside `arm2_attack`'s --workers."""
    if not n:
        return plan
    if n > plan.workers:
        raise SystemExit(f"--workers {n} exceeds the planned {plan.workers} "
                         f"({plan.reason}); that bound is memory or cores")
    return R.Plan(workers=n, per_job_bytes=plan.per_job_bytes,
                  headroom_bytes=plan.headroom_bytes,
                  reason=f"CAPPED by --workers {n} (planned {plan.workers}: "
                         f"{plan.reason})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    # ⚠ 25 draws x 10 targets = 250 trials per (tau, condition) per hypothesis,
    # which is A1's REALISED design. --draws 250 would give 2,500 per cell and
    # 70,000 releases; that overshoot was launched once and killed.
    ap.add_argument("--draws", type=int, default=25)
    ap.add_argument("--targets", type=int, default=10)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--controls-only", action="store_true")
    ap.add_argument("--workers", type=int, default=0,
                    help="cap the worker pool (0 = the memory/core plan). May only LOWER\n"
                         "the plan; use it to share the box with another run.")
    ap.add_argument("--weights", choices=("uniform", "inv_support"), default="uniform",
                    help="beta_v weighting; 'uniform' is the released tool's choice as "
                         "of 2026-10-24. Must match the weighting of every metric it "
                         "shares a table with (constraint 5).")
    ap.add_argument("--from-json", type=Path, default=None)
    # the `tau` index in the seed is its place in SEED_GRID, NOT in `--taus`, so a
    # subset run draws the same streams as the full grid (see arm2_attack.SEED_GRID)
    ap.add_argument("--taus", type=float, nargs="+", default=list(TAUS))
    ap.add_argument("--ceiling-grid", type=Path, default=None,
                    help="a `_verified.json` from scripts/arm2_chr21_kl.py --verify. "
                         "Each target runs at ITS OWN `eta_tau` with `acc_ceiling` as "
                         "the shared column -- constraint 7 matches on the ceiling, "
                         "never on eta.")
    ap.add_argument("--etas", type=float, nargs="+", default=None,
                    help="drive the tilt by `eta_tau` directly (the knob since "
                         "2026-10-24). ⚠ p_succ_bound is vacuous at these points; "
                         "the meaningful bound is the exact acc_ceiling. Must come "
                         "from arm2_attack.ETA_GRID so each cell has a fixed seed "
                         "index disjoint from the tau grid's.")
    a = ap.parse_args(argv)

    # ⚠ BEFORE THE FIRST `load_state()`. `load_state` lives in A1's module, so the
    # global that it reads is A1's, not this module's -- rebinding a local name here
    # would leave the attack running at the stale default. See the 2026-10-25 defect
    # note beside `arm2_attack.WEIGHTS`.
    A1.WEIGHTS = a.weights

    off = [t for t in a.taus if t not in SEED_GRID]
    if off:
        raise SystemExit(f"tau {off} is in neither the pinned grid nor "
                         f"EXPLORATORY_TAUS, so it has no fixed seed index")
    if any(t in EXPLORATORY_TAUS for t in a.taus):
        print(f"⚠⚠ EXPLORATORY `tau` {[t for t in a.taus if t in EXPLORATORY_TAUS]}: "
              f"NOT in the signed grid {list(TAUS)}. Report apart from the "
              f"pre-registered result.\n")

    rev = store.git_rev(REPO)
    prereg = "docs/plan.md direction block, 2026-10-24"
    st = load_state()
    tids = sorted(st["paths"])[: a.targets]
    a3_sites(st)

    if a.from_json:
        d = json.loads(a.from_json.read_text())
        # ⚠⚠ THE MODE MUST COME BACK WITH THE DATA. A ceiling run's column labels
        # (0.50 .. 0.80) are legal `tau` values too, so a re-report that guesses will
        # happily test a ceiling run against `p_succ_bound(label)` and move the bound
        # line by up to 0.25 without raising anything. That is how the wrong line
        # shipped in `a34_ceil.md`; refusing to guess is the fix.
        saved_knob = d.get("knob")
        saved_ceil = d.get("ceiling_etas")
        if saved_knob is None:
            if not a.ceiling_grid:
                raise SystemExit(
                    f"{a.from_json} records no `knob`, so its mode is not recoverable "
                    f"from the file. Its labels {list(d['taus'])} are legal in BOTH "
                    f"modes and the two bound lines differ by up to 0.25, so this will "
                    f"not be guessed. Pass --ceiling-grid <verified fit json> to "
                    f"re-report it as a CEILING run.")
            saved_ceil = ceiling_grid(a.ceiling_grid, d["targets"])
            saved_knob = "acc_ceiling"
            print(f"⚠ mode supplied from --ceiling-grid: this JSON predates the `knob` "
                  f"field and is being re-reported as a CEILING run.")
        if saved_knob == "acc_ceiling" and not saved_ceil:
            raise SystemExit(f"{a.from_json} says knob=acc_ceiling but carries no "
                             f"`ceiling_etas`, so the per-target eta cannot be shown.")
        _t, n_viol = _report(
            d["trials"], d["controls"], tuple(d["taus"]),
            d.get("methodology", d.get("prereg", "unrecorded")),
            d["git_rev"], d["targets"], d["n_per"], pw=d.get("power"),
            ceilings=saved_ceil if saved_knob == "acc_ceiling" else None,
            prior_table=d.get("table"))
        # ⚠⚠ a re-table must carry the same verdict as the run. Before 2026-10-29 this
        # driver returned 0 unconditionally, so a SIGNIFICANT excess exited clean and
        # the queue's `fail=0` completeness line recorded it as a success.
        return 1 if n_viol else 0

    print(f"methodology: {prereg} (git_rev {rev})")
    print(f"A3 scores on {_A3['n_bi']:,} bi-allelic cohort sites of "
          f"{len(st['util'].beta):,} chain positions.")
    print("computing the positive controls first -- a failed control voids the "
          "nulls...", flush=True)
    ctl = controls(st, tids)
    print(f"  A3 unobfuscated link rate: {ctl['a3']['rate']:.4f} "
          f"({ctl['a3']['n_correct']}/{ctl['a3']['n']})")
    print(f"  A4 planted-block detection: {ctl['a4']['rate']:.4f} "
          f"({ctl['a4']['n_correct']}/{ctl['a4']['n']})", flush=True)
    pw_jobs = [(t, POWER_SEED_BASE + i * 100 + k)
               for i, t in enumerate(tids) for k in range(POWER_SEEDS)]
    pplan = _cap(R.plan_workers(per_job_bytes=600_000_000, n_jobs=len(pw_jobs)), a.workers)
    print(f"\nA3 power curve: {len(pw_jobs)} jobs x 4 untilted releases "
          f"= {4*len(pw_jobs)} releases ({pplan.workers} workers)", flush=True)
    pw = R.run_jobs(_power_trial, pw_jobs, workers=pplan.workers, label="power",
                    progress_every=max(1, len(pw_jobs) // 5))

    if a.controls_only:
        _report([], ctl, (), prereg, rev, tids, a.draws, pw=pw)
        if a.json:
            a.json.write_text(json.dumps({"controls": ctl, "power": pw}, indent=1) + "\n")
            print(f"JSON -> {a.json}")
        return 0

    # ⚠⚠⚠ OCCURRENCE SIX OF THE COMMON-RANDOM-NUMBERS DEFECT, FOUND 2026-10-28 AND
    #     FIXED HERE. The seed was ARITHMETIC:
    #         SEED_BASE + i * 10_000 + (ETA_SEED_OFFSET + ETA_GRID.index(eta)) * 1_000 + k
    #     with `i` the TARGET index and `k` the draw. Appending the reported grid to
    #     ETA_GRID on 2026-10-28 pushed `index(eta)` to 10 and 11, and at that point
    #     `i * 10_000` and `index * 1_000` ALIAS: (target i, eta 378.5, draw k) lands on
    #     exactly the seed of (target i+1, eta 0.0, draw k). MEASURED on the shipped
    #     configuration: 225 colliding seeds, 450 of 1,750 cells = 25.71% sharing a
    #     release stream with another cell.
    #     ⚠ The append-only discipline on ETA_GRID protected the INDICES; it could not
    #       protect the ARITHMETIC that consumes them. That is the lesson.
    #     ⚠ `arm2_attack.py` was already immune: it was converted to a SeedSequence over
    #       the whole key in 2026-10-25. This driver was not, and inherited the gap.
    #     The fix is the same one: a SeedSequence over the full cell key, which cannot
    #     alias whatever the index ranges become.
    #     ⚠ The `34` tag keeps this driver's streams disjoint from `arm2_attack.py`'s,
    #       which uses the same SEED_BASE and the same (target, point, draw) key. Without
    #       it the two attacks would silently score the SAME releases.
    if a.ceiling_grid:
        if a.etas:
            raise SystemExit("--ceiling-grid and --etas are two different knobs. Pass "
                             "one: --etas puts every target on the SAME eta (and so on "
                             "different ceilings), --ceiling-grid on the same CEILING.")
        cg = ceiling_grid(a.ceiling_grid, tids)
        a.ceiling_etas = cg
        a.taus = [round(c, 2) for c in REPORTED_CEILINGS]
        # ⚠ a PRE-RUN confirmation only -- the per-label eta ranges are printed by
        # `_report`, so that a `--from-json` re-table carries them too. Before
        # 2026-10-29 they lived only here and a re-tabled md lost its own grid.
        print(f"⚠⚠ CEILING MODE: {len(tids)} targets, each at ITS OWN `eta_tau`, from "
              f"{a.ceiling_grid}. The column is `acc_ceiling` and IS the bound.\n")
        jobs = [(t, round(acc, 2), cg[t][f"{acc:.2f}"],
                 _cell_seed(i, CEIL_SEED_OFFSET + REPORTED_CEILINGS.index(acc), k))
                for i, t in enumerate(tids)
                for acc in REPORTED_CEILINGS
                for k in range(a.draws)]
    elif a.etas:
        off = [e for e in a.etas if e not in ETA_GRID]
        if off:
            raise SystemExit(f"eta {off} is not in arm2_attack.ETA_GRID "
                             f"{ETA_GRID}, so it has no fixed seed index")
        print(f"⚠⚠ ETA MODE: {list(a.etas)}. p_succ_bound is 1.0 to float64 at "
              f"every one of these, so no bound comparison is printed; compare "
              f"against the exact acc_ceiling in logs/B_eta_sweep/.\n")
        a.taus = list(a.etas)
        jobs = [(t, eta, eta,
                 _cell_seed(i, ETA_SEED_OFFSET + ETA_GRID.index(eta), k))
                for i, t in enumerate(tids)
                for eta in a.etas
                for k in range(a.draws)]
    else:
        jobs = [(t, tau, C.eta_tau(tau),
                 _cell_seed(i, SEED_GRID.index(tau), k))
                for i, t in enumerate(tids)
                for tau in a.taus
                for k in range(a.draws)]
    plan = _cap(R.plan_workers(per_job_bytes=600_000_000, n_jobs=len(jobs)), a.workers)
    print(f"\n{len(jobs):,} trials, each drawing 1 release under Q_p and 3 under "
          f"Q_q (one per condition) = {4*len(jobs):,} releases")
    print(f"{plan}\n", flush=True)
    res = R.run_jobs(_one_trial, jobs, workers=plan.workers, label="trials",
                     progress_every=max(1, len(jobs) // 20))

    table, n_viol = _report(res, ctl, tuple(a.taus), prereg, rev, tids, a.draws,
                            pw=pw, ceilings=getattr(a, "ceiling_etas", None))
    if a.json:
        # ⚠ the raw per-block statistics are DROPPED from the JSON: 1,750 trials x 3
        # conditions x 2 hypotheses x 743 blocks is ~8M numbers. `table["a4"]`
        # retains every per-block ACCURACY, which determines every reported quantity,
        # so the table can be rebuilt with --from-json without redrawing a release.
        lean = [{k: v for k, v in r.items() if k != "conditions"}
                | {"conditions": {c: {k: v for k, v in cell.items()
                                      if not k.startswith("a4_")}
                                  for c, cell in r["conditions"].items()}}
                for r in res]
        a.json.write_text(json.dumps(
            {"methodology": prereg, "git_rev": rev, "taus": list(a.taus), "targets": tids,
             # ⚠ THE KNOB IS PART OF THE DATA. Without it `--from-json` cannot tell a
             # ceiling run from a tau run and the bound line becomes a guess.
             "knob": ("acc_ceiling" if getattr(a, "ceiling_etas", None)
                      else ("eta_tau" if getattr(a, "etas", None) else "tau")),
             "ceiling_etas": getattr(a, "ceiling_etas", None),
             "ceiling_grid_source": (str(a.ceiling_grid) if a.ceiling_grid else None),
             "n_per": a.draws, "block_sizes": list(BLOCK_SIZES),
             "n_bi_sites": _A3["n_bi"], "controls": ctl, "power": pw, "table": table,
             "trials": lean,
             "note": ("raw per-block A4 statistics are not persisted; "
                      "table.a4[size][cond][tau].block_balanced holds every "
                      "per-block accuracy, which is what the report uses")},
            indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    # ⚠⚠ A SIGNIFICANT EXCESS MUST FAIL THE RUN. `arm2_attack.py` has returned
    # `1 if violations else 0` since it was written; this driver returned 0
    # unconditionally until 2026-10-29, so a falsifying A3 or A4 cell would have
    # exited clean and every queue that reads `$?` would have logged `fail=0`.
    return 1 if n_viol else 0


if __name__ == "__main__":
    sys.exit(main())
