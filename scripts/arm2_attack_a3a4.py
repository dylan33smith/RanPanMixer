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
from arm2_attack import (CHROM, DECOY_SEED, ETA_GRID, ETA_SEED_OFFSET,  # noqa: E402
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
            tids: list, n_per: int, pw: list | None = None) -> dict:
    table: dict = {"a3": {}, "a4": {}}
    print(f"methodology: {prereg} (git_rev {rev})\n")
    print(f"ARM 2. {len(tids)} targets x {len(taus)} tau x {n_per} draws = "
          f"{n_per * len(tids)} trials per (tau, condition) PER HYPOTHESIS -- matching "
          f"A1's realised design of 250. chr21, scope chain, beta_share = 1.\n")

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
                bal = _acc(lp, lq)[0]
                table.setdefault("power", {}).setdefault(cond, {})[str(_fr)] = bal
                cells.append(f"{bal:.4f}")
            print(f"| {cond} | " + " | ".join(cells) + " |")
        # the smallest declared leak A3 resolves at all
        for cond in CONDS:
            curve = [(x, table["power"][cond][str(x)]) for x in LEAK_F]
            hit = next((x for x, b in curve if b >= 0.75), None)
            print(f"\n- **{cond}:** A3 first reaches 0.75 accuracy at "
                  + (f"**f = {hit:g}**" if hit is not None
                     else "**no declared f** (never reaches 0.75)")
                  + f"; at f = 0 it sits at {table['power'][cond]['0.0']:.4f}.")
        print()

    if not res:
        return table

    # ---------------- A3
    print("### A3 — rarity-weighted genotype linkage, re-posed as 1-of-2\n")
    print("Metrics are rows; `tau` is the column. Balanced accuracy, ties split.\n")
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
            cells.append(f"{bal:.4f}")
        print(f"| {cond} balanced accuracy * | " + " | ".join(cells) + " |")
    print("| `p_succ_bound` * | " + " | ".join(
        ("~1.0 (vacuous)" if t >= 1.0 else f"{C.p_succ_bound(t):.4f}")
        for t in taus) + " |")

    # ---------------- A4
    print("\n### A4 — segment-restricted likelihood ratio, maximised over blocks\n")
    print("⚠ EVERY declared block size is reported. The MAX over blocks is an "
          "optimistically biased estimate by construction -- maximising over `B` "
          "noisy accuracies -- so the Bonferroni-corrected comparison below is the "
          "only claimable one, and the raw max is shown beside it to make the size "
          "of that bias visible.\n")
    print(hdr); print("|" + "---|" * (len(taus) + 1))
    n_cells = 0
    for size in BLOCK_SIZES:
        for cond in CONDS:
            cells = []
            for tau in taus:
                rows = [r["conditions"][cond] for r in res if r["tau"] == tau]
                lp = np.array([c["a4_under_p"][str(size)] for c in rows])  # (n, B)
                lq = np.array([c["a4_under_q"][str(size)] for c in rows])
                nb = lp.shape[1]
                accs = np.array([_acc(lp[:, b], lq[:, b])[0] for b in range(nb)])
                ses = np.array([_acc(lp[:, b], lq[:, b])[1] for b in range(nb)])
                b_max = int(accs.argmax())
                table["a4"].setdefault(str(size), {}).setdefault(cond, {})[tau] = {
                    "n_blocks": nb, "max_balanced": float(accs[b_max]),
                    "argmax_block": b_max, "se_at_max": float(ses[b_max]),
                    "mean_balanced": float(accs.mean()),
                    "n_per_arm": int(lp.shape[0]),
                    # ⚠ the per-block accuracies are KEPT so the table can be rebuilt
                    # without redrawing 7,000 releases. The raw per-trial block
                    # statistics are NOT persisted -- they are ~8M numbers -- and these
                    # determine every reported A4 quantity exactly.
                    "block_balanced": [round(float(x), 6) for x in accs]}
                n_cells += nb
                cells.append(f"{accs[b_max]:.4f}")
            print(f"| block {size}: {cond} MAX over blocks | " +
                  " | ".join(cells) + " |")
    print("| `p_succ_bound` * | " + " | ".join(
        ("~1.0 (vacuous)" if t >= 1.0 else f"{C.p_succ_bound(t):.4f}")
        for t in taus) + " |")

    # ---------------- the bound check, Bonferroni over every cell examined
    a3_cells = len(CONDS) * len(taus)
    z_a3 = _z_for_one_sided(0.05 / a3_cells)
    z_a4 = _z_for_one_sided(0.05 / n_cells)
    print(f"\n**THE BOUND CHECK.** ⚠ SKIPPED where the point is an `eta_tau` (p_succ_bound is 1.0 to float64 there and cannot fail; compare against the exact `acc_ceiling` in logs/B_eta_sweep/). Otherwise one-sided on the excess over `p_succ_bound`, "
          f"Bonferroni-corrected across the cells actually examined: A3 over "
          f"{a3_cells} cells (z_crit = {z_a3:.2f}), A4 over {n_cells:,} "
          f"block-cells (z_crit = {z_a4:.2f}).\n")
    viol, reported = [], []
    for cond in CONDS:
        for tau in taus:
            c = table["a3"][cond][tau]
            if tau >= 1.0:
                continue
            bound = C.p_succ_bound(tau)
            if c["balanced"] <= bound:
                continue
            z = (c["balanced"] - bound) / c["se"] if c["se"] > 0 else float("inf")
            msg = (f"A3 {cond} at tau={tau}: {c['balanced']:.4f} vs bound "
                   f"{bound:.4f} (+{c['balanced']-bound:.4f}, z = {z:+.2f})")
            (viol if z > z_a3 else reported).append(msg)
    for size in BLOCK_SIZES:
        for cond in CONDS:
            for tau in taus:
                c = table["a4"][str(size)][cond][tau]
                # ⚠ an eta_tau point: p_succ_bound is 1.0 to float64 and cannot fail
                if tau >= 1.0:
                    continue
                bound = C.p_succ_bound(tau)
                if c["max_balanced"] <= bound:
                    continue
                se = c["se_at_max"]
                z = (c["max_balanced"] - bound) / se if se > 0 else float("inf")
                msg = (f"A4 block {size} {cond} at tau={tau}: "
                       f"{c['max_balanced']:.4f} vs bound {bound:.4f} "
                       f"(+{c['max_balanced']-bound:.4f}, z = {z:+.2f}, "
                       f"block {c['argmax_block']} of {c['n_blocks']})")
                (viol if z > z_a4 else reported).append(msg)
    if reported:
        print(f"⚠ {len(reported)} cell(s) exceed `p_succ_bound` by POINT ESTIMATE but "
              f"NOT significantly. Reported, not claimed:")
        for m in reported[:12]:
            print(f"    {m}")
        if len(reported) > 12:
            print(f"    ... and {len(reported) - 12} more (see the JSON)")
    if viol:
        print(f"\n⚠⚠⚠ {len(viol)} cell(s) exceed `p_succ_bound` SIGNIFICANTLY. A "
              f"measured accuracy above the bound falsifies the IMPLEMENTATION, not "
              f"the theorem:")
        for m in viol:
            print(f"    {m}")
    else:
        print("\nno condition exceeded `p_succ_bound` significantly, in either attack")

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
          f"UNCORRECTED, so tau = 0.05 and 0.10 do not resolve below the bound in "
          f"either attack; those rows are reported, not claimed, and the PROOF "
          f"carries them. A4's corrected threshold is {z_a4/1.645:.2f}x wider still.")
    print("⚠ D2's decoy pool is the other prepared targets, not the 3,202-sample "
          "panel, so it is a LOWER BOUND on decoy difficulty and flatters us.")
    print("⚠ The three conditions share each release and must not be pooled.")
    print(f"⚠ A3 scores on {a3_sites(load_state())['n_bi']:,} bi-allelic cohort sites "
          f"of {len(load_state()['util'].beta):,} chain positions; the rest carry "
          f"weight the attack cannot see, which flatters us.")
    return table


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
        _report(d["trials"], d["controls"], tuple(d["taus"]), d.get("methodology", d.get("prereg", "unrecorded")),
                d["git_rev"], d["targets"], d["n_per"])
        return 0

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

    if a.etas:
        off = [e for e in a.etas if e not in ETA_GRID]
        if off:
            raise SystemExit(f"eta {off} is not in arm2_attack.ETA_GRID "
                             f"{ETA_GRID}, so it has no fixed seed index")
        print(f"⚠⚠ ETA MODE: {list(a.etas)}. p_succ_bound is 1.0 to float64 at "
              f"every one of these, so no bound comparison is printed; compare "
              f"against the exact acc_ceiling in logs/B_eta_sweep/.\n")
        a.taus = list(a.etas)
        jobs = [(t, eta, eta, SEED_BASE + i * 10_000
                 + (ETA_SEED_OFFSET + ETA_GRID.index(eta)) * 1_000 + k)
                for i, t in enumerate(tids)
                for eta in a.etas
                for k in range(a.draws)]
    else:
        jobs = [(t, tau, C.eta_tau(tau),
                 SEED_BASE + i * 10_000 + SEED_GRID.index(tau) * 1_000 + k)
                for i, t in enumerate(tids)
                for tau in a.taus
                for k in range(a.draws)]
    plan = _cap(R.plan_workers(per_job_bytes=600_000_000, n_jobs=len(jobs)), a.workers)
    print(f"\n{len(jobs):,} trials, each drawing 1 release under Q_p and 3 under "
          f"Q_q (one per condition) = {4*len(jobs):,} releases")
    print(f"{plan}\n", flush=True)
    res = R.run_jobs(_one_trial, jobs, workers=plan.workers, label="trials",
                     progress_every=max(1, len(jobs) // 20))

    table = _report(res, ctl, tuple(a.taus), prereg, rev, tids, a.draws, pw=pw)
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
             "n_per": a.draws, "block_sizes": list(BLOCK_SIZES),
             "n_bi_sites": _A3["n_bi"], "controls": ctl, "power": pw, "table": table,
             "trials": lean,
             "note": ("raw per-block A4 statistics are not persisted; "
                      "table.a4[size][cond][tau].block_balanced holds every "
                      "per-block accuracy, which is what the report uses")},
            indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
