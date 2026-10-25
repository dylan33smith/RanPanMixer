#!/usr/bin/env python
"""GENOME-WIDE `kl_exact` and `acc_ceiling`, by a route that does not cancel.

⚠⚠ WHY THE OBVIOUS ROUTE FAILS. `scripts/exact_kl_genome.py` sums the per-chromosome KL
computed the usual way and FAILS ITS OWN SIGN CHECK at genome scale: KL came back at
-1.96e-11 on chr19 at `tau` = 0.05, because `eta_tau * (E[u_p] - E[u_q])` = 6.768228233e-04
and `log Z_p - log Z_q` = 6.768228429e-04 agree to ten significant figures. The expected
genome-wide KL is ~3.4e-11 while the cancellation noise is ~2e-11: the signal sits at the
estimator's noise floor.

⚠⚠ THE CANCELLATION IS AN ARTEFACT OF HOW THE TERMS ARE GROUPED, NOT OF THE QUANTITY.
`ffbs.log_z_p`'s own docstring says its tolerance "covers accumulated floating-point error
in the difference of two large sums" -- each about -5.85e2 on chr21 and ~66x that
genome-wide. Two regroupings remove it:

  1. **`log Z_p - log Z_q` accumulated PER POSITION.** Both subtract the SAME
     `baseline_log_norm`, so it cancels EXACTLY and never enters the arithmetic:
         log Z_p - log Z_q = sum_t (log c_t^p - log c_t^q)
     Each term is a difference of two O(1e-3) numbers rather than of two O(1e2) sums, so
     its absolute error is ~1e-19 instead of ~1e-12.
  2. **`E[u_p] - E[u_q]` as ONE expectation of the DIFFERENCE.** The per-position factor
     is `match_p - match_q` in {-1, 0, +1} -- exact small integers -- so
         E_{Q_p}[u_p - u_q] = sum_t sum_k gamma_t(k) beta_t (match_p - match_q)
     has no cancellation at all, where differencing two separately-accumulated
     expectations of ~6.8e-4 does.

Both totals use `math.fsum`, which is exactly rounded, so the only error left is in the
inputs.

⚠ KL is ADDITIVE across chromosomes (the release factorises given the target, one budget
split by decision J shares), so the genome-wide KL is the sum of per-chromosome KLs, each
at its own `beta_share`. The Pinsker step is taken ONCE on the total: `sqrt` is concave,
so per-chromosome Pinsker would be wrong AND would flatter us.

⚠⚠ VALIDATION IS NOT OPTIONAL HERE. A new numerical route that reports a number the old
one could not is exactly the situation where a silent error is invisible. So `--validate`
reproduces the chr21-at-`beta_share`-1 figures, where the OLD route works and
`logs/B_exact_kl.md` has the answer, and refuses to report genome-wide numbers unless they
agree.

⚠ ARM 2: K = 88, nothing held out, D1_extremal decoy.

    python scripts/exact_kl_stable.py --validate
    python scripts/exact_kl_stable.py --targets 1 --workers 22
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

from cohort import loo                                     # noqa: E402
from mechanism import calibration as C                     # noqa: E402
from mechanism import ffbs, hmm, utility                   # noqa: E402

PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
AUTOSOMES = loo.AUTOSOMES


def kl_stable(m, hap, beta, p, q, tau: float | None = None, *,
              eta: float | None = None) -> dict:
    """KL(Q_p || Q_q) for one chromosome, grouped so nothing cancels.

    ⚠ DRIVEN BY `tau` OR BY `eta`, never both. The `eta` entry point exists because
    `tau` CANNOT express the high end of the tilt: `tanh(eta)` rounds to exactly 1.0
    in float64 beyond `eta` ~ 19, and `calibration.eta_tau` rightly refuses
    `tau` = 1. `scripts/eta_sweep.py` therefore parameterises by `eta` directly.
    Kept as ONE implementation with two entry points rather than a second copy,
    because two copies of a KL would drift and only one of them would be gated.
    """
    if (tau is None) == (eta is None):
        raise ValueError("pass exactly one of tau= or eta=")
    if eta is None:
        eta = C.eta_tau(tau)
    elif eta < 0.0:
        raise ValueError(f"eta must be non-negative, got {eta}")
    if eta == 0.0:
        # tau = 0: Q_p == Q_q == R_D definitionally, so the KL is exactly 0.
        return {"kl": 0.0, "d_logz": 0.0, "d_eu": 0.0}
    h = eta / 2.0
    # ⚠ `beta` may be a PER-HAPLOTYPE pair. `target_alt` weights each chain by the
    # positions where THAT haplotype carries a non-ref allele, and the two sets
    # differ, so sharing one array would weight each chain by the other's support.
    betas = beta if isinstance(beta, (tuple, list)) else (beta, beta)
    d_logz_terms, d_eu = [], 0.0
    for col in (0, 1):
        pc, qc = p[:, col], q[:, col]
        beta = betas[col]
        psi_p = utility.tilt_potentials(m.allowed, hap, pc, beta, h)
        psi_q = utility.tilt_potentials(m.allowed, hap, qc, beta, h)
        f_p = ffbs.forward(m, psi_p)
        f_q = ffbs.forward(m, psi_q)
        # (1) PER-POSITION difference. The baseline cancels exactly: it is the same
        #     constant inside both log Z values and never enters here.
        d_logz_terms.append(np.asarray(f_p.log_c) - np.asarray(f_q.log_c))
        # (2) ONE expectation of the DIFFERENCE. match_p - match_q is in {-1, 0, +1}.
        gam = ffbs.state_marginals(m, f_p, psi_p)
        mp = (hap == pc[None, :]) & (pc[None, :] >= 0)
        mq = (hap == qc[None, :]) & (qc[None, :] >= 0)
        diff = mp.astype(np.int8) - mq.astype(np.int8)
        d_eu += 0.5 * float((gam * diff.T).sum(axis=1) @ beta)
        del psi_p, psi_q, f_p, f_q, gam, mp, mq, diff
    d_logz = math.fsum(math.fsum(t.tolist()) for t in d_logz_terms)
    kl = eta * d_eu - d_logz
    return {"kl": kl, "d_logz": d_logz, "d_eu": d_eu}


def extremal(p: np.ndarray) -> np.ndarray:
    """D1_extremal, built exactly as the sampled attack builds it."""
    e = p.copy()
    called = e >= 0
    e[called] = np.where(e[called] == 0, 1, 0)
    return e


def knob_key(v: float) -> str:
    """The ONE spelling of a knob value used as a result-dict key, both sides.

    ⚠⚠ THE 2026-10-25 DEFECT. The producer (`one_chrom`) keyed `f"{v:g}"` and the
    consumer (`main`) asked for `str(v)`, so `eta` = 0.0 was WRITTEN '0' and READ
    '0.0'. `KeyError: '0.0'` after 3.4 hours of compute, with every per-chromosome
    result held in memory and the JSON write sitting downstream of the crash, so the
    whole run was lost. One helper called from both sides is what makes that drift
    impossible; a second spelling anywhere is the bug coming back.

    ⚠ `f"{v:g}"` was also LOSSY — 6 significant figures, so two grid points agreeing
    to 6 s.f. would silently collapse into ONE cell and one would overwrite the other.
    `repr(float)` is the shortest round-trippable form in Python 3, and the assertion
    makes the round-trip a gate rather than a hope. `src/releases.py:eta_dir` already
    refuses this for release DIRECTORIES; this refuses it for result KEYS.
    """
    k = repr(float(v))
    assert float(k) == float(v), f"knob key {k!r} does not round-trip to {v!r}"
    return k


def one_chrom(chrom: str, share: float | None, targets, taus, *,
              arm: str = "arm2", etas=None, weights: str = "inv_support") -> dict:
    """One chromosome's exact KL, for ARM 2 (K = 88) or ARM 1 (K = 86, leave-one-out).

    ⚠⚠ ARM 1 IS NOT ARM 2 WITH A FLAG. Its target IS a cohort member, so the
    guarantee holds ONLY with BOTH of the target's haplotypes out of the donor panel
    (standing constraint 1), which means a DIFFERENT chain, a different K and a
    different beta per target -- and `logs/B_loo_chain_survey.md` measured that the
    untraversable-block set is NOT monotone in the panel, so arm 1's chain can never
    be built by subsetting arm 2's. Everything is therefore rebuilt per target here.
    Arm 1's target path is its OWN two haplotype rows; no target VCF exists for it.
    """
    d = REPO / "data" / "cohort" / chrom
    points = [("eta", float(e)) for e in etas] if etas else [("tau", float(t)) for t in taus]
    out = {}
    if arm == "arm2":
        ic = np.load(d / "in_chain.npy").astype(bool)
        hap = np.load(d / "haplotypes.npy")[:, ic]
        m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                            np.load(d / "positions.npy")[ic], cohort_id=chrom)
        kw = {} if share is None else {"share": share, "chrom": chrom}
        beta = utility.beta_v(np.load(d / "support.npy"), in_chain=ic,
                              weights=weights, **kw).beta
        for t in targets:
            p = np.load(REPO / "data" / "targets" / t / chrom / "path.npy")[ic]
            q = extremal(p)
            out[t] = {knob_key(v): kl_stable(m, hap, beta, p, q,
                                          **({"eta": v} if k == "eta" else {"tau": v}))
                      for k, v in points}
        t_chain = int(ic.sum())
    else:
        hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
        ids = np.load(d / "haplotype_ids.npy", allow_pickle=True)
        idx = loo.sample_index(ids)
        t_chain = 0
        for t in targets:
            lc = loo.loo_chain(d, t, haplotypes=hap_all)
            ic = lc.in_chain
            panel = loo.masked_haplotypes(d, t, haplotypes=hap_all)[:, ic]
            m = hmm.from_arrays(panel, np.load(d / "genetic_pos.npy")[ic],
                                np.load(d / "positions.npy")[ic],
                                cohort_id=f"{chrom}_loo_{t}")
            kw = {} if share is None else {"share": share, "chrom": chrom}
            beta = utility.beta_v(lc.support, in_chain=ic,
                                  weights=weights, **kw).beta
            p = hap_all[list(idx[t])][:, ic].T          # arm 1's target is its own rows
            q = extremal(p)
            out[t] = {knob_key(v): kl_stable(m, panel, beta, p, q,
                                          **({"eta": v} if k == "eta" else {"tau": v}))
                      for k, v in points}
            t_chain = int(ic.sum())
            del panel, m, beta
    return {"chrom": chrom, "kl": out, "t_chain": t_chain, "arm": arm}


def _one(args, kw=None):
    try:
        return one_chrom(*args, **(kw or {}))
    except Exception as e:
        return {"chrom": args[0], "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def validate(taus) -> None:
    """Reproduce the chr21-at-`beta_share`-1 figures the OLD route already measured."""
    known = {"0.05": 6.894e-09, "0.1": 2.780e-08, "0.25": 1.801e-07,
             "0.5": 8.329e-07, "0.9": 5.983e-06, "0.99": 1.934e-05}
    tids = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                  if p.is_dir() and (p / "chr21" / "path.npy").exists())[:10]
    print(f"VALIDATION: chr21, `beta_share` = 1, K = 88, n = {len(tids)} targets, "
          f"D1_extremal — the configuration `logs/B_exact_kl.md` already measured.\n")
    r = one_chrom("chr21", None, tids, [float(t) for t in known])
    print(f"| `tau` | logged `kl_exact` (max) | STABLE ROUTE (max) | rel. diff |")
    print(f"|---|---|---|---|")
    worst = 0.0
    for k, v in known.items():
        got = max(r["kl"][t][knob_key(float(k))]["kl"] for t in tids)
        rel = abs(got - v) / v
        worst = max(worst, rel)
        print(f"| {k} | {v:.3e} | {got:.3e} | {rel:.2e} |")
    print(f"\nworst relative difference: {worst:.3e}")
    if worst > 0.01:
        raise SystemExit(
            f"the stable route does NOT reproduce the logged chr21 figures (worst "
            f"{worst:.3e}). It must agree where the old route works before any "
            f"genome-wide number from it is reported.")
    print(f"✅ the stable route reproduces the logged figures to {worst:.2e}. Its "
          f"genome-wide output may be reported.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--targets", type=int, default=1,
                    help="how many targets; arm1 draws from the 44 cohort members, "
                         "arm2 from the external test set")
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--etas", type=float, nargs="+", default=None,
                    help="drive by eta_tau instead of tau (the knob since 2026-10-24)")
    ap.add_argument("--arm", choices=("arm1", "arm2"), default="arm2",
                    help="arm1 = K = 86 per-target leave-one-out, target is its own "
                         "haplotype rows; arm2 = K = 88, external target")
    ap.add_argument("--weights", choices=("uniform", "inv_support"),
                    default="uniform")
    ap.add_argument("--workers", type=int, default=22)
    ap.add_argument("--resume", action="store_true",
                    help="reuse per-chromosome checkpoints beside --json; only\n"
                         "those written for THIS arm, weighting, grid and\n"
                         "target set are reused, the rest are recomputed")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.validate:
        # ⚠⚠ `validate()` reproduces logged chr21 ARM 2 / inv_support figures and
        # NOTHING ELSE, then prints "its genome-wide output may be reported". Under
        # --arm arm1 or --etas it was gating a configuration it never ran, which is a
        # gate that cannot fail — decoration, per CLAUDE.md. Refuse instead of
        # passing: a missing reference must be reported, never silently become a pass.
        if a.etas or a.arm != "arm2" or a.weights != "inv_support":
            raise SystemExit(
                f"--validate only has logged reference figures for arm2 / tau / "
                f"inv_support (chr21, beta_share = 1). You asked for arm={a.arm}, "
                f"weights={a.weights}, knob={'eta_tau' if a.etas else 'tau'}, which "
                f"has no anchor on disk — so a PASS here would mean nothing. Run "
                f"--validate bare to check the stable route, then run your real "
                f"configuration without it.")
        validate(a.taus)
        return 0

    T = {c: int(np.load(REPO / "data" / "cohort" / c / "in_chain.npy").sum())
         for c in AUTOSOMES}
    tot = sum(T.values())
    sh = {c: T[c] / tot for c in AUTOSOMES}
    assert abs(sum(sh.values()) - 1.0) < 1e-12
    if a.arm == "arm1":
        # ARM 1's targets are the 44 COHORT MEMBERS, not the external test set --
        # its target IS a cohort member and its path is its own haplotype rows.
        ids = np.load(REPO / "data" / "cohort" / "chr1" / "haplotype_ids.npy",
                      allow_pickle=True)
        tids = sorted(loo.sample_index(ids))[: a.targets]
    else:
        tids = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                      if p.is_dir() and (p / "chr1" / "path.npy").exists())[: a.targets]

    knob = "eta_tau" if a.etas else "tau"
    # ⚠⚠ ONE grid object, built once, used everywhere. The 2026-10-25 failure was not
    # one bug: `a.taus` was read at the aggregation loop, at the `p_succ_bound` call,
    # at the `_report` call and in the JSON while the run was driven by `--etas`, so
    # the aggregator walked a seven-point tau grid the workers had never evaluated.
    # Five crash and mislabel sites, one cause. Two names for one thing is the cause.
    grid = [float(e) for e in a.etas] if a.etas else [float(t) for t in a.taus]
    if len({knob_key(v) for v in grid}) != len(grid):
        raise SystemExit(f"duplicate knob values in {grid}: two grid points would "
                         f"share one result key and silently overwrite each other")
    print(f"GENOME-WIDE `kl_exact`, STABLE ROUTE. "
          f"{'ARM 1, K = 86, per-target leave-one-out' if a.arm == 'arm1' else 'ARM 2, K = 88'}, "
          f"22 autosomes, T = {tot:,}, n = {len(tids)} target(s), D1_extremal, "
          f"weights {a.weights}, knob {knob}.")
    if a.arm == "arm1":
        print(f"⚠ ARM 1 uses ARM 2's decision-J shares. Arm 1's per-target chain "
              f"lengths span 22,129,531..22,136,532 against arm 2's {tot:,} "
              f"(logs/B_loo_chain_survey.json), a 0.03164% spread, so the share a "
              f"chromosome takes differs by that much from its own.")
        print(f"  ⚠ [CORRECTION - 2026-10-25] The effect on kl is **0.00010%** at "
              f"worst over all 44 targets, MEASURED, not the ~0.06% this line used "
              f"to assert from 'kl scales with share^2'. That reasoning dropped the "
              f"1/T_c factor: kl_c goes as share_c^2 / T_c, and a chromosome with a "
              f"shorter chain has BOTH a smaller share and a smaller T, so the two "
              f"errors very nearly cancel. The old figure was wrong by ~600x in the "
              f"CONSERVATIVE direction (it overstated the error), so nothing "
              f"previously reported was flattered by it.")
        print(f"  ⚠ The old justification -- 'a per-chromosome worker cannot see "
              f"them' -- was also FALSE: all 22 x 44 leave-one-out chain lengths are "
              f"on disk in logs/B_loo_chain_survey.json under "
              f"per_chrom[chrom][target]['t_chain'] and have been since 2026-10-12. "
              f"The approximation is kept because it is accurate to 1e-6 relative, "
              f"NOT because the data is unavailable.")
    print()
    # ⚠⚠ CHECKPOINT PER CHROMOSOME. Each chromosome costs 1.5-3.5 hours of the ~3.4h
    # wall time; holding 22 of them only in memory means ANY later fault -- a reporting
    # bug, an OOM, a stray Ctrl-C -- throws all of it away, which is exactly what
    # happened on 2026-10-25. Appending each result as it lands makes the run resumable
    # with `--resume`, so a re-run pays only for what did not finish.
    ckpt = (a.json.with_name(a.json.stem + ".perchrom.jsonl") if a.json else None)
    per: dict = {}
    if ckpt and a.resume and ckpt.exists():
        for ln in ckpt.read_text().splitlines():
            if not ln.strip():
                continue
            r = json.loads(ln)
            # ⚠ Only reuse a checkpoint written for THIS configuration. A cached
            # chromosome from another arm, weighting or grid is a silent wrong answer,
            # which is worse than recomputing it.
            if (r.get("_arm") == a.arm and r.get("_weights") == a.weights
                    and r.get("_grid") == [knob_key(v) for v in grid]
                    and r.get("_targets") == list(tids)):
                per[r["chrom"]] = r
        print(f"  resumed {len(per)} chromosome(s) from {ckpt.name}: "
              f"{sorted(per)}\n", flush=True)

    todo = [c for c in AUTOSOMES if c not in per]
    t0 = time.time()
    failed: dict = {}
    if todo:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {ex.submit(_one, (c, sh[c], tids, grid if not a.etas else a.taus),
                              {"arm": a.arm, "etas": grid if a.etas else None,
                               "weights": a.weights}): c for c in todo}
            for n, f in enumerate(as_completed(futs), 1):
                r = f.result()
                if "error" in r:
                    # ⚠ flush: the SUCCESS line flushed and the FAILURE line did not,
                    # so a crash could leave the reason invisible in a buffered log.
                    print(f"  [{n}/{len(todo)}] ❌ {r['chrom']}: {r['error']}\n"
                          f"{r['traceback']}", flush=True)
                    failed[r["chrom"]] = r
                    # ⚠ Do NOT raise in here: `SystemExit` inside the `with` block
                    # makes the executor's __exit__ wait out every queued future
                    # first, so the abort takes as long as finishing would have.
                    break
                per[r["chrom"]] = r
                if ckpt:
                    with ckpt.open("a") as fh:
                        fh.write(json.dumps({**r, "_arm": a.arm,
                                             "_weights": a.weights,
                                             "_grid": [knob_key(v) for v in grid],
                                             "_targets": list(tids)}) + "\n")
                print(f"  [{len(per)}/22] ✅ {r['chrom']} "
                      f"[{time.time() - t0:.0f}s]", flush=True)
            if failed:
                ex.shutdown(wait=False, cancel_futures=True)
    if failed:
        raise SystemExit(f"{sorted(failed)} failed; the total would be wrong. "
                         f"{len(per)} of 22 are checkpointed"
                         + (f" in {ckpt.name} -- re-run with --resume." if ckpt else
                            " NOWHERE (no --json given, so nothing was saved)."))

    missing = sorted(set(AUTOSOMES) - set(per))
    assert not missing, f"missing {missing}; a partial sum UNDERSTATES the KL"

    rows = []
    for t in tids:
        for v in grid:
            # ⚠ The per-chromosome parts stay in hand so the sign check can see THEM.
            # Clamping only the 22-term SUM (what this did before) lets one genuinely
            # negative chromosome hide inside a positive total.
            parts = [per[c]["kl"][t][knob_key(v)]["kl"] for c in AUTOSOMES]
            bad_parts = {c: x for c, x in zip(AUTOSOMES, parts) if x < -1e-15}
            kl_raw = math.fsum(parts)
            neg = kl_raw < -1e-15 or bool(bad_parts)
            kl = max(kl_raw, 0.0)
            tv = min(1.0, math.sqrt(kl / 2.0))
            # ⚠⚠ `p_succ_bound` REQUIRES tau in [0, 1) and raises otherwise. In eta
            # mode the knob is an eta, and tanh(eta) is 1.0 to float64 past eta ~ 19,
            # so for every point on a real eta grid the bound is VACUOUS rather than
            # computable. Report that, do not crash on it -- and do not quietly pass
            # the eta in as if it were a tau, which is how `p_succ_bound(96.0)` took
            # down the A1 bound check once already.
            tau_eq = math.tanh(v) if a.etas else v
            psb = 1.0 if tau_eq >= 1.0 else C.p_succ_bound(tau_eq)
            rows.append({"target": t, "knob": knob, "knob_value": v,
                         # ⚠ "tau" is kept because every reader of these JSONs keys on
                         # it (the repo's convention since arm2_attack.py), but in eta
                         # mode it holds an ETA. "knob"/"knob_value" are the honest
                         # fields; read those.
                         "tau": v, "tau_equiv": tau_eq,
                         "kl_exact": kl, "kl_raw": kl_raw, "tv_pinsker": tv,
                         "acc_ceiling": 0.5 * (1.0 + tv),
                         "p_succ_bound": psb, "bound_vacuous": tau_eq >= 1.0,
                         "was_negative": neg,
                         "negative_chroms": sorted(bad_parts)})

    # ⚠⚠ THE JSON IS WRITTEN BEFORE THE REPORT, AND THAT ORDERING IS THE POINT.
    # 2026-10-25: all 22 chromosomes finished, then `_report` raised, and because the
    # write was downstream of it 3.4 hours of compute went in the bin. A report is a
    # formatter; it must never be able to destroy a completed computation. `per_chrom`
    # is included for the same reason -- the expensive part is then recoverable from
    # disk whatever happens next.
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": a.arm, "scope": "genome", "route": "stable",
             "k_states": 86 if a.arm == "arm1" else 88,
             "weights": a.weights, "knob": knob, "knob_values": grid,
             "taus": grid if not a.etas else None,
             "etas": grid if a.etas else None,
             "T": tot, "T_is_arm2_reference": a.arm == "arm1",
             "t_chain": {c: per[c]["t_chain"] for c in AUTOSOMES},
             "shares": sh, "rows": rows,
             "per_chrom": {c: per[c]["kl"] for c in AUTOSOMES}}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")

    _report(rows, grid, len(tids), tot, knob=knob, weights=a.weights, arm=a.arm)
    return 0


def _report(rows, grid, n, tot, *, knob="tau", weights="inv_support",
            arm="arm2") -> None:
    # ⚠ The arm, K, weighting and knob were HARDCODED here as "ARM 2, K = 88" with no
    # weighting at all, so an `--arm arm1 --weights uniform` run printed a table
    # stamped as arm 2 — the reporting contract in CLAUDE.md requires every number to
    # carry its arm, knob, denominator and n, and a wrong stamp makes a figure
    # uninterpretable rather than imprecise.
    k_states = 86 if arm == "arm1" else 88
    print(f"\nGENOME-WIDE privacy ceiling, EXACT. "
          f"**{arm.upper().replace('ARM', 'ARM ')}**, K = {k_states}, "
          f"{'per-target leave-one-out, ' if arm == 'arm1' else ''}22 autosomes, "
          f"T = {tot:,}{' (ARM 2 REFERENCE — arm 1 chains differ per target)' if arm == 'arm1' else ''}, "
          f"shares summing to 1.000000, n = {n} target(s), weights {weights}, "
          f"knob {knob}, D1_extremal, zero Monte-Carlo error.\n")
    hdr = ["metric"] + [f"{knob} = {v:g}" for v in grid]
    print("| " + " | ".join(hdr) + " |"); print("|" + "---|" * len(hdr))
    for lab, k, f in (("`kl_exact` (max over targets)", "kl_exact", ".3e"),
                      ("`tv_pinsker` upper bound", "tv_pinsker", ".3e"),
                      ("`acc_ceiling` *", "acc_ceiling", ".9f"),
                      ("`p_succ_bound` *", "p_succ_bound", ".4f")):
        cells = []
        for v in grid:
            # ⚠ `max()` over an empty generator raises, and that is precisely how the
            # grid mismatch destroyed a finished run. Name the missing cell instead.
            vals = [r[k] for r in rows if r["knob_value"] == v]
            if not vals:
                raise AssertionError(
                    f"no rows at {knob} = {v!r} for metric {k}; the grid passed to "
                    f"_report does not match the grid the rows were built from")
            cells.append(f"{max(vals):{f}}")
        print(f"| {lab} | " + " | ".join(cells) + " |")
    if any(r.get("bound_vacuous") for r in rows):
        print(f"| bound comparison meaningful? | " + " | ".join(
            "**no (vacuous)**" if any(r["knob_value"] == v and r.get("bound_vacuous")
                                      for r in rows) else "yes"
            for v in grid) + " |")
    bad = [r for r in rows if r["acc_ceiling"] > r["p_succ_bound"] + 1e-12]
    negs = [r for r in rows if r["was_negative"]]
    print(f"\n{'❌ GATE FAILED' if bad else '✅ gate passes'}: `acc_ceiling` <= "
          f"`p_succ_bound` in {len(rows) - len(bad)} of {len(rows)} cells.")
    if negs:
        print(f"⚠ {len(negs)} cell(s) came back marginally NEGATIVE and were clamped to "
              f"0 — at those `tau` the KL is still at the noise floor even for this "
              f"route, and the figure is an upper bound on a quantity too small to "
              f"resolve, not a measurement.")
    print(f"⚠⚠ THIS IS THE GENOME-WIDE RELEASE — the same object the utility movement is "
          f"measured on — so the two may be read together. The chr21 figures in "
          f"`logs/B_exact_kl.md` are at `beta_share` = 1, 66.6x more generous per "
          f"position, and are NOT this quantity.")
    print(f"⚠ Pinsker is an inequality: the true ceiling is closer to 0.5 than "
          f"`acc_ceiling` prints, and the looseness flatters us.")
    print(f"⚠ This bounds EVERY attacker, including unwritten ones. A measured attacker's "
          f"accuracy is a LOWER bound on distinguishability (constraint 7).")


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
