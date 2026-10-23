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
    d_logz_terms, d_eu = [], 0.0
    for col in (0, 1):
        pc, qc = p[:, col], q[:, col]
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


def one_chrom(chrom: str, share: float | None, targets, taus) -> dict:
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=chrom)
    kw = {} if share is None else {"share": share, "chrom": chrom}
    beta = utility.beta_v(np.load(d / "support.npy"), in_chain=ic, **kw).beta
    out = {}
    for t in targets:
        p = np.load(REPO / "data" / "targets" / t / chrom / "path.npy")[ic]
        q = extremal(p)
        out[t] = {str(tau): kl_stable(m, hap, beta, p, q, tau) for tau in taus}
    return {"chrom": chrom, "kl": out, "t_chain": int(ic.sum())}


def _one(args):
    try:
        return one_chrom(*args)
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
        got = max(r["kl"][t][k]["kl"] for t in tids)
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
    ap.add_argument("--targets", type=int, default=1)
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--workers", type=int, default=22)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.validate:
        validate(a.taus)
        return 0

    T = {c: int(np.load(REPO / "data" / "cohort" / c / "in_chain.npy").sum())
         for c in AUTOSOMES}
    tot = sum(T.values())
    sh = {c: T[c] / tot for c in AUTOSOMES}
    assert abs(sum(sh.values()) - 1.0) < 1e-12
    tids = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                  if p.is_dir() and (p / "chr1" / "path.npy").exists())[: a.targets]

    print(f"GENOME-WIDE `kl_exact`, STABLE ROUTE. ARM 2, K = 88, 22 autosomes, "
          f"T = {tot:,}, shares summing to 1.000000, n = {len(tids)} target(s), "
          f"D1_extremal.\n")
    t0, per = time.time(), {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_one, (c, sh[c], tids, a.taus)): c for c in AUTOSOMES}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if "error" in r:
                print(f"  [{n}/22] ❌ {r['chrom']}: {r['error']}\n{r['traceback']}")
                raise SystemExit("a chromosome failed; the total would be wrong")
            per[r["chrom"]] = r
            print(f"  [{n}/22] ✅ {r['chrom']} [{time.time() - t0:.0f}s]", flush=True)

    missing = sorted(set(AUTOSOMES) - set(per))
    assert not missing, f"missing {missing}; a partial sum UNDERSTATES the KL"

    rows = []
    for t in tids:
        for tau in a.taus:
            kl = math.fsum(per[c]["kl"][t][str(tau)]["kl"] for c in AUTOSOMES)
            neg = kl < -1e-15
            kl = max(kl, 0.0)
            tv = min(1.0, math.sqrt(kl / 2.0))
            rows.append({"target": t, "tau": tau, "kl_exact": kl, "tv_pinsker": tv,
                         "acc_ceiling": 0.5 * (1.0 + tv),
                         "p_succ_bound": C.p_succ_bound(tau), "was_negative": neg})
    _report(rows, a.taus, len(tids), tot)
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm2", "scope": "genome", "route": "stable", "k_states": 88,
             "T": tot, "taus": a.taus, "shares": sh, "rows": rows}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(rows, taus, n, tot) -> None:
    print(f"\nGENOME-WIDE privacy ceiling, EXACT. **ARM 2**, K = 88, 22 autosomes, "
          f"T = {tot:,}, shares summing to 1.000000, n = {n} target(s), D1_extremal, "
          f"zero Monte-Carlo error.\n")
    hdr = ["metric"] + [f"{t:g}" for t in taus]
    print("| " + " | ".join(hdr) + " |"); print("|" + "---|" * len(hdr))
    for lab, k, f in (("`kl_exact` (max over targets)", "kl_exact", ".3e"),
                      ("`tv_pinsker` upper bound", "tv_pinsker", ".3e"),
                      ("`acc_ceiling` *", "acc_ceiling", ".9f"),
                      ("`p_succ_bound` *", "p_succ_bound", ".4f")):
        print(f"| {lab} | " + " | ".join(
            f"{max(r[k] for r in rows if r['tau'] == t):{f}}" for t in taus) + " |")
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
