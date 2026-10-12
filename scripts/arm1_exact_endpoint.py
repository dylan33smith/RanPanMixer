#!/usr/bin/env python
"""ARM 1's PRIMARY endpoint, exactly: `target_fidelity` from the posterior marginals.

⚠⚠ THIS IS THE PRE-REGISTERED PRIMARY (prereg_2026-10-10.md section 3): `E[differ]`, the
expected number of released cells differing from the target, computed EXACTLY from the
marginals — `ffbs.state_marginals` then `ffbs.expected_mismatches`. No sampling, no
differencing. `target_fidelity` is that count over the target's CALLED cells, and it is a
DIFFER RATE: **lower is more faithful** (`docs/terms.md`, direction corrected 2026-10-08).

⚠⚠ THE FINITE-DIFFERENCE ROUTE IS FORBIDDEN HERE and the pre-registration says so: at
genome-wide scope the second derivative of `log Z_p` returns a NEGATIVE `Var_Q(u)`. The
marginals route is the only admissible one, and it needs no step size.

⚠⚠ ARM 1, K = 86, LEAVE-ONE-OUT, AND EVERY NUMBER IS CONDITIONAL. The target is a cohort
member, so the guarantee holds only because both its haplotypes leave the donor panel
(decision E) — and `G` was built from all HPRC assemblies including this target's, so the
guarantee is conditional ON THE GRAPH. Never state an arm 1 number without its arm.

⚠⚠ THE 44 TARGETS CARRY 44 DIFFERENT UTILITY FUNCTIONS, so standing constraint 8 forbids
averaging their utilities. Leave-one-out gives each target its own chain, its own `T`, its
own 22 `beta_share` values and its own `u_ceiling`. This script therefore reports a
DISTRIBUTION across targets (min / median / max) and never a pooled mean, and the spread
is a spread of distinct quantities, NOT a sampling error.

⚠ It recomputes each chain rather than reading the survey, because the
`untraversable_block` set is not monotone in the donor panel and a stale chain would be
silently wrong (`logs/B_loo_chain_survey.md`).

    python scripts/arm1_exact_endpoint.py --workers 24
"""

from __future__ import annotations

import argparse
import json
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


def run_target(target: str, taus) -> dict:
    """Every `tau` for one held-out individual, chromosome-outer."""
    t0 = time.time()
    chains = {c: loo.loo_chain(REPO / "data" / "cohort" / c, target)
              for c in AUTOSOMES}
    shares = loo.loo_shares({c: chains[c].t_chain for c in AUTOSOMES})

    acc = {t: {"differ": 0.0, "e_q_u": 0.0, "called": 0, "ceiling": 0.0}
           for t in taus}
    for c in AUTOSOMES:
        d = REPO / "data" / "cohort" / c
        lc = chains[c]
        ic = lc.in_chain
        hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
        panel = loo.masked_haplotypes(d, target, haplotypes=hap_all)[:, ic]
        rows = loo.sample_index(
            np.load(d / "haplotype_ids.npy", allow_pickle=True))[target]
        tgt = hap_all[list(rows)][:, ic].T                 # (T, 2), its OWN rows
        m = hmm.from_arrays(panel, np.load(d / "genetic_pos.npy")[ic],
                            np.load(d / "positions.npy")[ic],
                            cohort_id=f"{c}_loo_{target}")
        util = utility.beta_v(lc.support, in_chain=ic, share=shares[c], chrom=c)
        # ⚠ BOTH OF THESE ARE `tau`-INDEPENDENT, so they are accumulated ONCE per
        # chromosome and then added to every `tau`. Computing them inside the `tau`
        # loop would multiply each by len(taus) and still print a plausible table.
        called_c = int((tgt[:, 0] >= 0).sum() + (tgt[:, 1] >= 0).sum())
        ceil_c = utility.diploid_utility(
            utility.u_ceiling(tgt[:, 0], util.beta),
            utility.u_ceiling(tgt[:, 1], util.beta))
        for t2 in taus:
            acc[t2]["called"] += called_c
            acc[t2]["ceiling"] += ceil_c

        for tau in taus:
            eta = C.eta_tau(tau)
            for col in (0, 1):
                t1 = tgt[:, col]
                # ⚠ at eta = 0 the potentials are all-ones, so this IS the baseline
                # posterior: the tau = 0 floor needs no special case.
                psi = utility.tilt_potentials(m.allowed, panel, t1, util.beta,
                                              eta / 2.0)
                gam = ffbs.state_marginals(m, ffbs.forward(m, psi), psi)
                acc[tau]["differ"] += ffbs.expected_mismatches(m, panel, t1, gam)
                known = t1 >= 0
                match = (panel == t1[None, :]) & known[None, :]
                # ⚠ 0.5 per haplotype: `diploid_utility` is the MEAN, matching the
                # density exp(eta * u_mean) the sampler draws from.
                acc[tau]["e_q_u"] += 0.5 * float(
                    (gam * match.T).sum(axis=1) @ util.beta)
        del panel, hap_all, m, util

    out = {"target": target, "t_total": sum(chains[c].t_chain for c in AUTOSOMES),
           "k_states": chains[AUTOSOMES[0]].k_states, "seconds": time.time() - t0,
           "cells": {}}
    for tau in taus:
        a = acc[tau]
        out["cells"][str(tau)] = {
            "e_differ": a["differ"], "called_cells": a["called"],
            "target_fidelity": a["differ"] / a["called"] if a["called"] else float("nan"),
            "e_q_u": a["e_q_u"], "u_ceiling": a["ceiling"],
            "e_q_u_over_ceiling": a["e_q_u"] / a["ceiling"] if a["ceiling"] else float("nan"),
        }
    return out


def _one(args):
    target, taus = args
    try:
        return run_target(target, taus)
    except Exception as e:
        return {"target": target, "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1500:]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--from-json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.from_json:
        d = json.loads(a.from_json.read_text())
        print(f"re-reporting {len(d['results'])} saved target(s) from {a.from_json}. "
              f"NO marginals computed.")
        _report(d["results"], [float(t) for t in d["taus"]])
        return 0

    ids = np.load(REPO / "data" / "cohort" / "chr1" / "haplotype_ids.npy",
                  allow_pickle=True)
    members = sorted(loo.sample_index(ids))
    targets = a.targets or members
    print(f"ARM 1 primary endpoint, EXACT. {len(targets)} held-out individuals x "
          f"{len(a.taus)} `tau`, K = 86, 22 autosomes, marginals route, zero "
          f"Monte-Carlo error.\n")

    t0, out = time.time(), []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_one, (t, a.taus)): t for t in targets}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            out.append(r)
            if "error" in r:
                print(f"  [{n}/{len(targets)}] ❌ {r['target']}: {r['error']}", flush=True)
            else:
                print(f"  [{n}/{len(targets)}] ✅ {r['target']}: T = {r['t_total']:,} "
                      f"({r['seconds'] / 60:.1f} min) "
                      f"[{(time.time() - t0) / 60:.0f} min elapsed]", flush=True)

    errs = [r for r in out if "error" in r]
    ok = [r for r in out if "error" not in r]
    # ⚠ assert completeness: a fan-out that filters cannot tell "found nothing" from
    # "ran nothing".
    assert len(ok) + len(errs) == len(targets), "results lost between submit and collect"
    if ok:
        _report(ok, a.taus)
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm1", "taus": a.taus, "k_states": 86, "results": out},
            indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    if errs:
        raise SystemExit(f"INCOMPLETE: {len(errs)} target(s) failed; re-run to finish.")
    return 0


def _report(rows, taus) -> None:
    assert rows, "no targets measured; the report would be silent"
    n = len(rows)
    print(f"\nARM 1 primary endpoint, EXACT. **ARM 1**, K = 86, leave-one-out, diploid, "
          f"22 autosomes, n = {n} held-out individual(s), zero Monte-Carlo error.\n")

    def cells(key):
        return [[r["cells"][str(t)][key] for r in rows] for t in taus]

    hdr = ["metric"] + [f"`tau` = {t:g}" for t in taus]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))

    def line(label, key, fmt, agg="median"):
        vals = cells(key)
        f = {"median": np.median, "min": np.min, "max": np.max}[agg]
        print(f"| {label} | " + " | ".join(f"{f(v):{fmt}}" for v in vals) + " |")

    # primary -> gates -> structure -> context, per the reporting contract
    line("`target_fidelity` * [median of 44]", "target_fidelity", ".6f")
    line("`target_fidelity` [min]", "target_fidelity", ".6f", "min")
    line("`target_fidelity` [max]", "target_fidelity", ".6f", "max")
    line("`e_q_u` [median of 44]", "e_q_u", ".9f")
    line("`u_ceiling` [median of 44]", "u_ceiling", ".9f")
    line("`e_q_u_over_ceiling` [median of 44]", "e_q_u_over_ceiling", ".9f")
    print(f"| `k_states` * | " + " | ".join("86" for _ in taus) + " |")

    lo, hi = str(min(taus)), str(max(taus))
    fid_lo = np.median([r["cells"][lo]["target_fidelity"] for r in rows])
    fid_hi = np.median([r["cells"][hi]["target_fidelity"] for r in rows])
    eu_lo = np.median([r["cells"][lo]["e_q_u"] for r in rows])
    eu_hi = np.median([r["cells"][hi]["e_q_u"] for r in rows])
    print(f"\n⚠⚠ **THE MOVEMENT IS THE FIGURE, NOT THE LEVEL.** Median "
          f"`target_fidelity` goes {fid_lo:.6f} -> {fid_hi:.6f} across `tau` = "
          f"{min(taus):g} to {max(taus):g}, a change of {fid_hi - fid_lo:+.3e}; median "
          f"`e_q_u` goes {eu_lo:.9f} -> {eu_hi:.9f}, a change of {eu_hi - eu_lo:+.3e}. "
          f"The LEVEL at `tau` = 0 is reached by an untilted draw carrying ZERO target "
          f"information, so it is not attributable to the mechanism.")
    print(f"⚠⚠ **THESE ARE {n} DIFFERENT UTILITY FUNCTIONS, NOT {n} SAMPLES OF ONE.** "
          f"Leave-one-out gives each target its own chain, `T`, 22 `beta_share` values "
          f"and `u_ceiling`, so standing constraint 8 forbids averaging their utilities. "
          f"The min/median/max rows are a spread of DISTINCT quantities and the spread is "
          f"NOT a sampling error. Pool measured attack outcomes, never utilities.")
    print(f"⚠⚠ **THE JOINT BOUND DOES NOT HOLD.** These {n} releases have donor panels "
          f"that are each other's TRUE haplotypes. The per-release bound holds; the joint "
          f"bound over the set does not (standing constraint 5).")
    print(f"⚠ `target_fidelity` is a DIFFER rate: LOWER is more faithful.")
    print(f"⚠ ARM 1's `u_ceiling` is NOT comparable to arm 2's — a cohort member's "
          f"alleles are largely carried by the other 86 haplotypes and an external "
          f"target's are not. Compare within an arm, never across.")
    print(f"⚠ The `tau` grid is a declared EVALUATION DEVICE, never a release procedure.")


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
