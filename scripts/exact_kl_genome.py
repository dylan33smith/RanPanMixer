#!/usr/bin/env python
"""The GENOME-WIDE privacy ceiling, exactly. The measurement the paper was missing.

⚠⚠ WHY THIS EXISTS. Every privacy figure on record came from `scripts/exact_kl.py`, which
hardcodes `CHROM = "chr21"` and calls `utility.beta_v` with NO `share`, so `share` takes
its default of 1.0 — the whole budget on one chromosome. That configuration is **66.6x
more generous per position than the genome-wide release**, while every utility figure IS
the genome-wide release. So the paper's two headline results were measured on two
different objects, and the sentence joining them could not be supported. This closes that
gap: one release, its own privacy ceiling beside its own utility movement.

ROUTE. KL is ADDITIVE across chromosomes. The release factorises given the target — one
budget split into per-chromosome shares by decision J, each chromosome sampled
independently — so `KL(Q_p || Q_q) = sum_c KL_c`, with each chromosome's KL computed at
its own `beta_share`. The Pinsker step is taken ONCE on the total, never per chromosome:
`sqrt` is concave, so summing per-chromosome `tv_pinsker` values would be wrong AND would
flatter us.

⚠ The per-chromosome KL is computed by `scripts/exact_kl.py:exact_kl`, imported unchanged
rather than reimplemented, so the genome-wide figure inherits its sign check (KL >= 0 is a
theorem, not a convention) and its diploid half-`eta_tau` handling.

⚠ ARM 2: the target is external, nothing is held out, K = 88. The decoy is D1_extremal,
built exactly as the sampled attack builds it — every called allele flipped.

    python scripts/exact_kl_genome.py --targets 1 --workers 22
"""

from __future__ import annotations

import argparse
import importlib.util
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

_spec = importlib.util.spec_from_file_location("exact_kl", REPO / "scripts" / "exact_kl.py")
EK = importlib.util.module_from_spec(_spec)
sys.modules["exact_kl"] = EK
_spec.loader.exec_module(EK)

PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
AUTOSOMES = loo.AUTOSOMES


def shares() -> dict[str, float]:
    """`beta_share` by variant count, decision J — the same rule the release used."""
    T = {c: int(np.load(REPO / "data" / "cohort" / c / "in_chain.npy").sum())
         for c in AUTOSOMES}
    tot = sum(T.values())
    sh = {c: T[c] / tot for c in AUTOSOMES}
    s = sum(sh.values())
    assert abs(s - 1.0) < 1e-12, f"shares sum to {s!r}, not 1; the calibration is void"
    return sh


def kl_one_chrom(chrom: str, share: float, targets: list[str], taus) -> dict:
    """Every target's per-chromosome KL contribution, at this chromosome's share."""
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=chrom)
    # ⚠⚠ THE SHARE IS THE WHOLE POINT. exact_kl.py omits it and so measures chr21 at the
    # whole budget; passing it is what makes this the genome-wide release.
    util = utility.beta_v(np.load(d / "support.npy"), in_chain=ic, share=share,
                          chrom=chrom)
    st = {"hap": hap, "model": m, "util": util,
          "base": ffbs.baseline_log_norm(m), "ic": ic}
    out = {}
    for t in targets:
        p = np.load(REPO / "data" / "targets" / t / chrom / "path.npy")[ic]
        ext = p.copy()
        called = ext >= 0
        ext[called] = np.where(ext[called] == 0, 1, 0)      # D1_extremal
        out[t] = {str(tau): EK.exact_kl(st, p, ext, tau)["kl"] for tau in taus}
    return {"chrom": chrom, "share": share, "t_chain": int(ic.sum()), "kl": out}


def _one(args):
    try:
        return kl_one_chrom(*args)
    except Exception as e:
        return {"chrom": args[0], "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=1)
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--workers", type=int, default=22)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    sh = shares()
    tids = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                  if p.is_dir() and (p / "chr1" / "path.npy").exists())[: a.targets]
    assert tids, "no prepared targets with a chr1 path"
    print(f"GENOME-WIDE exact KL. ARM 2, K = 88, nothing held out, 22 autosomes, "
          f"shares summing to 1.000000, n = {len(tids)} target(s), D1_extremal decoy.")
    print(f"⚠ KL is additive across chromosomes; the Pinsker step is taken ONCE on the "
          f"total, because sqrt is concave and per-chromosome Pinsker would flatter us.\n")

    t0, per = time.time(), {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_one, (c, sh[c], tids, a.taus)): c for c in AUTOSOMES}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if "error" in r:
                print(f"  [{n}/22] ❌ {r['chrom']}: {r['error']}", flush=True)
                print(r["traceback"])
                raise SystemExit("a chromosome failed; the total would be wrong")
            per[r["chrom"]] = r
            print(f"  [{n}/22] ✅ {r['chrom']} (share {r['share']:.6f}, "
                  f"T = {r['t_chain']:,}) [{time.time() - t0:.0f}s]", flush=True)

    # ⚠ assert completeness: a KL summed over 21 chromosomes is not the release's KL,
    # and it would be SMALLER, which flatters us.
    missing = sorted(set(AUTOSOMES) - set(per))
    assert not missing, f"missing chromosomes {missing}; the total would understate the KL"

    rows = []
    for t in tids:
        for tau in a.taus:
            kl = sum(per[c]["kl"][t][str(tau)] for c in AUTOSOMES)
            tv = min(1.0, math.sqrt(kl / 2.0))
            rows.append({"target": t, "tau": tau, "kl_exact": kl, "tv_pinsker": tv,
                         "acc_ceiling": 0.5 * (1.0 + tv),
                         "p_succ_bound": C.p_succ_bound(tau)})
    _report(rows, a.taus, len(tids))
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm2", "scope": "genome", "k_states": 88, "taus": a.taus,
             "shares": sh, "rows": rows,
             "per_chrom": {c: per[c]["kl"] for c in AUTOSOMES}}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(rows, taus, n) -> None:
    assert rows, "no rows; the report would be silent"
    print(f"\nGENOME-WIDE privacy ceiling, EXACT. **ARM 2**, K = 88, 22 autosomes, "
          f"T = 22,134,562, shares summing to 1.000000, n = {n} target(s), "
          f"D1_extremal, zero Monte-Carlo error.\n")
    hdr = ["metric"] + [f"{t:g}" for t in taus]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for lab, k, f in (("`kl_exact` (max over targets)", "kl_exact", ".3e"),
                      ("`tv_pinsker` upper bound", "tv_pinsker", ".3e"),
                      ("`acc_ceiling` *", "acc_ceiling", ".6f"),
                      ("`p_succ_bound` *", "p_succ_bound", ".4f")):
        cells = [max(r[k] for r in rows if r["tau"] == t) for t in taus]
        print(f"| {lab} | " + " | ".join(f"{c:{f}}" for c in cells) + " |")

    worst = max(rows, key=lambda r: r["acc_ceiling"])
    bad = [r for r in rows if r["acc_ceiling"] > r["p_succ_bound"] + 1e-12]
    print(f"\n{'❌ GATE FAILED' if bad else '✅ gate passes'}: `acc_ceiling` <= "
          f"`p_succ_bound` in {len(rows) - len(bad)} of {len(rows)} cells.")
    print(f"⚠⚠ This is the GENOME-WIDE release — the same object the utility movement is "
          f"measured on — so it may be read beside it. The chr21 figures in "
          f"`logs/B_exact_kl.md` are at `beta_share` = 1, 66.6x more generous per "
          f"position, and are NOT this quantity.")
    print(f"⚠ Pinsker is an inequality: the true ceiling is closer to 0.5 than "
          f"`acc_ceiling` prints, and the looseness flatters us.")
    print(f"⚠ An attacker that fails is not evidence the bound holds; this is an UPPER "
          f"bound on every attacker, which is the stronger statement.")
    print(f"⚠ Worst cell: `acc_ceiling` {worst['acc_ceiling']:.6f} at `tau` = "
          f"{worst['tau']:g} against a bound of {worst['p_succ_bound']:.4f}.")


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
