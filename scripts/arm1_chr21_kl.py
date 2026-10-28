"""ARM 1's chr21 `kl_exact` and `acc_ceiling`, at each target's OWN leave-one-out share.

⚠⚠ WHY THIS EXISTS. `scripts/exact_kl_stable.py`'s CLI is GENOME-WIDE ONLY: it builds
the decision-J shares over all 22 autosomes and fans out one worker per CHROMOSOME.
Arm 1 is chr21 ONLY from 2026-10-27 (Dylan; arm 2 carries genome-wide alone), so what
is needed is ONE chromosome at that target's own `beta_share` for all 44 targets --
which `one_chrom()` already computes. This is a DRIVER over that function and adds no
arithmetic of its own: `kl` comes from `kl_stable`, and `tv_pinsker` / `acc_ceiling`
are the same `min(1, sqrt(kl/2))` and `0.5 * (1 + tv)` as `exact_kl_stable.py:399,415`.
One definition, not two -- a second copy is how `acc_ceiling` would drift.

⚠ THE SHARE IS PER TARGET AND IT MOVES. Leave-one-out changes every `T_c`, so each of
the 44 targets has its own 22 shares (`loo.loo_shares`). chr21's share is taken from
THAT target's shares, never from arm 2's. Using arm 2's share would mean `sum_v beta_v`
is not the scope the bound covers.

⚠ FANS OUT OVER TARGETS, not chromosomes. `exact_kl_stable` parallelises per
chromosome, which on a one-chromosome run would pin everything to a single core.

⚠ `kl_exact = c * eta^2` EXACTLY, so one non-zero `eta` determines `c` and the second
only CONFIRMS the quadratic. Both are computed and the residual is reported, because a
fit that cannot fail is decoration.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from cohort import loo                                        # noqa: E402
from exact_kl_stable import knob_key, one_chrom                # noqa: E402

CHROM = "chr21"
SURVEY = REPO / "logs" / "B_loo_chain_survey.json"
# The released regime, fixed by Dylan 2026-10-27: 1.0 is essentially a guaranteed
# attack, so the grid stops at 0.80.
CEILINGS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80)


def loo_share_for(target: str, survey: dict) -> float:
    """chr21's `beta_share` for ONE held-out individual, from the measured chains."""
    per = survey["per_chrom"]
    t_by_chrom = {}
    for c in loo.AUTOSOMES:
        row = per[c][target]
        # ⚠ REAL targets carry 't_chain'; only the '__arm2__' pseudo-row carries 't'.
        # Reading the field name off that pseudo-row gives KeyError on every real
        # target, which is how this was first written. Assert rather than .get() --
        # a missing chain length must not become a silently wrong share.
        assert "t_chain" in row, (
            f"{c}/{target} has no 't_chain' (keys {sorted(row)}); a share computed "
            f"from the wrong field would be silently wrong, not absent")
        t_by_chrom[c] = int(row["t_chain"])
    return loo.loo_shares(t_by_chrom)[CHROM]


def one_target(target: str, etas, weights: str) -> dict:
    try:
        share = loo_share_for(target, json.loads(SURVEY.read_text()))
        out = one_chrom(CHROM, share, [target], None,
                        arm="arm1", etas=list(etas), weights=weights)
        # ⚠ `one_chrom` returns {"chrom":..., "kl": {target: {knob_key: {...}}}, ...}
        # -- the per-target map is UNDER "kl", as `exact_kl_stable.py:395` reads it.
        # Indexing out[target] directly raises KeyError on the target's own name.
        per_target = out["kl"]
        assert target in per_target, (
            f"{target} absent from one_chrom's result (keys {sorted(per_target)})")
        return {"target": target, "share": share,
                "kl": {knob_key(e): float(per_target[target][knob_key(e)]["kl"])
                       for e in etas}}
    except Exception as e:
        return {"target": target, "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def ceiling_from_kl(kl: float) -> tuple[float, float]:
    kl = max(float(kl), 0.0)
    tv = min(1.0, math.sqrt(kl / 2.0))
    return tv, 0.5 * (1.0 + tv)


def eta_for_ceiling(acc: float, c: float) -> float:
    """The `eta` realising `acc` on THIS target, from `kl = c * eta^2`."""
    tv = 2.0 * acc - 1.0
    return math.sqrt(2.0 * tv * tv / c) if c > 0 else float("nan")


def verify_target(target: str, etas_by_ceiling: dict, weights: str,
                  share: float) -> dict:
    """MEASURE `acc_ceiling` at the solved `eta`, then refine once and measure again.

    ⚠ This is the standard both arm 2 grids were held to: all seven of their points
    land within 6e-06 (genome) and 1.2e-05 (chr21) of target BY MEASUREMENT. A grid
    quoted from a fit is a derived number wearing a measured number's clothes.
    """
    try:
        labels = [k for k in etas_by_ceiling if float(k) > 0.5]
        etas = [etas_by_ceiling[k] for k in labels]
        out = one_chrom(CHROM, share, [target], None,
                        arm="arm1", etas=etas, weights=weights)["kl"][target]
        first = {}
        refined = {}
        for lab, e in zip(labels, etas):
            kl = float(out[knob_key(e)]["kl"])
            _, acc = ceiling_from_kl(kl)
            first[lab] = {"eta": e, "kl_exact": kl, "acc_ceiling": acc,
                          "error": acc - float(lab)}
            # one multiplicative correction: kl is locally ~ c * eta^2
            kl_t = 2.0 * (2.0 * float(lab) - 1.0) ** 2
            refined[lab] = e * math.sqrt(kl_t / kl) if kl > 0 else e
        out2 = one_chrom(CHROM, share, [target], None, arm="arm1",
                         etas=[refined[l] for l in labels], weights=weights)["kl"][target]
        second = {}
        for lab in labels:
            e = refined[lab]
            kl = float(out2[knob_key(e)]["kl"])
            _, acc = ceiling_from_kl(kl)
            second[lab] = {"eta": e, "kl_exact": kl, "acc_ceiling": acc,
                           "error": acc - float(lab)}
        return {"target": target, "first": first, "refined": second}
    except Exception as e:
        return {"target": target, "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def run_verify(path: Path, workers: int) -> int:
    fit = json.loads(path.read_text())
    per = fit["per_target"]
    weights = fit["weights"]
    tids = sorted(per)
    print(f"VERIFY: ARM 1, {CHROM}, measuring `acc_ceiling` AT the solved `eta` for "
          f"{len(tids)} targets, then one refinement and a re-measure. weights "
          f"{weights}, K = 86.\n")
    res, failed = {}, []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_target, t, per[t]["grid"], weights,
                          per[t]["share"]): t for t in tids}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if "error" in r and "first" not in r:
                failed.append(r["target"])
                print(f"  [{n}/{len(tids)}] ❌ {r['target']}: {r['error']}", flush=True)
                continue
            res[r["target"]] = r
            w = max(abs(v["error"]) for v in r["refined"].values())
            print(f"  [{n}/{len(tids)}] ✅ {r['target']} worst |error| after "
                  f"refinement {w:.2e}", flush=True)
    if failed:
        raise SystemExit(f"{sorted(failed)} failed; the grid is NOT verified.")
    assert set(res) == set(tids), f"missing {sorted(set(tids) - set(res))}"

    labels = sorted(next(iter(res.values()))["refined"], key=float)
    print(f"\n**ARM 1 chr21 GRID, VERIFIED BY MEASUREMENT.** Metrics are ROWS, "
          f"`acc_ceiling` is the COLUMN. n = {len(tids)} targets.\n")
    hdr = ["metric"] + labels
    print("| " + " | ".join(hdr) + " |"); print("|" + "---|" * len(hdr))
    for lab, key, fn, fmt in (
            ("`eta_tau` median [refined]", "refined",
             lambda vs: sorted(vs)[len(vs) // 2], ",.1f"),
            ("`eta_tau` min [refined]", "refined", min, ",.1f"),
            ("`eta_tau` max [refined]", "refined", max, ",.1f")):
        cells = [f"{fn([res[t][key][l]['eta'] for t in tids]):{fmt}}" for l in labels]
        print(f"| {lab} | " + " | ".join(cells) + " |")
    for lab, key in (("worst |error| FROM THE FIT", "first"),
                     ("worst |error| AFTER REFINEMENT *", "refined")):
        cells = [f"{max(abs(res[t][key][l]['error']) for t in tids):.2e}"
                 for l in labels]
        print(f"| {lab} | " + " | ".join(cells) + " |")
    worst = max(abs(res[t]["refined"][l]["error"]) for t in tids for l in labels)
    print(f"\n{'✅' if worst < 1e-4 else '❌'} worst |`acc_ceiling` - target| over "
          f"{len(tids)} targets x {len(labels)} points after one refinement: "
          f"{worst:.3e}")
    print(f"⚠ Pinsker is an inequality: the TRUE ceiling is closer to 0.5 than every "
          f"figure here.")
    outp = path.with_name(path.stem + "_verified.json")
    outp.write_text(json.dumps({"arm": "arm1", "scope": CHROM, "k_states": 86,
                                "weights": weights, "worst_abs_error": worst,
                                "per_target": res}, indent=1) + "\n")
    print(f"\nJSON -> {outp}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", type=Path, default=None,
                    help="a fit JSON: measure `acc_ceiling` at its solved `eta`, "
                         "refine once, re-measure. The grid is quotable only after "
                         "this passes.")
    ap.add_argument("--targets", type=int, default=44)
    ap.add_argument("--etas", type=float, nargs="+", default=[15000.0, 30000.0],
                    help="two non-zero points: one fixes `c`, the second confirms "
                         "the quadratic. chr21 at share ~0.0150 needs eta ~66.6x the "
                         "share-1 values, so 0.80 sits near 30,000.")
    ap.add_argument("--weights", default="uniform",
                    choices=("uniform", "inv_support"))
    ap.add_argument("--workers", type=int, default=22)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.verify:
        return run_verify(a.verify, a.workers)

    survey = json.loads(SURVEY.read_text())
    members = [m for m in sorted(survey["per_chrom"][CHROM]) if m != "__arm2__"]
    tids = members[: a.targets]
    assert tids, "no targets; the survey's per_chrom is empty for " + CHROM
    if any(e <= 0 for e in a.etas):
        raise SystemExit("every `eta` must be > 0: `c` is undefined at eta = 0, where "
                         "kl is identically 0 and the fit carries no information.")

    print(f"ARM 1, {CHROM} ONLY, `kl_exact` at each target's OWN leave-one-out "
          f"`beta_share`. K = 86, n = {len(tids)} targets, weights {a.weights}, "
          f"knob `eta_tau`, {len(a.etas)} points, {a.workers} workers.")
    print(f"⚠ ARM 1 IS chr21 ONLY from 2026-10-27 -- arm 2 carries genome-wide alone.")
    print(f"⚠ `acc_ceiling` is a PINSKER bound: the true ceiling is closer to 0.5 "
          f"than every figure below.\n")

    res: dict[str, dict] = {}
    failed: list[str] = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one_target, t, a.etas, a.weights): t for t in tids}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if "error" in r:
                failed.append(r["target"])
                print(f"  [{n}/{len(tids)}] ❌ {r['target']}: {r['error']}\n"
                      f"{r['traceback']}", flush=True)
                continue
            res[r["target"]] = r
            print(f"  [{n}/{len(tids)}] ✅ {r['target']} share {r['share']:.10f}",
                  flush=True)

    # ⚠ A fan-out that filters its own results cannot tell "found nothing" from
    # "ran nothing". Assert completeness and fail loudly.
    if failed:
        raise SystemExit(f"{sorted(failed)} failed; a partial grid would be solved "
                         f"from an incomplete fit. Nothing is written.")
    missing = sorted(set(tids) - set(res))
    assert not missing, f"missing {missing} with no error recorded -- worse than a crash"

    rows = []
    worst_resid = 0.0
    for t in tids:
        kls = res[t]["kl"]
        cs = [kls[knob_key(e)] / (e * e) for e in a.etas]
        c = sum(cs) / len(cs)
        resid = (max(cs) / min(cs) - 1.0) if min(cs) > 0 else float("inf")
        worst_resid = max(worst_resid, resid)
        for e in a.etas:
            kl = kls[knob_key(e)]
            tv, acc = ceiling_from_kl(kl)
            rows.append({"target": t, "arm": "arm1", "scope": CHROM,
                         "beta_share": res[t]["share"], "k_states": 86,
                         "knob": "eta_tau", "knob_value": e,
                         "kl_exact": kl, "tv_pinsker": tv, "acc_ceiling": acc})
        res[t]["c"] = c
        res[t]["c_residual"] = resid
        res[t]["grid"] = {f"{acc:.2f}": eta_for_ceiling(acc, c) for acc in CEILINGS}

    # ⚠ `kl = c * eta^2` is the LEADING term, not an identity: the genome grid's own
    # `c` drifts 0.0068% monotonically across its range, and this one drifts more
    # because the tilt is pushed harder. So the spread below is a CONDITIONING figure
    # for the fit, NOT a pass/fail gate -- the grid is established by MEASURING
    # `acc_ceiling` at the solved `eta` (--verify), exactly as both arm 2 grids were.
    print(f"\n⚠ FIT CONDITIONING: spread in `c` = kl/eta^2 across the "
          f"{len(a.etas)} points, worst over all {len(tids)} targets: "
          f"{worst_resid:.3e} ({100 * worst_resid:.4f}%). That propagates to about "
          f"{100 * worst_resid / 2:.4f}% in `eta`, so the solved grid below is a "
          f"STARTING POINT and must be verified by measurement, never quoted from "
          f"the fit alone.")

    cvals = [res[t]["c"] for t in tids]
    print(f"\n`c` across targets: min {min(cvals):.6e}  max {max(cvals):.6e}  "
          f"spread {100 * (max(cvals) / min(cvals) - 1):.4f}%")

    print(f"\n**THE ARM 1 chr21 GRID — one `eta_tau` PER TARGET, so every target "
          f"sits EXACTLY on the ceiling.** Metrics are ROWS, `acc_ceiling` is the "
          f"COLUMN.\n")
    hdr = ["`acc_ceiling`"] + [f"{a_:.2f}" for a_ in CEILINGS]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for lab, fn in (("`eta_tau` median", lambda vs: sorted(vs)[len(vs) // 2]),
                    ("`eta_tau` min", min), ("`eta_tau` max", max)):
        cells = [f"{fn([res[t]['grid'][f'{a_:.2f}'] for t in tids]):,.1f}"
                 for a_ in CEILINGS]
        print(f"| {lab} | " + " | ".join(cells) + " |")
    print(f"| `kl_exact` at the ceiling | "
          + " | ".join(f"{2 * (2 * a_ - 1) ** 2:.6f}" for a_ in CEILINGS) + " |")

    print(f"\n⚠ Each target's `eta_tau` differs because its leave-one-out chain and "
          f"share differ. That is REPORTABLE, not a defect: the comparison axis is "
          f"`acc_ceiling`, and matching on it is what standing constraint 7 requires.")

    if a.json:
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(
            {"arm": "arm1", "scope": CHROM, "k_states": 86, "weights": a.weights,
             "knob": "eta_tau", "etas": list(a.etas), "ceilings": list(CEILINGS),
             "n_targets": len(tids), "worst_c_residual": worst_resid,
             "per_target": res, "rows": rows}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
