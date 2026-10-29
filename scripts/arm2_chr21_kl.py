"""ARM 2's chr21 `kl_exact` and `acc_ceiling`, PER TARGET, at `beta_share` = 1.

⚠⚠ WHY THIS EXISTS. Arm 2's reported grid -- the seven `eta` [0.0, 75.7, 151.4, 227.1,
302.8, 378.5, 454.19] that label `logs/B_acc_arm2/a1_acc.md` and `a34_acc.md` 0.50 to
0.80 -- was fitted on **HG00096 ALONE** (`logs/B_eta_grid_v2/chr21_uniform.json`, whose
`target` field says so), and then used for all TEN Arm 2 targets. `scripts/arm2_attack.py`
takes `--acc-ceiling` as ETA=VAL pairs keyed on the eta label only (`:454-456`, `:479`),
so its bound check compares accuracy POOLED OVER 10 TARGETS against ONE target's ceiling.

⚠ ARM 1 MEASURED THAT THIS IS UNSAFE THERE: `c` = kl/eta^2 spans 8.664432e-10 to
1.176051e-09 across its 44 targets, a 35.7331% spread, "so a single shared `eta` could
NOT have placed every target on one ceiling, and the per-target grid is required rather
than preferred" (`docs/memory.md` 2026-10-27, B-EVL-arm1-chr21-grid).

⚠ BUT ARM 1'S SPREAD IS NOT ARM 2'S, AND THE DIFFERENCE CUTS OUR WAY. Arm 1 compounds
leave-one-out: K = 86, the chain rebuilt per target, each target's own `beta_share` from
`loo.loo_shares`. Arm 2 holds nothing out -- K = 88 and share 1 for every target -- so the
only target-dependence left is WHICH POSITIONS the target agrees with. The spread is
therefore expected to be smaller. Expected is not measured, which is why this runs.

⚠ This is a DRIVER over `exact_kl_stable.one_chrom`, with the ceiling arithmetic IMPORTED
from `arm1_chr21_kl` rather than copied. One definition of `acc_ceiling`, not three -- a
second copy is how it drifts, and this file would have been the third.

⚠ `kl_exact = c * eta^2` is the LEADING term, not an identity, so the solved grid is a
STARTING POINT: `--verify` measures `acc_ceiling` AT the solved `eta`, refines once and
re-measures. A grid quoted from a fit is a derived number wearing a measured number's
clothes.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from exact_kl_stable import knob_key, one_chrom                      # noqa: E402
# ⚠ IMPORTED, never reimplemented -- see the module docstring.
from arm1_chr21_kl import CEILINGS, ceiling_from_kl, eta_for_ceiling  # noqa: E402

CHROM = "chr21"
K_STATES = 88
# ⚠ Arm 2's chr21 release is at share 1.0, which `one_chrom` spells as share=None
#   (`exact_kl_stable.py:216-220` prints "`beta_share` = 1, K = 88" for exactly this
#   call). Passing 1.0 is NOT the same thing and would re-scale beta.
SHARE = None
# The grid the Arm 2 attack tables were run at, so the fit provably covers the table it
# will relabel. Read from the run rather than retyped.
ATTACK_RUN = REPO / "logs" / "B_acc_arm2" / "a1_acc.json"


def refine_eta(eta: float, kl: float, label: float) -> float:
    """One multiplicative correction, `kl` being locally ~ `c * eta^2`.

    ⚠ The same two lines as `arm1_chr21_kl.py:124-125`; held to it by
    `tests/test_arm2_chr21_kl.py`, because two copies of a refinement that silently
    disagree would move a reported ceiling without failing anything.
    """
    kl_t = 2.0 * (2.0 * float(label) - 1.0) ** 2
    return eta * math.sqrt(kl_t / kl) if kl > 0 else eta


def one_target(target: str, etas, weights: str) -> dict:
    try:
        out = one_chrom(CHROM, SHARE, [target], None,
                        arm="arm2", etas=list(etas), weights=weights)
        per_target = out["kl"]
        assert target in per_target, (
            f"{target} absent from one_chrom's result (keys {sorted(per_target)})")
        return {"target": target,
                "kl": {knob_key(e): float(per_target[target][knob_key(e)]["kl"])
                       for e in etas}}
    except Exception as e:
        return {"target": target, "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def verify_target(target: str, etas_by_ceiling: dict, weights: str) -> dict:
    try:
        labels = [k for k in etas_by_ceiling if float(k) > 0.5]
        etas = [etas_by_ceiling[k] for k in labels]
        out = one_chrom(CHROM, SHARE, [target], None,
                        arm="arm2", etas=etas, weights=weights)["kl"][target]
        first, refined = {}, {}
        for lab, e in zip(labels, etas):
            kl = float(out[knob_key(e)]["kl"])
            _, acc = ceiling_from_kl(kl)
            first[lab] = {"eta": e, "kl_exact": kl, "acc_ceiling": acc,
                          "error": acc - float(lab)}
            refined[lab] = refine_eta(e, kl, float(lab))
        out2 = one_chrom(CHROM, SHARE, [target], None, arm="arm2",
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


def attack_targets() -> list[str]:
    """The ten targets the Arm 2 attack tables were actually run on."""
    tids = sorted(json.loads(ATTACK_RUN.read_text())["targets"])
    on_disk = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                     if p.is_dir() and (p / CHROM / "path.npy").exists())
    missing = sorted(set(tids) - set(on_disk))
    if missing:
        raise SystemExit(
            f"{missing} were attacked but have no {CHROM} path on disk, so their "
            f"ceiling cannot be fitted and the table cannot be relabelled.")
    return tids


def run_verify(path: Path, workers: int) -> int:
    fit = json.loads(path.read_text())
    per, weights = fit["per_target"], fit["weights"]
    tids = sorted(per)
    print(f"VERIFY: ARM 2, {CHROM} at `beta_share` = 1, K = {K_STATES}. Measuring "
          f"`acc_ceiling` AT the solved `eta` for {len(tids)} targets, then one "
          f"refinement and a re-measure. weights {weights}.\n")
    res, failed = {}, []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_target, t, per[t]["grid"], weights): t for t in tids}
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
    print(f"\n**ARM 2 chr21 GRID, PER TARGET, VERIFIED BY MEASUREMENT.** Metrics are "
          f"ROWS, `acc_ceiling` is the COLUMN. n = {len(tids)} targets, K = "
          f"{K_STATES}, `beta_share` = 1.\n")
    hdr = ["metric"] + labels
    print("| " + " | ".join(hdr) + " |"); print("|" + "---|" * len(hdr))
    for lab, fn, fmt in (("`eta_tau` median [refined]",
                          lambda vs: sorted(vs)[len(vs) // 2], ",.2f"),
                         ("`eta_tau` min [refined]", min, ",.2f"),
                         ("`eta_tau` max [refined]", max, ",.2f")):
        cells = [f"{fn([res[t]['refined'][l]['eta'] for t in tids]):{fmt}}"
                 for l in labels]
        print(f"| {lab} | " + " | ".join(cells) + " |")
    # ⚠⚠ THE ROW THIS RUN EXISTS FOR: how far the single-target label was off for the
    #    nine targets it was never fitted on.
    shared = json.loads(Path(fit["shared_grid_source"]).read_text())
    shared_eta = {f"{r['acc_ceiling']:.2f}": r["eta"] for r in shared["rows"]}
    cells = []
    for l in labels:
        e0 = shared_eta.get(l)
        cells.append("n/a (no shared eta at this label)" if e0 is None else
                     f"{max(abs(ceiling_from_kl(res[t]['first'][l]['kl_exact'] * (e0 / res[t]['first'][l]['eta']) ** 2)[1] - float(l)) for t in tids):.2e}")
    print(f"| worst |`acc_ceiling` - label| AT THE SHARED `eta` * | "
          + " | ".join(cells) + " |")
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
    print(f"⚠ This grid labels ARM 2's chr21 attack at share 1.0. Arm 2's UTILITY "
          f"numbers are GENOME-WIDE and need their own fit -- "
          f"`logs/B_eta_grid_v2/genome_uniform.json` is HG00096-only too.")
    outp = path.with_name(path.stem + "_verified.json")
    outp.write_text(json.dumps({"arm": "arm2", "scope": CHROM, "beta_share": 1.0,
                                "k_states": K_STATES, "weights": weights,
                                "worst_abs_error": worst, "per_target": res},
                               indent=1) + "\n")
    print(f"\nJSON -> {outp}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", type=Path, default=None,
                    help="a fit JSON: measure `acc_ceiling` at its solved `eta`, "
                         "refine once, re-measure. Quotable only after this passes.")
    ap.add_argument("--etas", type=float, nargs="+", default=[227.1, 454.19],
                    help="two non-zero points: one fixes `c`, the second confirms the "
                         "quadratic. Arm 2 at share 1 reaches 0.80 near 454, i.e. "
                         "~61x lower than arm 1's chr21 share of ~0.0150 needs.")
    ap.add_argument("--weights", default="uniform",
                    choices=("uniform", "inv_support"))
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--shared-grid",
                    default=str(REPO / "logs" / "B_eta_grid_v2" / "chr21_uniform.json"),
                    help="the HG00096-only grid whose labels this run is checking.")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.verify:
        return run_verify(a.verify, a.workers)

    if any(e <= 0 for e in a.etas):
        raise SystemExit("every `eta` must be > 0: `c` is undefined at eta = 0, where "
                         "kl is identically 0 and the fit carries no information.")
    tids = attack_targets()

    print(f"ARM 2, {CHROM} at `beta_share` = 1, `kl_exact` PER TARGET. K = "
          f"{K_STATES}, n = {len(tids)} targets, weights {a.weights}, knob "
          f"`eta_tau`, {len(a.etas)} points, {a.workers} workers.")
    print(f"⚠ Nothing is held out in arm 2, so every target shares one panel and one "
          f"share; the only target-dependence is which positions it agrees with.")
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
            print(f"  [{n}/{len(tids)}] ✅ {r['target']}", flush=True)

    # ⚠ A fan-out that filters its own results cannot tell "found nothing" from
    # "ran nothing". Assert completeness and fail loudly.
    if failed:
        raise SystemExit(f"{sorted(failed)} failed; a partial grid would be solved "
                         f"from an incomplete fit. Nothing is written.")
    missing = sorted(set(tids) - set(res))
    assert not missing, f"missing {missing} with no error recorded -- worse than a crash"

    rows, worst_resid = [], 0.0
    for t in tids:
        kls = res[t]["kl"]
        cs = [kls[knob_key(e)] / (e * e) for e in a.etas]
        c = sum(cs) / len(cs)
        resid = (max(cs) / min(cs) - 1.0) if min(cs) > 0 else float("inf")
        worst_resid = max(worst_resid, resid)
        for e in a.etas:
            kl = kls[knob_key(e)]
            tv, acc = ceiling_from_kl(kl)
            rows.append({"target": t, "arm": "arm2", "scope": CHROM,
                         "beta_share": 1.0, "k_states": K_STATES,
                         "knob": "eta_tau", "knob_value": e,
                         "kl_exact": kl, "tv_pinsker": tv, "acc_ceiling": acc})
        res[t]["c"] = c
        res[t]["c_residual"] = resid
        res[t]["grid"] = {f"{acc:.2f}": eta_for_ceiling(acc, c) for acc in CEILINGS}

    print(f"\n⚠ FIT CONDITIONING: spread in `c` = kl/eta^2 across the {len(a.etas)} "
          f"points, worst over all {len(tids)} targets: {worst_resid:.3e} "
          f"({100 * worst_resid:.4f}%). The solved grid is a STARTING POINT and must "
          f"be verified by measurement (--verify), never quoted from the fit alone.")

    # ⚠⚠ THE NUMBER THIS RUN EXISTS FOR, and it is a SPREAD ACROSS TARGETS, not across
    #    the fit's etas. Arm 1's is 35.7331%; if arm 2's is of that order then nine of
    #    its ten columns carry a ceiling that is not theirs.
    cvals = {t: res[t]["c"] for t in tids}
    lo_t, hi_t = min(cvals, key=cvals.get), max(cvals, key=cvals.get)
    lo, hi = cvals[lo_t], cvals[hi_t]
    spread = (hi / lo - 1.0) if lo > 0 else float("inf")
    print(f"\n**`c` = kl/eta^2 ACROSS THE {len(tids)} ARM 2 TARGETS: {lo:.6e} "
          f"({lo_t}) to {hi:.6e} ({hi_t}), a {100 * spread:.4f}% spread.**")
    print(f"  ARM 1's, for comparison: 8.664432e-10 to 1.176051e-09, 35.7331% "
          f"(44 targets, K = 86, leave-one-out).")
    print(f"  ⚠ A spread here is NOT a defect in the mechanism. It is how far the "
          f"HG00096-only label is from the other targets' true ceiling, and the "
          f"direction matters: a target whose `c` is LARGER than HG00096's sits at a "
          f"HIGHER true ceiling at the shared `eta` (its bound is looser than printed, "
          f"so its null is safe), and one whose `c` is SMALLER sits LOWER (its bound is "
          f"TIGHTER than printed, which is the case that could hide an exceedance).")
    c0 = cvals.get("HG00096")
    if c0:
        below = sorted(t for t in tids if cvals[t] < c0)
        print(f"  Targets with `c` BELOW HG00096's ({c0:.6e}) -- the ones whose "
              f"printed ceiling is too generous: {below or 'none'}")

    payload = {"arm": "arm2", "scope": CHROM, "beta_share": 1.0,
               "k_states": K_STATES, "weights": a.weights, "etas": list(a.etas),
               "ceilings": list(CEILINGS), "shared_grid_source": a.shared_grid,
               "c_spread_across_targets": spread,
               "c_min": {"target": lo_t, "c": lo}, "c_max": {"target": hi_t, "c": hi},
               "worst_c_residual_within_target": worst_resid,
               "rows": rows,
               "per_target": {t: {k: v for k, v in res[t].items() if k != "kl"}
                              | {"kl": res[t]["kl"]} for t in tids}}
    outp = a.json or (REPO / "logs" / "B_acc_arm2" / "arm2_chr21_kl.json")
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"\nJSON -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
