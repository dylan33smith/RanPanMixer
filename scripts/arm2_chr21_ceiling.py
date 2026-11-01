#!/usr/bin/env python
"""ARM 2's exact `kl_exact` on chr21 for all FOUR decoy conditions.

⚠⚠ WHAT THIS EXISTS TO TEST. Every Arm 2 chr21 table labels its columns with the
`acc_ceiling` of the D1_extremal pair alone (`logs/B_acc_arm2/a1_ceil_n1000_d4.log`,
`a34_ceil_n1000_d4.log`), and that label is a valid bound for the three REAL decoys
only because D1 is the maximum-divergence decoy. Arm 1 asserts the property per cell
(`scripts/arm1_attack.py:342-351`); `scripts/arm2_genome_ceiling.py` MEASURED it at
genome scope on 2026-10-31. ⚠ chr21 was never measured, and chr21 is the scope the
head-to-head and every Arm 2 attack table actually report. This closes it.

⚠⚠ WHAT IT DOES AND DOES NOT CLOSE. `kl_stable` gives KL(Q_p || Q_q) exactly, and
`acc_ceiling` is Pinsker's 0.5*(1 + sqrt(kl/2)) -- an UPPER BOUND on TV, not TV. So
this settles the envelope the COLUMN LABEL rests on: whether D1's exact divergence
dominates every real decoy's at the same `eta`. It does NOT re-rank the conditions by
distinguishability -- that is a TV question, answered on the sampled releases by
E_p[(1 - e^-L)_+], where D4 beats D2 in 6 of 6 tilted chr21 columns. Two different
quantities; do not quote one as the other.

⚠ THE `eta` GRID IS READ, NEVER RE-SOLVED. It comes from
`logs/B_acc_arm2/arm2_chr21_kl_verified.json`'s `refined` block -- the grid whose
`acc_ceiling` was measured at the solved `eta` and refined once, and the same file
`scripts/arm2_attack.py` reads its bound from. A grid re-solved here could not be
compared to the attack cells it is meant to bound.

⚠ THE DECOYS ARE THE ATTACK'S OWN. `load_state`, `concordance`, `DECOY_SEED` and
`pick_decoys` are IMPORTED from `scripts/arm2_attack.py` and the sha256 rng is built
the same way, so every pair here is the pair the sampled cells were scored against.
A local copy of any of them would let this measurement and the attack disagree about
which genome D3 is, and neither would fail.

⚠ D1 IS RE-MEASURED AT EVERY LABEL AS AN INSTRUMENT CHECK, not for its value: the
verified grid already records `kl_exact` at each refined `eta` by construction, so a
run that cannot reproduce it to ~1e-7 is wired wrong -- wrong `beta` scope, wrong
weighting, wrong chain mask -- and would still return plausible numbers for the three
real conditions. That failure is made loud rather than inferred.

    python scripts/arm2_chr21_ceiling.py --workers 30
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src"), str(REPO / "scripts")]

import numpy as np                                              # noqa: E402

from exact_kl_stable import extremal, kl_stable                 # noqa: E402
from arm1_chr21_kl import ceiling_from_kl                       # noqa: E402
from arm2_attack import (CHROM, DECOY_SEED, WEIGHTS, CONDS,     # noqa: E402
                         load_state, pick_decoys)

REAL = tuple(c for c in CONDS if c != "D1_extremal")
VERIFIED = REPO / "logs" / "B_acc_arm2" / "arm2_chr21_kl_verified.json"
# The instrument tolerance. The verified grid's own worst refinement error is
# 3.03e-08 in `acc_ceiling`; this bounds the KL we must reproduce for the SAME eta,
# which is a recomputation of one number by one code path and so should be exact to
# float accumulation order alone.
D1_ATOL = 1e-7


def load_grid() -> dict:
    """Each target's chr21 (label -> (eta, D1 kl_exact)) from the VERIFIED grid."""
    d = json.loads(VERIFIED.read_text())
    if d.get("arm") != "arm2" or d.get("scope") != CHROM:
        raise SystemExit(f"{VERIFIED} is {d.get('arm')}/{d.get('scope')}, not arm2/{CHROM}")
    if d.get("weights") != WEIGHTS:
        raise SystemExit(
            f"{VERIFIED} was solved under weights {d.get('weights')!r} but this run "
            f"uses {WEIGHTS!r}. Two weightings' eta are different quantities wearing "
            f"one name (CLAUDE.md constraint 5) -- do not mix them.")
    if float(d.get("beta_share", -1)) != 1.0:
        raise SystemExit(
            f"{VERIFIED} carries beta_share {d.get('beta_share')!r}, not 1.0, so it is "
            f"not the chr21-scope grid `load_state` builds (scope 'chain', share 1).")
    out = {}
    for t, rec in d["per_target"].items():
        if "refined" not in rec:
            raise SystemExit(f"{t} has no `refined` block in {VERIFIED}: it is a FIT.")
        out[t] = {lab: (float(r["eta"]), float(r["kl_exact"]))
                  for lab, r in rec["refined"].items()}
    return out


def decoys_for(target_id: str, st: dict) -> dict:
    """The attack's own decoy identities for one target. Same rng construction."""
    stable = int.from_bytes(hashlib.sha256(target_id.encode()).digest()[:4], "big")
    rng = np.random.default_rng([DECOY_SEED, stable])
    return {c: d["id"] for c, d in pick_decoys(target_id, st, rng).items()}


def one_cell(job: tuple) -> dict:
    """ONE (target, label): the exact KL of all four conditions at that `eta`."""
    target_id, label, eta, kl_d1, ids = job
    try:
        st = load_state()
        p = st["paths"][target_id]
        qs = {"D1_extremal": extremal(p)}
        for cond in REAL:
            qs[cond] = st["paths"][ids[cond]]
        kl = {}
        for cond in CONDS:
            kl[cond] = float(kl_stable(st["model"], st["hap"], st["util"].beta,
                                       p, qs[cond], eta=eta)["kl"])
        # ⚠ THE INSTRUMENT CHECK. Fails the CELL, not the run, so the report can say
        # which cells are trustworthy instead of losing every other cell's compute.
        d1_err = kl["D1_extremal"] - kl_d1
        return {"target": target_id, "label": label, "eta": eta, "kl": kl,
                "d1_expected": kl_d1, "d1_err": d1_err,
                "d1_ok": bool(abs(d1_err) <= D1_ATOL)}
    except Exception:
        return {"target": target_id, "label": label, "error": traceback.format_exc()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=30)
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--json", type=Path,
                    default=REPO / "logs" / "B_acc_arm2" / "arm2_chr21_ceiling.json")
    a = ap.parse_args(argv)

    grid = load_grid()
    targets = sorted(a.targets or grid)
    missing = sorted(set(targets) - set(grid))
    if missing:
        raise SystemExit(f"{missing} have no verified chr21 ceiling in {VERIFIED}")

    print(f"methodology: chr21-scope exact KL for all {len(CONDS)} decoy conditions, "
          f"arm 2, weights {WEIGHTS!r}, beta_share 1.0, K = 88")
    print(f"⚠ `eta` READ from {VERIFIED.relative_to(REPO)} (refined block), not re-solved.")
    print(f"⚠ `acc_ceiling` is Pinsker on the exact KL -- an UPPER BOUND on TV, not TV. "
          f"This closes the COLUMN-LABEL envelope, not the TV ranking of the conditions.")

    st = load_state()
    ids = {t: decoys_for(t, st) for t in targets}
    for t in targets:
        print(f"  {t}: " + ", ".join(f"{c}={ids[t][c]}" for c in CONDS))

    jobs = [(t, lab, grid[t][lab][0], grid[t][lab][1], ids[t])
            for t in targets for lab in sorted(grid[t], key=float)]
    n_cells = len(jobs) * len(CONDS)
    print(f"\n{len(targets)} targets x {len(jobs) // len(targets)} labels x "
          f"{len(CONDS)} conditions = {n_cells} exact KL evaluations "
          f"({len(jobs)} jobs, {a.workers} workers)")

    t0 = time.time()
    rows, errs = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(one_cell, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            (errs if "error" in r else rows).append(r)
            if i % 10 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} jobs  ({time.time() - t0:.0f}s)", flush=True)

    # ⚠ COMPLETENESS, ASSERTED. A script that filters results cannot tell "found
    # nothing" from "ran nothing" (CLAUDE.md, Execution).
    if errs:
        for e in errs[:3]:
            print(f"\n⚠ {e['target']} @ {e['label']}:\n{e['error']}")
        raise SystemExit(f"{len(errs)} of {len(jobs)} jobs FAILED; nothing is reported.")
    got = {(r["target"], r["label"]) for r in rows}
    want = {(j[0], j[1]) for j in jobs}
    if got != want:
        raise SystemExit(f"fan-out returned {len(got)} of {len(want)} cells; "
                         f"missing {sorted(want - got)}")

    bad = [r for r in rows if not r["d1_ok"]]
    if bad:
        for r in bad[:5]:
            print(f"⚠ {r['target']} @ {r['label']}: D1 kl {r['kl']['D1_extremal']:.12g} "
                  f"vs verified {r['d1_expected']:.12g} (err {r['d1_err']:+.3g})")
        raise SystemExit(
            f"INSTRUMENT CHECK FAILED on {len(bad)} of {len(rows)} cells: this run "
            f"cannot reproduce the D1 `kl_exact` its own `eta` was solved for, so its "
            f"real-decoy numbers are not trustworthy either. Nothing is reported.")
    print(f"\n✅ instrument check: D1 reproduces the verified `kl_exact` on all "
          f"{len(rows)} cells, worst |err| = {max(abs(r['d1_err']) for r in rows):.3g} "
          f"(atol {D1_ATOL:g})")

    rows.sort(key=lambda r: (r["target"], float(r["label"])))
    labels = sorted({r["label"] for r in rows}, key=float)

    print("\nExact chr21 `kl_exact` by decoy condition. Metrics are ROWS, "
          "`acc_ceiling` labels are COLUMNS.")
    print("\n| metric | " + " | ".join(labels) + " |")
    print("|---|" + "---|" * len(labels))
    for cond in CONDS:
        cells = []
        for lab in labels:
            v = [r["kl"][cond] for r in rows if r["label"] == lab]
            cells.append(f"{np.median(v):.6f}")
        print(f"| {cond} `kl_exact` median [n={len(targets)}] | " + " | ".join(cells) + " |")
    for cond in CONDS:
        cells = []
        for lab in labels:
            v = [ceiling_from_kl(r["kl"][cond])[1] for r in rows if r["label"] == lab]
            cells.append(f"{np.median(v):.4f}")
        print(f"| {cond} `acc_ceiling` median [n={len(targets)}] | " + " | ".join(cells) + " |")
    cells = []
    for lab in labels:
        v = [max(r["kl"][c] for c in REAL) / r["kl"]["D1_extremal"]
             for r in rows if r["label"] == lab]
        cells.append(f"{max(v):.6f}")
    print("| worst real/D1 `kl_exact` ratio [max over targets] * | " + " | ".join(cells) + " |")

    # The envelope itself, per cell, named by its arm and its scope.
    viol = [(r["target"], r["label"], c, r["kl"][c], r["kl"]["D1_extremal"])
            for r in rows for c in REAL
            if r["kl"][c] > r["kl"]["D1_extremal"] + 1e-9]
    print()
    if viol:
        print(f"⚠⚠ THE ENVELOPE DOES NOT HOLD in {len(viol)} of {len(rows) * len(REAL)} "
              f"(target, label, real condition) cells on chr21, arm 2:")
        for t, lab, c, kc, k1 in viol:
            print(f"    {t} @ {lab}: {c} kl {kc:.9g} > D1 {k1:.9g} "
                  f"(+{kc - k1:.3g}, ceiling {ceiling_from_kl(kc)[1]:.6f} vs "
                  f"{ceiling_from_kl(k1)[1]:.6f})")
        print("  ⚠ Every chr21 column label is D1's, so in these cells the label "
              "UNDERSTATES the bound for that real decoy. Report it; do not relabel.")
    else:
        print(f"✅ D1_extremal carries the HIGHEST exact `kl_exact` in all "
              f"{len(rows) * len(REAL)} (target, label, real condition) cells on chr21, "
              f"arm 2, so every chr21 column label bounds all three real decoys.")

    rank = {c: 0 for c in REAL}
    for r in rows:
        rank[max(REAL, key=lambda c: r["kl"][c])] += 1
    print("\nmost-divergent REAL condition, by cell: " +
          ", ".join(f"{c} {rank[c]}/{len(rows)}" for c in REAL))
    print("⚠ This ranks the three real decoys by exact KL, which is an upper-bound "
          "input and NOT the TV that the sampled attack ranks them on.")

    out = {"arm": "arm2", "scope": CHROM, "weights": WEIGHTS, "beta_share": 1.0,
           "k_states": 88, "conds": list(CONDS), "d1_atol": D1_ATOL,
           "decoy_ids": ids, "envelope_violations": len(viol), "rows": rows}
    a.json.write_text(json.dumps(out, indent=1))
    print(f"\nJSON -> {a.json.relative_to(REPO)}")
    print(f"elapsed {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
