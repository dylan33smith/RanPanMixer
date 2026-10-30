#!/usr/bin/env python
"""ARM 2's `acc_ceiling` GENOME-WIDE for all THREE decoy conditions.

⚠⚠ WHAT THIS EXISTS TO TEST. Every Arm 2 table labels its columns with the
`acc_ceiling` of the D1_extremal pair alone, and the label is a VALID bound for
D2_nearest and D3_random only because D1 is the maximum-TV decoy. Arm 1 ASSERTS
that per cell (`scripts/arm1_attack.py`, `cc <= c1 + 1e-9`). ARM 2 NEVER CHECKED IT,
and never checked it at genome scope in either arm. The property is an argument
about the construction -- D1 flips every called allele, so it maximises the utility
gap and therefore `kl` at any `eta` -- but `kl` is only approximately quadratic in
`eta`, the per-chromosome `share` splits the budget unevenly, and a decoy's own
chain differs from the target's, so "it must hold" is a prediction and this measures
it.

⚠ The `eta` grid is NOT re-solved here. It is read from the genome-wide sweep
(`logs/B_acc_arm2/genome_sweep/<target>.json`), which solved it against D1 at genome
scope -- the same calibration the chr21 grid uses, so the two are parallel.

⚠ D2/D3 are picked on GENOME-WIDE concordance, which need not agree with chr21's
nearest. That is the point: D2 is defined as the nearest, and at genome scope the
nearest is a genome-scope question.

⚠ `concordance` and `DECOY_SEED` are IMPORTED from `scripts/arm2_attack.py`, and the
D3 != D2 rejection is applied here the same way, so the two cannot drift silently.

    python scripts/arm2_genome_ceiling.py --chroms chr21 chr22 --workers 4
    python scripts/arm2_genome_ceiling.py --workers 6 --json out.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src"), str(REPO / "scripts")]

import numpy as np                                          # noqa: E402

from mechanism import hmm, utility                           # noqa: E402
from exact_kl_stable import AUTOSOMES, extremal, kl_stable   # noqa: E402
from arm2_attack import DECOY_SEED, concordance              # noqa: E402

CONDS = ("D1_extremal", "D2_nearest", "D3_random")
SWEEP = REPO / "logs" / "B_acc_arm2" / "genome_sweep"


def acc_ceiling(kl: float) -> float:
    """0.5 * (1 + min(1, sqrt(kl/2))) -- Pinsker, the project's one definition."""
    return 0.5 * (1.0 + min(1.0, math.sqrt(max(0.0, kl) / 2.0)))


def _target_path(t: str, chrom: str, ic: np.ndarray) -> np.ndarray:
    return np.load(REPO / "data" / "targets" / t / chrom / "path.npy")[ic]


def load_grid(targets: list[str]) -> dict:
    """Each target's genome-wide (acc_ceiling label -> eta), from the sweep."""
    grid = {}
    for t in targets:
        f = SWEEP / f"{t}.json"
        if not f.exists():
            raise SystemExit(f"no genome-wide sweep for {t} at {f} -- run the sweep first")
        d = json.loads(f.read_text())
        if d.get("scope") != "genome":
            raise SystemExit(f"{f} is scope {d.get('scope')!r}, not 'genome'")
        grid[t] = {"T": int(d["T"]), "shares": d["shares"],
                   "points": [(float(r["acc_ceiling"]), float(r["eta"]),
                               float(r["kl_exact"])) for r in d["aggregate"]]}
    return grid


def concordance_one_chrom(chrom: str, targets: list[str]) -> dict:
    """Per-chromosome (match, known) counts for every ordered target pair.

    Counted rather than averaged, because a genome-wide concordance is the ratio of
    GENOME totals and a mean of per-chromosome means would weight chr21 like chr1.
    """
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    paths = {t: _target_path(t, chrom, ic) for t in targets}
    out = {}
    for a in targets:
        for b in targets:
            if a >= b:
                continue
            pa, pb = paths[a], paths[b]
            known = (pa >= 0).all(axis=1) & (pb >= 0).all(axis=1)
            out[f"{a}|{b}"] = (int((pa[known] == pb[known]).all(axis=1).sum()),
                               int(known.sum()))
    return out


def kl_one_cell(chrom: str, t: str, grid: dict, decoys: dict,
                weights: str) -> dict:
    """ONE (chromosome, target)'s exact KL for every eta point and condition.

    ⚠ THE FAN-OUT UNIT IS (chromosome, target), NOT chromosome, and the reason is
    measured: one `kl_stable` call costs ~17.0 s at chr21's 332,521 positions, so a
    whole chromosome is ~18 calls per target and chr2 (5.59x chr21) would be a ~4.75 h
    CRITICAL PATH that no number of workers could shorten. Splitting by target drops
    the longest job to ~28 min. The per-chromosome state load this duplicates was
    measured at 0.1 s against 17 s per call, so it is noise.
    """
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=chrom)
    support = np.load(d / "support.npy")
    out = {}
    beta = utility.beta_v(support, in_chain=ic, weights=weights,
                          share=float(grid[t]["shares"][chrom]), chrom=chrom).beta
    p = _target_path(t, chrom, ic)
    qs = {"D1_extremal": extremal(p)}
    for cond in ("D2_nearest", "D3_random"):
        qs[cond] = _target_path(decoys[t][cond], chrom, ic)
    for label, eta, _ in grid[t]["points"]:
        for cond in CONDS:
            key = f"{t}|{label!r}|{cond}"
            if eta == 0.0:
                out[key] = 0.0            # Q_p = Q_q identically; no work to do
                continue
            out[key] = float(kl_stable(m, hap, beta, p, qs[cond], eta=eta)["kl"])
    del hap, m, beta
    return {"chrom": chrom, "target": t, "t_chain": int(ic.sum()), "kl": out}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=list(AUTOSOMES),
                    help="subset, for timing a single chromosome before the full run")
    ap.add_argument("--targets", nargs="+", default=None,
                    help="default: every target with a genome-wide sweep on disk")
    ap.add_argument("--weights", choices=("uniform", "inv_support"), default="uniform")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args()

    targets = a.targets or sorted(f.stem for f in SWEEP.glob("*.json"))
    if not targets:
        raise SystemExit(f"no sweep JSONs in {SWEEP}")
    grid = load_grid(targets)
    partial = len(a.chroms) != len(AUTOSOMES)

    print(f"ARM 2 GENOME-WIDE `acc_ceiling`, ALL THREE CONDITIONS. "
          f"{len(targets)} targets x {len(grid[targets[0]]['points'])} labels x "
          f"{len(CONDS)} conditions, {len(a.chroms)} chromosome(s), "
          f"weights {a.weights}, {a.workers} workers.")
    if partial:
        print(f"⚠⚠ PARTIAL: {' '.join(a.chroms)} only, so every `kl` below is a PARTIAL "
              f"SUM and no ceiling here is the genome-wide ceiling. Timing run only.")

    # ---- pass 1: genome-wide concordance, then D2/D3 with the D3 != D2 rejection
    t0 = time.time()
    match, known = {}, {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(concordance_one_chrom, c, targets): c for c in a.chroms}
        for fu in as_completed(futs):
            for k, (mm, kk) in fu.result().items():
                match[k] = match.get(k, 0) + mm
                known[k] = known.get(k, 0) + kk
    pairs = {k: match[k] / known[k] for k in match if known[k]}
    assert len(pairs) == len(targets) * (len(targets) - 1) // 2, (
        f"concordance fan-out produced {len(pairs)} pairs, expected "
        f"{len(targets) * (len(targets) - 1) // 2} -- a missing chromosome would look "
        f"like a different nearest neighbour, not like an error")

    decoys = {}
    for t in targets:
        others = [o for o in targets if o != t]
        conc = {o: pairs[f"{min(t, o)}|{max(t, o)}"] for o in others}
        nearest = max(conc, key=conc.get)
        # ⚠ SAME RULE AS `arm2_attack.pick_decoys`: D3 must never BE D2, and the
        # rejection (not a narrowed pool) keeps a non-colliding draw where it landed.
        if len(others) < 2:
            raise SystemExit(f"{t}: pool of {len(others)}, D3 cannot differ from D2")
        stable = int.from_bytes(__import__("hashlib").sha256(t.encode()).digest()[:4],
                                "big")
        rng = np.random.default_rng([DECOY_SEED, stable])
        rand = others[int(rng.integers(len(others)))]
        while rand == nearest:
            rand = others[int(rng.integers(len(others)))]
        decoys[t] = {"D2_nearest": nearest, "D3_random": rand,
                     "concordance": conc, "conc_nearest": conc[nearest]}
    print(f"  concordance pass: {time.time() - t0:.1f}s")
    for t in targets:
        print(f"    {t}: D2 {decoys[t]['D2_nearest']} "
              f"(conc {decoys[t]['conc_nearest']:.6f})  D3 {decoys[t]['D3_random']}")

    # ---- pass 2: exact KL per chromosome, summed
    t1 = time.time()
    tot: dict[str, float] = {}
    done, per_chrom_tchain = [], {}
    jobs = [(c, t) for c in a.chroms for t in targets]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(kl_one_cell, c, t, grid, decoys, a.weights): (c, t)
                for c, t in jobs}
        for fu in as_completed(futs):
            r = fu.result()
            for k, v in r["kl"].items():
                tot[k] = tot.get(k, 0.0) + v
            done.append((r["chrom"], r["target"]))
            # ⚠ ONCE per chromosome: every target of a chromosome reports the SAME
            # t_chain, and summing per job would multiply the genome length by 10.
            per_chrom_tchain[r["chrom"]] = r["t_chain"]
            print(f"  {r['chrom']}/{r['target']}: t_chain {r['t_chain']:,} "
                  f"({len(done)}/{len(jobs)})", flush=True)
    assert sorted(done) == sorted(jobs), (
        f"KL fan-out ran {len(done)} of {len(jobs)} (chromosome, target) jobs -- a "
        f"filtered result cannot tell 'found nothing' from 'ran nothing'")
    t_chain = sum(per_chrom_tchain.values())
    print(f"  KL pass: {time.time() - t1:.1f}s, t_chain total {t_chain:,}")

    # ---- the property under test
    labels = [lab for lab, _, _ in grid[targets[0]]["points"]]
    rows, viol = [], []
    for t in targets:
        for lab, eta, kl_d1_sweep in grid[t]["points"]:
            cs = {c: acc_ceiling(tot[f"{t}|{lab!r}|{c}"]) for c in CONDS}
            kl1 = tot[f"{t}|{lab!r}|D1_extremal"]
            row = {"target": t, "label": lab, "eta": eta, "t_chain": t_chain,
                   "kl": {c: tot[f"{t}|{lab!r}|{c}"] for c in CONDS},
                   "acc_ceiling": cs,
                   "kl_d1_sweep": kl_d1_sweep,
                   "d1_is_max": all(cs["D1_extremal"] >= cs[c] - 1e-9 for c in CONDS)}
            rows.append(row)
            if eta > 0.0 and not row["d1_is_max"]:
                viol.append(row)
            # D1 here must reproduce the sweep's own kl, or the two disagree about
            # the same quantity and one of them is wrong.
            if eta > 0.0 and kl_d1_sweep > 0:
                rel = abs(kl1 - kl_d1_sweep) / kl_d1_sweep
                row["d1_vs_sweep_rel"] = rel

    print("\n| target | label | eta | `acc_ceiling` D1 | D2 | D3 | D1 is max |")
    print("|---|---|---|---|---|---|---|")
    for r in rows:
        c = r["acc_ceiling"]
        print(f"| {r['target']} | {r['label']:.2f} | {r['eta']:.3f} | "
              f"{c['D1_extremal']:.6f} | {c['D2_nearest']:.6f} | "
              f"{c['D3_random']:.6f} | {'yes' if r['d1_is_max'] else '⚠⚠ NO'} |")

    rels = [r["d1_vs_sweep_rel"] for r in rows if "d1_vs_sweep_rel" in r]
    if rels and partial:
        # ⚠ THE COMPARISON IS MEANINGLESS ON A PARTIAL CHROMOSOME SET and must not be
        # printed as though it were a disagreement. `kl` is a SUM over positions while
        # the sweep's `kl_exact` is the GENOME total, so a chr21-only run differs by
        # ~(1 - 332,521/22,134,562) = 0.985 BY CONSTRUCTION. The first partial run
        # printed 9.883e-01 under a warning saying nothing here is interpretable,
        # which is exactly the wrong reading of a correct number.
        print(f"\n⚠ D1 `kl` vs the sweep's `kl_exact`: NOT COMPARED — this is a "
              f"PARTIAL chromosome set, so these `kl` are partial sums and the "
              f"sweep's are genome totals. The observed ratio would be the position "
              f"share, not an error. (Partial/total differed by {max(rels):.3e}, "
              f"against an expected ~{1 - t_chain / grid[targets[0]]['T']:.3e}.)")
    elif rels:
        print(f"\n⚠ D1 `kl` here vs the sweep's own `kl_exact`: worst relative "
              f"difference {max(rels):.3e} over {len(rels)} tilted cells. A large "
              f"value means this script and the sweep disagree about D1, and then "
              f"NOTHING here is interpretable.")
    tested = [r for r in rows if r["eta"] > 0.0]
    print(f"\n{'⚠⚠ ' if viol else ''}D1 IS THE MAXIMUM-TV DECOY IN "
          f"{len(tested) - len(viol)} OF {len(tested)} TILTED CELLS"
          f"{'' if viol else ' — the label is a valid bound for D2 and D3 genome-wide'}.")
    for r in viol:
        c = r["acc_ceiling"]
        print(f"  ⚠⚠ {r['target']} at {r['label']:.2f}: D1 {c['D1_extremal']:.6f} < "
              f"D2 {c['D2_nearest']:.6f} / D3 {c['D3_random']:.6f} — THE COLUMN LABEL "
              f"IS NOT CONSERVATIVE FOR THIS CELL.")
    print("\n⚠ The 0.50 label is eta = 0, where Q_p = Q_q for EVERY decoy, so all "
          "three ceilings are exactly 0.5 and the cell is a harness check.")
    if partial:
        print("⚠⚠ PARTIAL CHROMOSOME SET — not a genome-wide result.")

    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm2", "scope": "genome" if not partial else "partial",
             "chroms": sorted(per_chrom_tchain), "t_chain": t_chain, "weights": a.weights,
             "targets": targets, "labels": labels, "conds": list(CONDS),
             "decoys": {t: {k: v for k, v in decoys[t].items()} for t in targets},
             "rows": rows, "n_violations": len(viol)}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 1 if viol else 0


if __name__ == "__main__":
    raise SystemExit(main())
