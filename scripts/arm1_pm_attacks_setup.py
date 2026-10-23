#!/usr/bin/env python
"""Stage all 308 ARM 1 releases for PanMixer's three runnable LINKAGE/MIA evaluators.

This extends the Beagle attack's mirror root with one experiment per evaluator, over the SAME
308 rows (44 targets x 7 `tau`), so every attack in their paper that can run here runs on one
staging with one row order.

⚠⚠ ONE EXPERIMENT PER EVALUATOR, because `gap_score.py` and `diploid_gap_score.py` write the
SAME output filenames (`gap_score_all.json`, `genotypes_scores.npy`, `genotypes_scores_self.npy`)
into the row directory. Sharing a row would let the second silently overwrite the first — the
same hazard already recorded for the 15-row run.

  exp_1 = Beagle reconstruction (already staged by scripts/beagle_attack_setup.py)
  exp_2 = gap_score            (rarity-weighted GENOTYPE linkage)
  exp_3 = diploid_gap_score    (the HAPLOTYPE version)
  exp_4 = MIA_privacy          (membership inference)

⚠ `new_haplotypes.npy` is SYMLINKED to exp_1's copy rather than duplicated: 308 arrays x 1.4 MB
per experiment is 420 MB of identical bytes per evaluator, and three copies of one array is three
chances for them to drift apart. The evaluators only read it.

⚠ `sampled_gap_score.py` is NOT staged and cannot run: it reads a PRE-stacker object from
`starting_data/chr{C}/subjects/{id}/new_haplotypes.npy`, and that directory is absent from the
shipped tree. Three evaluators, never four.

    python scripts/arm1_pm_attacks_setup.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
CHROM = 21
EVAL = {"2": "gap_score", "3": "diploid_gap_score", "4": "MIA_privacy"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    # exp_5 = PanMixer's af_loss/ld_loss (EXPLORATORY, decision H), staged by hand on
    # 2026-10-14: the same symlinks, and a `capacity` column holding the literal "n/a"
    ap.add_argument("--af-ld", action="store_true", help="also stage exp_5")
    a = ap.parse_args(argv)
    root = a.root

    man = json.loads((root / "manifest.json").read_text())
    rows = man["rows"]
    src_base = root / "experiments" / "exp_1" / "data" / f"chr{CHROM}"
    if not src_base.exists():
        raise SystemExit(f"{src_base} missing -- run scripts/beagle_attack_setup.py first")

    for exp, name in (EVAL | ({"5": "af_ld"} if a.af_ld else {})).items():
        ed = root / "experiments" / f"exp_{exp}"
        (ed / "data" / f"chr{CHROM}").mkdir(parents=True, exist_ok=True)
        with (ed / f"data_chr{CHROM}.csv").open("w", newline="") as f:
            w = csv.writer(f)
            cap = exp == "5"
            w.writerow(["subject", "tau"] + (["capacity"] if cap else []))
            for r in rows:
                w.writerow([r["subject"], r["tau"]] + (["n/a"] if cap else []))
        n_link = 0
        for r in rows:
            i = r["row_id"]
            rd = ed / "data" / f"chr{CHROM}" / str(i)
            rd.mkdir(parents=True, exist_ok=True)
            src = (src_base / str(i) / "new_haplotypes.npy").resolve()
            if not src.exists():
                raise SystemExit(f"row {i}: {src} missing")
            dst = rd / "new_haplotypes.npy"
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src)
            n_link += 1
        print(f"  exp_{exp} ({name}): {n_link} rows staged")

    man["evaluators"] = {"1": "beagle_reconstruction", **EVAL}
    (root / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(f"\n{len(rows)} rows x {len(EVAL)} evaluators = {len(rows)*len(EVAL)} jobs to run")
    print(f"manifest updated -> {root / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
