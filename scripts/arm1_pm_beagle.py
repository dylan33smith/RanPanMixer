#!/usr/bin/env python
"""Beagle reconstruction on PANMIXER'S OWN releases, so the Arm 1 table has a cross-tool
Beagle number instead of `n/a`.

⚠⚠ THE `xsol` MASK IS A COMPARABILITY DECISION, NOT A MISSING FILE, and validating one row
before running eight is what surfaced it. `accuracy_stats.py:162` scores every metric ONLY
over cells inside `xsol`-selected blocks:

    xsol_expanded[block_index, j] = True  for each selected (block, haplotype)
    subsets_for_just_obfuscated[name] = subset AND xsol_expanded

Our releases were scored with an ALL-ONES `xsol` -- the declared substitution in
`scripts/beagle_attack_setup.py`, because a path release has no selection step and may
alter any position it covers. PanMixer's releases carry a REAL `xsol` from their optimiser
(14,342 blocks selected at `capacity` 0.01, rising to 177,865). Scoring their releases
under their mask while ours used all-ones would put DIFFERENT DENOMINATORS on the two sides
of one table, which is exactly the defect this comparison exists to avoid. So both:

    exp_6   all-ones mask, matched to ours   -> THE COMPARABLE COLUMN
    exp_7   their own xsol, native protocol  -> reported separately, NOT comparable

⚠ `reconstruction_accuracy` is defined over CHANGED cells only, so its denominator differs
between the tools by construction (597 SNP cells at `capacity` 0.01 against ~81,851 for
ours). That ratio belongs to the mechanisms, not to the attack, so the changed-cell count
is printed beside every accuracy and the bare numbers are never compared.

    python scripts/arm1_pm_beagle.py --stage     # stage exp_6 + exp_7
    python scripts/arm1_pm_beagle.py             # tabulate (after the runs)

The runs themselves go through the shared per-row runner, which needs `--experiment` and
`--subject` because PanMixer's rows are NOT in our manifest (whose row_ids 0..307 are ours):

    python scripts/beagle_attack_one.py --row N --experiment 6 --subject HG00438
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
SWEEP = Path("/data/ds85/RanPanMixer/runs/capacity_sweep")
PM = REPO / "external" / "PanMixer" / "starting_data" / "chr21"
CAPS_DIR = ["0p01", "0p05", "0p1", "0p2", "0p3", "0p5", "0p7", "0p9"]
CAPS = [c.replace("p", ".") for c in CAPS_DIR]
TARGET = "HG00438"
KEYS = ("obfuscated_accuracy", "beagle_accuracy", "major_allele_accuracy",
        "reconstruction_accuracy")


def stage() -> None:
    n_blocks = len(json.loads((PM / "blocks_dict.json").read_text()))
    ones = np.ones((n_blocks, 2), dtype=np.int64)
    for exp, own in ((6, False), (7, True)):
        ed = ROOT / "experiments" / f"exp_{exp}"
        (ed / "data" / "chr21").mkdir(parents=True, exist_ok=True)
        with (ed / "data_chr21.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["subject", "capacity"])
            for c in CAPS:
                w.writerow([TARGET, c])
        for i, c in enumerate(CAPS_DIR):
            rd = ed / "data" / "chr21" / str(i)
            rd.mkdir(parents=True, exist_ok=True)
            src = (SWEEP / f"cap{c}" / "chr21" / "new_haplotypes.npy").resolve()
            dst = rd / "new_haplotypes.npy"
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src)
            if own:
                xs = (SWEEP / f"cap{c}" / "chr21" / "xsol.npy").resolve()
                xd = rd / "xsol.npy"
                if xd.is_symlink() or xd.exists():
                    xd.unlink()
                xd.symlink_to(xs)
            else:
                np.save(rd / "xsol.npy", ones)
        print(f"exp_{exp} staged: 8 capacities, "
              f"{'THEIR OWN xsol' if own else f'ALL-ONES xsol ({n_blocks:,} blocks)'}")


def _changed_snp_cells() -> list[int]:
    orig = np.load(PM / "pangenome.npy")
    subs = np.load(PM / "pangenome_subjects.npy").tolist()
    o = orig[subs.index(TARGET)]
    snp = np.zeros(orig.shape[1], bool)
    snp[np.load(PM / "pangenome_mask.npy")] = True
    s2 = np.column_stack((snp, snp))
    out = []
    for c in CAPS_DIR:
        rel = np.load(SWEEP / f"cap{c}" / "chr21" / "new_haplotypes.npy")
        out.append(int((((o >= 0) & (rel >= 0) & (o != rel)) & s2).sum()))
    return out


def _load(exp: int, i: int):
    # ⚠ the runner leaves the non-exclude-self result under its original name; the
    # _protocol copy only exists where an exclude-self pass stashed it.
    for nm in ("accuracy_results.json", "accuracy_results_protocol.json"):
        f = ROOT / f"experiments/exp_{exp}/data/chr21" / str(i) / nm
        if f.exists():
            return json.loads(f.read_text())["snp_only"]
    return None


def report() -> int:
    changed = _changed_snp_cells()
    missing = []
    for exp, lab in ((6, "MATCHED all-ones mask -- COMPARABLE TO OURS"),
                     (7, "THEIR OWN xsol -- native protocol, NOT comparable")):
        rows = {k: [] for k in KEYS}
        for i in range(8):
            d = _load(exp, i)
            if d is None:
                missing.append((exp, i))
            for k in KEYS:
                rows[k].append(d[k] if d else float("nan"))
        print(f"\n### PanMixer's own releases, {lab}\n")
        print("| variable | " + " | ".join(CAPS) + " |")
        print("|---|" + "---|" * 8)
        for k in KEYS[:3]:
            print(f"| `{k}` | " + " | ".join(f"{x:.6f}" for x in rows[k]) + " |")
        diff = [a - b for a, b in zip(rows["beagle_accuracy"],
                                      rows["major_allele_accuracy"])]
        print("| **beagle − major** | " + " | ".join(f"{x:+.6f}" for x in diff) + " |")
        print("| beagle beats baseline? | "
              + " | ".join("**YES**" if x > 0 else "no" for x in diff) + " |")
        print(f"| `reconstruction_accuracy` | "
              + " | ".join(f"{x:.6f}" for x in rows["reconstruction_accuracy"]) + " |")
        print("| changed cells (snp_only) | " + " | ".join(f"{c:,}" for c in changed) + " |")
    if missing:
        raise SystemExit(f"\n❌ INCOMPLETE: {len(missing)} missing results {missing}. "
                         f"A table built from a partial fan-out cannot tell "
                         f"'found nothing' from 'ran nothing'.")
    print("\n✅ 16 of 16 results present.")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", action="store_true")
    a = ap.parse_args(argv)
    if a.stage:
        stage()
        return 0
    return report()


if __name__ == "__main__":
    raise SystemExit(main())
