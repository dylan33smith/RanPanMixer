#!/usr/bin/env python
"""Stage both tools' releases so PanMixer's OWN privacy evaluators run on them unchanged.

⚠⚠ WHY A MIRROR ROOT. `external/PanMixer/constants.py:4` sets
`BASE_PATH = os.path.dirname(os.path.abspath(__file__))`, with `STARTING_DATA_PATH` and
`EXPERIMENT_PATH` derived from it and NO environment override anywhere in their tree.
`abspath` does not resolve symlinks, so BASE_PATH is whatever path string Python used to
reach `constants.py`. Copying `constants.py` into a writable root and symlinking `tools`
and `starting_data` back to theirs therefore redirects EXPERIMENT_PATH into our root while
their code stays byte-identical and untouched. Nothing under `external/` is written.

⚠ THEIR `main.py` IS NOT USABLE: every evaluator verb builds a SLURM array and shells out
to `sbatch`, which is not installed here. The per-task CLI
`<evaluator>.py <experiment> <chromosome> <row_id>` is the real entry point.

⚠⚠ ONE EXPERIMENT PER EVALUATOR, because `gap_score.py` and `diploid_gap_score.py` write
the SAME filenames (`gap_score_all.json`, `genotypes_scores.npy`,
`genotypes_scores_self.npy`) into the SAME row directory. Sharing a row would let the
second silently overwrite the first.

⚠ `sampled_gap_score.py` is NOT runnable here: it reads the PRE-stacker object from
`starting_data/chr{C}/subjects/{id}/new_haplotypes.npy`, and `starting_data/chr21/subjects`
does not exist. Three evaluators run, not four.

⚠ THEIR side is taken from their own `obfuscated.vcf.gz`, not from the `new_haplotypes.npy`
in their sweep output, because that file is the pre-stacker resampled array while the VCF
is the released object. The VCF extraction is already validated: it reproduces their
reported `changed_alleles` and `total_alleles` EXACTLY at all 8 capacities.

    python scripts/pm_eval_setup.py
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PM = REPO / "external" / "PanMixer"
PM_SWEEP = Path("/data/ds85/RanPanMixer/runs/capacity_sweep")
ROOT = Path("/data/ds85/RanPanMixer/runs/pm_eval")
OURS = REPO / "data" / "panmixer_integration" / "arm1" / "HG00438"
TARGET = "HG00438"
CHROM = 21
TAUS = ("0p00", "0p05", "0p10", "0p25", "0p50", "0p90", "0p99")
# one experiment per evaluator -- they collide on output filenames
EXPERIMENTS = {"1": "gap_score", "2": "diploid_gap_score", "3": "MIA_privacy"}

import numpy as np                                        # noqa: E402


def their_release(cap_dir: Path, n_sites: int) -> np.ndarray:
    """Their FINAL released object, from their own obfuscated VCF."""
    out = subprocess.run(
        ["bcftools", "query", "-s", TARGET, "-f", "[%GT]\\n",
         str(cap_dir / "chr21" / "obfuscated.vcf.gz")],
        capture_output=True, text=True, check=True)
    g = out.stdout.split()
    assert len(g) == n_sites, f"{cap_dir.name}: {len(g)} genotypes vs {n_sites} sites"
    a = np.full((n_sites, 2), -1, dtype=np.int16)
    for i, s in enumerate(g):
        parts = s.replace("/", "|").split("|")
        for j in (0, 1):
            x = parts[j] if j < len(parts) else "."
            if x not in (".", ""):
                a[i, j] = int(x)
    return a


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    # release directory labels in `releases.tau_dir` form, so a second mirror can hold
    # the exploratory high-`tau` rows (0p999, 0p9999, 0p99999) beside their own copy
    # of THEIR capacity sweep -- the two sides must always be staged together
    ap.add_argument("--taus", nargs="+", default=list(TAUS))
    a = ap.parse_args(argv)
    root = a.root

    pos = np.load(PM / "starting_data" / "chr21" / "pangenome_positions.npy")
    na = np.load(PM / "starting_data" / "chr21" / "num_alleles.npy")
    n_sites = len(pos)

    # ---- the mirror. constants.py is a COPY (its directory becomes BASE_PATH);
    #      tools and starting_data are symlinks back to theirs, read-only.
    root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PM / "constants.py", root / "constants.py")
    for name in ("tools", "starting_data"):
        link = root / name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to((PM / name).resolve())
    print(f"mirror root: {root}")
    print(f"  constants.py copied, tools and starting_data symlinked (read-only)")

    # ---- the 15 release objects, in one row order shared by every experiment
    rows = []
    for t in a.taus:
        rows.append({"kind": "ours", "label": f"tau_{t}",
                     "array": np.load(OURS / f"tau_{t}" / "chr21" / "new_haplotypes.npy")})
    for cd in sorted(PM_SWEEP.glob("cap*")):
        st = json.loads((cd / "chr21" / "stats.json").read_text())
        rows.append({"kind": "theirs", "label": f"capacity_{st['capacity']:g}",
                     "array": their_release(cd, n_sites),
                     "capacity": float(st["capacity"]),
                     "their_changed": int(st["changed_alleles"]),
                     "their_total": int(st["total_alleles"])})

    # ⚠ every array must satisfy THEIR OWN pre-flight before it is staged:
    # tools/common/verify_genotypes.py asserts code < num_alleles[site].
    for i, r in enumerate(rows):
        arr = r["array"]
        assert arr.shape == (n_sites, 2) and arr.dtype == np.int16, \
            f"row {i} ({r['label']}): shape {arr.shape} dtype {arr.dtype}"
        bad = int(((arr >= na[:, None]) & (arr >= 0)).sum())
        assert bad == 0, (
            f"row {i} ({r['label']}): {bad} allele code(s) out of range for their "
            f"encoding — their verify_genotypes would reject this file")
    print(f"  {len(rows)} release objects staged, all passing their verify_genotypes "
          f"contract ({sum(r['kind'] == 'ours' for r in rows)} ours, "
          f"{sum(r['kind'] == 'theirs' for r in rows)} theirs)")

    # ---- one experiment per evaluator, identical rows in each
    import csv
    for exp, name in EXPERIMENTS.items():
        ed = root / "experiments" / f"exp_{exp}"
        (ed / "data" / f"chr{CHROM}").mkdir(parents=True, exist_ok=True)
        with (ed / f"data_chr{CHROM}.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["subject", "label", "kind"])
            for r in rows:
                w.writerow([TARGET, r["label"], r["kind"]])
        for i, r in enumerate(rows):
            rd = ed / "data" / f"chr{CHROM}" / str(i)
            rd.mkdir(parents=True, exist_ok=True)
            np.save(rd / "new_haplotypes.npy", r["array"])
        print(f"  exp_{exp} ({name}): {len(rows)} rows")

    man = {"root": str(root), "target": TARGET, "chrom": CHROM,
           "experiments": EXPERIMENTS,
           "rows": [{"row_id": i, "kind": r["kind"], "label": r["label"],
                     **({k: r[k] for k in ("capacity", "their_changed", "their_total")}
                        if r["kind"] == "theirs" else {})}
                    for i, r in enumerate(rows)]}
    (root / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(f"\nmanifest -> {root / 'manifest.json'}")
    print(f"⚠ run each evaluator from {root} with PYTHONDONTWRITEBYTECODE=1 and the "
          f"`panmixer` env — a 3.11 import would otherwise write .pyc under external/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
