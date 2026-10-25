#!/usr/bin/env python
"""Stage our ARM 1 releases so PanMixer's BEAGLE RECONSTRUCTION ATTACK runs on them unchanged.

⚠⚠ WHICH BEAGLE PATH IS THE REAL ONE. Their tree contains two. `tools/beagle/beagle_stats.py`
(`wged`, `pmi`) is **DEAD CODE and cannot run as shipped** — audited 2026-10-14: it uses
`PANGENOME_SUBJECTS` and `EXPERIMENT_PATH` without importing them (NameError at its lines 89
and 128), its `main()` reads `new_haplotypes_beagle.vcf.gz.vcf.gz` while `beagle_refinement.py`
writes `beagle_refined.vcf.gz`, its paths omit the `chr{N}` component, it compares the Beagle
array to the pangenome array by INDEX when the two axes differ (950,597 Beagle records against
340,824 pangenome records, measured), and `main.py` registers a `beagle_stats` subparser it
never dispatches. Its `pmi` term also needs `pop_blocks_dict.json` and `pmi/{subject}_all.npy`,
neither of which is in the shipped tree.
**The LIVE, published path is `beagle_refinement.py` -> `accuracy_stats.py`**, whose headline
metric is `reconstruction_accuracy`, and `accuracy_stats.py` aligns Beagle to the pangenome axis
BY POSITION (`align_beagle_to_positions`) rather than by index. That is the one we run.

WHAT THE ATTACK IS. Beagle is a statistical phasing and imputation tool. The adversary takes the
released VCF, hands it to Beagle together with a 1000 Genomes reference panel, and lets Beagle's
population model pull the genotypes back toward what a real human haplotype would look like.
`reconstruction_accuracy` then asks, over exactly the cells the release CHANGED:

    among cells where truth, release and Beagle output are all called, and the release differs
    from the truth, how often did Beagle put the TRUE allele back?

0.0 means the attack never recovers a changed cell; 1.0 means it recovers all of them.

⚠⚠ WE HAVE NO `xsol`, AND THE SUBSTITUTION IS DECLARED HERE. `accuracy_stats.py` restricts every
metric to the (site, strand) cells PanMixer's OPTIMIZER selected for obfuscation, read from
`xsol.npy` of shape (n_blocks, 2) in `blocks_dict.json` key order. Our mechanism has no selection
step at all — a path release may alter any position it covers — so the faithful analogue is an
ALL-ONES `xsol`: every block selected. `xsol_expanded` is then all True and each metric is
restricted by its variant-class subset alone. ⚠ This is the one substitution in the whole run,
it is an input we supply rather than a change to their code, and it makes OUR denominator the
larger one: their numbers are computed over a selected subset of cells and ours over every
released cell.

⚠ 39 OF OUR 44 TARGETS ARE IN THE BEAGLE REFERENCE PANEL. Their `beagle_refinement.py` passes no
`excludesamples=`, so for those 39 the panel contains the target's own true row and Beagle can
copy from it. That is a property of their published protocol, not a defect of it, and it is why
`--exclude-self` exists below as a labelled diagnostic: the 5 out-of-panel targets and the
exclude-self runs are what separate what the RELEASE leaks from what the PANEL leaks.

⚠ Nothing under `external/` is written. The mirror root holds a COPY of `constants.py` (their
`BASE_PATH` is `dirname(abspath(__file__))` with no env override, and `abspath` does not resolve
symlinks, so copying that one file redirects `EXPERIMENT_PATH` into our tree) with `tools` and
`starting_data` symlinked back read-only.

    python scripts/beagle_attack_setup.py
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
PM = REPO / "external" / "PanMixer"
OURS = REPO / "data" / "panmixer_integration" / "arm1"
ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
CHROM = 21
EXPERIMENT = "1"
TAUS = ("0p00", "0p05", "0p10", "0p25", "0p50", "0p90", "0p99")


def point_label(label: str) -> tuple:
    """A `--taus` entry -> (directory name, knob, numeric value).

    ⚠ ACCEPTS THREE FORMS because the release tree now has two namespaces:
      "0p99"      legacy bare suffix -> tau_0p99, knob "tau",  0.99
      "tau_0p99"  explicit           -> tau_0p99, knob "tau",  0.99
      "eta_250"   the eta namespace  -> eta_250,  knob "eta",  250.0
    The old code hardcoded an "tau_" prefix, so an eta label became `tau_eta_0`
    and the stage aborted on a missing release -- measured 2026-10-25.
    """
    if label.startswith("eta_"):
        return label, "eta", float(label[4:].replace("p", "."))
    name = label if label.startswith("tau_") else f"tau_{label}"
    return name, "tau", float(name[4:].replace("p", "."))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--targets", type=int, default=0,
                    help="limit to the first N targets (0 = all); for a smoke test")
    ap.add_argument("--taus", nargs="+", default=list(TAUS),
                    help="release directory labels, e.g. 0p999 (releases.tau_dir form)")
    a = ap.parse_args(argv)
    root = a.root

    # ---- the mirror
    root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PM / "constants.py", root / "constants.py")
    for name in ("tools", "starting_data", "downloaded_tools"):
        link = root / name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to((PM / name).resolve())
    print(f"mirror root: {root}")

    sub = np.load(PM / "starting_data" / "chr21" / "pangenome_subjects.npy")
    na = np.load(PM / "starting_data" / "chr21" / "num_alleles.npy")
    n_sites = len(np.load(PM / "starting_data" / "chr21" / "pangenome_positions.npy"))
    blocks = json.load(open(PM / "starting_data" / "chr21" / "blocks_dict.json"))
    n_blocks = len(blocks)
    print(f"  {n_sites:,} sites, {n_blocks:,} blocks -> xsol shape ({n_blocks}, 2), ALL ONES")

    targets = sorted(d.name for d in OURS.iterdir() if d.is_dir())
    if a.targets:
        targets = targets[: a.targets]
    # ⚠ their evaluators look the subject up in pangenome_subjects.npy and would
    # IndexError on a miss; check up front rather than 300 rows in.
    missing = [t for t in targets if t not in set(sub.tolist())]
    if missing:
        raise SystemExit(f"not in their pangenome_subjects.npy: {missing}")

    xsol = np.ones((n_blocks, 2), dtype=np.int64)
    rows = []
    ed = root / "experiments" / f"exp_{EXPERIMENT}"
    (ed / "data" / f"chr{CHROM}").mkdir(parents=True, exist_ok=True)

    for t in targets:
        for tau in a.taus:
            pname, knob, pval = point_label(tau)
            src = OURS / t / pname / f"chr{CHROM}" / "new_haplotypes.npy"
            if not src.exists():
                raise SystemExit(f"missing release: {src}")
            arr = np.load(src)
            assert arr.shape == (n_sites, 2) and arr.dtype == np.int16, \
                f"{t}/{tau}: shape {arr.shape} dtype {arr.dtype}"
            # ⚠ their own pre-flight: an allele code must be < num_alleles[site].
            bad = int(((arr >= na[:, None]) & (arr >= 0)).sum())
            assert bad == 0, f"{t}/{tau}: {bad} out-of-range allele code(s)"
            i = len(rows)
            rd = ed / "data" / f"chr{CHROM}" / str(i)
            rd.mkdir(parents=True, exist_ok=True)
            np.save(rd / "new_haplotypes.npy", arr)
            np.save(rd / "xsol.npy", xsol)
            rows.append({"row_id": i, "subject": t, "tau": tau,
                         "knob": knob, "point_value": pval,
                         "tau_value": pval if knob == "tau" else None})

    with (ed / f"data_chr{CHROM}.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subject", "tau"])
        for r in rows:
            w.writerow([r["subject"], r["tau"]])

    man = {"root": str(root), "experiment": EXPERIMENT, "chrom": CHROM,
           "n_sites": n_sites, "n_blocks": n_blocks,
           "xsol": "ALL ONES -- declared substitution, see the module docstring",
           "taus": list(a.taus), "targets": targets, "rows": rows}
    (root / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(f"  staged {len(rows)} rows ({len(targets)} targets x {len(a.taus)} tau)")
    print(f"manifest -> {root / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
