#!/usr/bin/env python
"""Stage BOTH tools' chr21 releases for PanMixer's `af_loss`/`ld_loss` ON ONE BUILD.

⚠⚠ THIS EXISTS BECAUSE THE EARLIER af/ld RUN WAS CROSS-BUILD. `logs/B_eta_arm1/af_ld.md`
compared our releases (measured against PanMixer's vg-1.36 `starting_data`) with their own
published numbers, so decision F was unsatisfied and the comparison carried a build
confound on top of an operating-point confound. Both halves now sit on OUR rebuilt
GRCh38 tree -- 340,849 records -- and the evaluators run UNMODIFIED against it.

WHAT IS STAGED
  exp_5  OURS    44 targets x 7 `acc_ceiling` (0.50..0.80) = 308 rows
  exp_6  THEIRS  44 subjects x 8 `capacity` (0.01..0.9)    = 352 rows
Both read `new_haplotypes.npy` of shape (n_sites, 2) int16 and both are measured against
the SAME `pangenome.npy`, so a difference between the two tables is a difference between
the two mechanisms and nothing else.

⚠⚠ THE AXIS IS THE WHOLE POINT, SO IT IS ASSERTED PER ROW. `data/panmixer_integration/`
(the 2026-10-25 preview) holds 7 HG00438 cells on PanMixer's ORIGINAL 340,824 axis; this
script's source is `data/panmixer_integration_ourbuild/`, on our 340,849 axis. Pointing
`--ours` at the preview tree must RAISE rather than evaluate two builds as one -- the
`n_sites` assertion below is that guard, and it is why the shape check is not just a
sanity check.

⚠ `pangenome_mask.npy` carries 258,610 indices on BOTH builds -- the same COUNT. Every
one of the 258,610 VALUES differs (measured), because the mask is re-indexed onto each
build's own axis. A count-only check would pass on the wrong build; the axis check is on
`pangenome_positions.npy`, whose length differs.

⚠ The arrays are SYMLINKED, not copied: 660 x 1.4 MB of identical bytes buys nothing and
two copies of one array is one chance for them to drift. Both evaluators only read it.
`xsol.npy` is NOT staged -- neither `af_loss.py` nor `ld_loss.py` reads it (that is the
Beagle path's input), so staging it would imply a selection step neither evaluator applies.

⚠ The mirror root holds a COPY of `constants.py`: their `BASE_PATH` is
`dirname(abspath(__file__))` with no env override and `abspath` does not resolve symlinks,
so copying that one file is what redirects `STARTING_DATA_PATH` and `EXPERIMENT_PATH` into
our tree. Nothing under `external/` is written.

⚠ `af_loss.py:156` reads `row["capacity"]` and never uses it, so the CSV must carry that
column or every row raises `KeyError`. Ours writes the literal "n/a" there rather than a
number a later reader could mistake for a real PanMixer setting; theirs writes the real
capacity and "n/a" in the `tau` column.

⚠ The root lives on `/` (550 GB free), NOT on `/data` (91 GB free).

    python scripts/arm1_afld_ourbuild.py
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
OURS = Path("/data/ds85/RanPanMixer/data/panmixer_integration_ourbuild/arm1")
THEIRS = Path("/home/ds85/pm_ourbuild/obf_out")
STARTING = Path("/home/ds85/pm_ourbuild/starting_data")
ROOT = Path("/home/ds85/pm_eval_ourbuild")
CHROM = 21
CEILINGS = ("0p50", "0p55", "0p60", "0p65", "0p70", "0p75", "0p80")
CAPS = ("0p01", "0p05", "0p1", "0p2", "0p3", "0p5", "0p7", "0p9")


def _label_value(label: str) -> float:
    """`'0p80'` -> 0.80. The inverse of `releases.ceiling_dir`'s suffix."""
    return float(label.replace("p", "."))


def _check(arr: np.ndarray, n_sites: int, na: np.ndarray, who: str) -> None:
    """Every pre-flight their evaluators would fail on, asserted here instead.

    ⚠ The `n_sites` clause is the BUILD guard, not a sanity check: a (340824, 2)
    array is the stale preview tree and must never be measured against our
    (44, 340849, 2) `pangenome.npy`.
    """
    if arr.shape != (n_sites, 2):
        raise SystemExit(
            f"{who}: shape {arr.shape}, expected ({n_sites}, 2). An array on a "
            f"DIFFERENT BUILD'S AXIS -- refusing to evaluate two builds as one.")
    if arr.dtype != np.int16:
        raise SystemExit(f"{who}: dtype {arr.dtype}, expected int16")
    bad = int(((arr >= na[:, None]) & (arr >= 0)).sum())
    if bad:
        raise SystemExit(f"{who}: {bad} allele code(s) >= num_alleles[site]")


def _stage(ed: Path, rows: list[dict], n_sites: int, na: np.ndarray) -> None:
    """One experiment: the CSV their `load_data` reads, plus one symlink per row."""
    (ed / "data" / f"chr{CHROM}").mkdir(parents=True, exist_ok=True)
    with (ed / f"data_chr{CHROM}.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subject", "tau", "capacity"])
        for r in rows:
            w.writerow([r["subject"], r["tau"], r["capacity"]])
    for r in rows:
        src = Path(r["source"]).resolve()
        _check(np.load(src), n_sites, na, f"row {r['row_id']} {r['subject']}/{r['point']}")
        rd = ed / "data" / f"chr{CHROM}" / str(r["row_id"])
        rd.mkdir(parents=True, exist_ok=True)
        dst = rd / "new_haplotypes.npy"
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--ours", type=Path, default=OURS)
    ap.add_argument("--theirs", type=Path, default=THEIRS)
    ap.add_argument("--starting-data", type=Path, default=STARTING)
    a = ap.parse_args(argv)
    root, sd = a.root, a.starting_data

    # ---- the mirror. `starting_data` is OURS; `tools` is theirs, unmodified.
    root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PM / "constants.py", root / "constants.py")
    for name, target in (("tools", PM / "tools"),
                         ("downloaded_tools", PM / "downloaded_tools"),
                         ("starting_data", sd)):
        link = root / name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(target.resolve())
    print(f"mirror root     : {root}")
    print(f"  starting_data -> {sd.resolve()}   (OURS)")
    print(f"  tools         -> {(PM / 'tools').resolve()}   (theirs, unmodified)")

    pos = np.load(sd / f"chr{CHROM}" / "pangenome_positions.npy")
    pan = np.load(sd / f"chr{CHROM}" / "pangenome.npy", mmap_mode="r")
    sub = np.load(sd / f"chr{CHROM}" / "pangenome_subjects.npy")
    na = np.load(sd / f"chr{CHROM}" / "num_alleles.npy")
    mask = np.load(sd / f"chr{CHROM}" / "pangenome_mask.npy")
    n_sites = len(pos)
    if pan.shape != (len(sub), n_sites, 2):
        raise SystemExit(f"pangenome.npy {pan.shape} is not ({len(sub)}, {n_sites}, 2)")
    if int(mask.max()) >= n_sites:
        raise SystemExit(f"pangenome_mask max {mask.max()} >= n_sites {n_sites}")
    print(f"  {n_sites:,} sites, {len(sub)} subjects, "
          f"{len(mask):,} SNP-mask indices, pangenome {pan.shape}")

    known = set(sub.tolist())
    ours_t = sorted(d.name for d in a.ours.iterdir() if d.is_dir())
    theirs_t = sorted(d.name for d in a.theirs.iterdir() if d.is_dir())
    # ⚠ their evaluators do `np.where(subjects == subject)[0][0]` and would IndexError
    # on a miss, 300 rows in. And the two tools must cover the SAME individuals or the
    # comparison is not paired -- that is the whole reason this build exists.
    for lbl, ts in (("ours", ours_t), ("theirs", theirs_t)):
        miss = [t for t in ts if t not in known]
        if miss:
            raise SystemExit(f"{lbl}: not in pangenome_subjects.npy: {miss}")
    if ours_t != theirs_t:
        raise SystemExit(
            f"the two tools do not cover the same individuals: "
            f"ours-only {sorted(set(ours_t) - set(theirs_t))}, "
            f"theirs-only {sorted(set(theirs_t) - set(ours_t))}")
    print(f"  {len(ours_t)} individuals, identical set on both tools")

    rows_ours, rows_theirs = [], []
    for t in ours_t:
        for c in CEILINGS:
            src = a.ours / t / f"acc_{c}" / f"chr{CHROM}" / "new_haplotypes.npy"
            if not src.exists():
                raise SystemExit(f"missing our release: {src}")
            rows_ours.append({"row_id": len(rows_ours), "subject": t,
                              "tau": f"acc_{c}", "capacity": "n/a",
                              "point": f"acc_{c}", "knob": "acc_ceiling",
                              "point_value": _label_value(c), "source": str(src)})
    for t in theirs_t:
        for c in CAPS:
            src = a.theirs / t / f"cap{c}" / f"chr{CHROM}" / "new_haplotypes.npy"
            if not src.exists():
                raise SystemExit(f"missing their run: {src}")
            rows_theirs.append({"row_id": len(rows_theirs), "subject": t,
                                "tau": "n/a", "capacity": _label_value(c),
                                "point": f"cap{c}", "knob": "capacity",
                                "point_value": _label_value(c), "source": str(src)})

    _stage(root / "experiments" / "exp_5", rows_ours, n_sites, na)
    print(f"  exp_5 (OURS)  : {len(rows_ours)} rows "
          f"({len(ours_t)} targets x {len(CEILINGS)} acc_ceiling)")
    _stage(root / "experiments" / "exp_6", rows_theirs, n_sites, na)
    print(f"  exp_6 (THEIRS): {len(rows_theirs)} rows "
          f"({len(theirs_t)} subjects x {len(CAPS)} capacity)")

    man = {"root": str(root), "chrom": CHROM, "n_sites": n_sites,
           "starting_data": str(sd.resolve()),
           "build": "OURS -- vg 1.68 GRCh38 cohort44, 340,849 records",
           "experiments": {"5": "ours_acc_ceiling", "6": "theirs_capacity"},
           "points": list(CEILINGS), "caps": list(CAPS),
           "taus": [f"acc_{c}" for c in CEILINGS],
           "targets": ours_t, "rows": rows_ours, "theirs_rows": rows_theirs}
    (root / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(f"\n{len(rows_ours) + len(rows_theirs)} rows x 2 evaluators = "
          f"{2 * (len(rows_ours) + len(rows_theirs))} jobs")
    print(f"manifest -> {root / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
