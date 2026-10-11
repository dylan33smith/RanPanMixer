#!/usr/bin/env python
"""Chunk 6(c): map all 10 targets onto all 22 autosomes' site axes.

⚠ WHY A DRIVER AND NOT A LOOP. Two things here are easy to get wrong silently and
both have bitten this project already:

  * **The per-chromosome target VCF has no producer on record.** chr21's
    `target.chr21.vcf.gz` carries NO provenance sidecar — it was made by hand, so the
    one target set we have is not reproducible from anything written down. This script
    is that producer, and it records the extraction in the outputs' provenance.
  * **The extracted VCF is deleted after use, and the DIGEST IS TAKEN FIRST.** It is
    about 6.8 GB per target across 22 autosomes, 68 GB for ten, and it is a pure
    `bcftools view -s` of a panel we already keep checksummed. ⚠ Taking the digest
    before the unlink is not an incidental detail: doing it the other way round is
    exactly the bug that cost nine chromosomes on 2026-10-10, where a cleanup was
    inserted above the line that hashed the file.

⚠ AND IT ASSERTS ITS OWN COMPLETENESS over every (target, chromosome) pair and every
expected artifact, because a filter cannot tell "found nothing" from "ran nothing".
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import runner as R                                        # noqa: E402
import store                                              # noqa: E402

PANEL = Path("/data/ds85/RanPanMixer/data/raw")
EXPECTED = ("path.npy", "readable.npy", "target_reason.npy", "representability.json")


def panel_for(chrom: str) -> Path:
    """chr21's panel predates the bulk fetch and is on disk under another name."""
    return (PANEL / "1000g_30x_phased.vcf.gz" if chrom == "chr21"
            else PANEL / f"1000g_30x_{chrom}.vcf.gz")


def _one(job: tuple) -> dict:
    target, chrom = job
    out = REPO / "data" / "targets" / target / chrom
    out.mkdir(parents=True, exist_ok=True)
    tv = out / f"target.{chrom}.vcf.gz"
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    t0 = time.time()
    log = REPO / "logs" / "targets" / f"{target}.{chrom}.log"
    log.parent.mkdir(parents=True, exist_ok=True)

    # ⚠ RESUME. 220 jobs is enough that one failure must not cost the other 219, and
    # chr21 is already prepared for all ten targets by hand — re-running it would
    # overwrite artifacts that recorded numbers were measured against. A pair is done
    # when every expected artifact exists and is non-empty.
    if all((out / f).exists() and (out / f).stat().st_size > 0 for f in EXPECTED):
        return {"target": target, "chrom": chrom, "seconds": 0.0,
                "freed_bytes": 0, "extracted": False, "skipped": True}

    made_here = False
    with open(log, "w") as fh:
        # ⚠⚠ A LEFTOVER VCF IS EVIDENCE OF AN INTERRUPTED RUN, NOT OF PROGRESS.
        # This function DELETES the extracted VCF on success, so finding one here with
        # no `path.npy` beside it means a previous attempt died mid-write -- and a
        # half-written bgzf file is NON-EMPTY, so an `st_size == 0` test calls it
        # complete and feeds a truncated VCF to prep-target. The target path would then
        # be silently missing its tail, with no error anywhere.
        # ⚠ MEASURED: interrupting this driver on 2026-10-11 left nine such files,
        # one of 581,931,908 B, every one of which the previous check would have reused.
        # ⚠ chr21 is the exception and must be left alone: its VCFs are the hand-made
        # originals that every recorded chr21 number was measured against, and they
        # legitimately sit beside a completed `path.npy`.
        if tv.exists() and chrom != "chr21":
            fh.write(f"discarding a leftover {tv.name} ({tv.stat().st_size} B) -- no "
                     f"path.npy beside it, so a previous attempt died mid-write\n")
            tv.unlink()
            for side in (Path(str(tv) + ".csi"), Path(str(tv) + ".tbi")):
                if side.exists():
                    side.unlink()

        if not tv.exists() or tv.stat().st_size == 0:
            src = panel_for(chrom)
            if not src.exists():
                raise FileNotFoundError(f"{src} is not on disk; cannot extract {target}")
            fh.write(f"extracting {target} from {src.name}\n"); fh.flush()
            p = subprocess.run(["bcftools", "view", "-s", target, "-Oz", "-o", str(tv),
                                str(src)], stdout=fh, stderr=subprocess.STDOUT)
            if p.returncode != 0:
                raise RuntimeError(f"{target}/{chrom}: bcftools view exited {p.returncode}")
            p = subprocess.run(["bcftools", "index", "-f", str(tv)],
                               stdout=fh, stderr=subprocess.STDOUT)
            if p.returncode != 0:
                raise RuntimeError(f"{target}/{chrom}: bcftools index exited {p.returncode}")
            made_here = True

        # ⚠ THE DIGEST BEFORE ANY CLEANUP, for the reason in the module docstring.
        vcf_digest = store.sha256_file(tv)
        fh.write(f"vcf sha256 {vcf_digest}\n"); fh.flush()

        p = subprocess.run([sys.executable, "-c",
                            "import cli, sys; sys.exit(cli.target_main())",
                            "--target", target, "--chrom", chrom,
                            "--target-vcf", str(tv)],
                           stdout=fh, stderr=subprocess.STDOUT, cwd=REPO, env=env)
        if p.returncode != 0:
            raise RuntimeError(f"{target}/{chrom}: prep-target exited {p.returncode}; "
                               f"see {log}")

    missing = [f for f in EXPECTED
               if not (out / f).exists() or (out / f).stat().st_size == 0]
    if missing:
        raise RuntimeError(f"{target}/{chrom}: prep-target exited 0 but "
                           f"{', '.join(missing)} is missing or empty")

    # record the extraction on an artifact that SURVIVES, then reclaim the VCF
    store.attach_prov(out / "representability.json",
                      inputs={"panel": store.sha256_file(panel_for(chrom)),
                              "target_vcf": vcf_digest},
                      params={"target": target, "chrom": chrom,
                              "extracted_by": f"bcftools view -s {target} "
                                              f"{panel_for(chrom).name}",
                              "target_vcf_deleted": made_here},
                      note="⚠ the single-sample target VCF this was built from is "
                           "DELETED after use (about 6.8 GB per target genome-wide) "
                           "and is regenerable by the recorded bcftools command from "
                           "the panel digest above. The digest was taken BEFORE the "
                           "deletion.")
    freed = 0
    if made_here and chrom != "chr21":        # chr21's was hand-made; leave it alone
        for f in (tv, Path(str(tv) + ".csi"), Path(str(tv) + ".tbi")):
            if f.exists():
                freed += f.stat().st_size
                f.unlink()
    return {"target": target, "chrom": chrom, "seconds": round(time.time() - t0, 1),
            "freed_bytes": freed, "extracted": made_here}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--chroms", type=int, nargs="+", default=list(range(1, 23)))
    ap.add_argument("--dry-run", action="store_true")
    # ⚠ A CAP, because this is not the only thing on the machine. `plan_workers` sizes
    # the pool from free RAM and idle cores, which is right when a script owns the box
    # and wrong when another job already has 46 workers on it: both then run at half
    # speed and finish later than if they had been serialised. Cap deliberately when
    # sharing, and say what the cap was.
    ap.add_argument("--workers", type=int, default=None,
                    help="cap the pool (use when another job holds the cores)")
    a = ap.parse_args(argv)

    targets = a.targets or sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                                  if p.is_dir())
    jobs = [(t, f"chr{c}") for t in targets for c in a.chroms]
    print(f"{len(targets)} target(s) x {len(a.chroms)} chromosome(s) = {len(jobs)} jobs")
    if a.dry_run:
        for t, c in jobs[:5]:
            print(f"  {t}/{c}  panel={panel_for(c).name}  "
                  f"exists={(REPO / 'data' / 'targets' / t / c / 'path.npy').exists()}")
        return 0

    # the cohort arrays for one chromosome plus a parsed single-sample VCF
    plan = R.plan_workers(2 << 30, n_jobs=len(jobs))
    workers = min(plan.workers, a.workers) if a.workers else plan.workers
    print(f"{plan}")
    if a.workers:
        print(f"⚠ CAPPED to {workers} workers by --workers {a.workers} "
              f"(another job holds the cores)")
    print()
    res = R.run_jobs(_one, jobs, workers=workers, label="target-chromosomes",
                     progress_every=max(1, len(jobs) // 20))

    freed = sum(r["freed_bytes"] for r in res)
    skipped = sum(1 for r in res if r.get("skipped"))
    print(f"\n{len(res)} pair(s): {len(res) - skipped} built, {skipped} already done; "
          f"{freed / 1e9:.1f} GB of extracted VCFs reclaimed")

    # ⚠ completeness over the FULL grid, asserted independently of what ran
    bad = []
    for t in targets:
        for c in a.chroms:
            d = REPO / "data" / "targets" / t / f"chr{c}"
            for f in EXPECTED:
                if not (d / f).exists() or (d / f).stat().st_size == 0:
                    bad.append(f"{t}/chr{c}: missing or empty {f}")
    if bad:
        print(f"\nINCOMPLETE — {len(bad)} problem(s):")
        for b in bad[:20]:
            print(f"  {b}")
        print("\nDo NOT quote a genome-wide target number from this tree.")
        return 1
    print(f"\nall {len(targets)} target(s) complete on all {len(a.chroms)} chromosome(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
