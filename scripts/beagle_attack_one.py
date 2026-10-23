#!/usr/bin/env python
"""One row of the Beagle reconstruction attack: release array -> VCF -> Beagle -> their stats.

⚠ EVERY SCORING STEP IS THEIR CODE, CALLED NOT COPIED. `produce_vcf_from_np` comes from
`tools/common/convert_2_vcf.py`; `accuracy_computer` comes from `tools/beagle/accuracy_stats.py`.
The only things this script owns are the Beagle invocation (their `beagle_refinement.py` builds a
SLURM array and shells to `sbatch`, which is not installed here, so the java command is lifted
from its `get_command` verbatim) and the file plumbing around it.

⚠ THE OUTPUT FILENAME IS A DOUBLE EXTENSION ON PURPOSE. `accuracy_stats.py:182` reads
`beagle_refined.vcf.gz.vcf.gz`, because Beagle appends `.vcf.gz` to whatever `out=` prefix it is
given and their `out=` already ends in `.vcf.gz`. Passing `out=.../beagle_refined.vcf.gz`
reproduces the name their reader expects.

⚠ `--exclude-self` is a LABELLED DIAGNOSTIC, not their protocol. 39 of our 44 targets are in the
1000G reference panel and their `beagle_refinement.py` passes no `excludesamples=`, so by default
Beagle may copy the target's own true row out of the panel. The default here is THEIR protocol
(no exclusion); `--exclude-self` reruns with the target removed from the panel, and the difference
between the two is how much of the "reconstruction" came from the panel rather than from our
release.

    python scripts/beagle_attack_one.py --row 6 [--exclude-self] [--threads 8]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
PANEL = Path("/data/ds85/RanPanMixer/runs/beagle_attack/panel/"
             "1000g_snps_no_missing_renamed_norm.vcf.gz")
JAR = ROOT / "downloaded_tools" / "beagle.27Feb25.75f.jar"
JRE = Path("/home/ds85/miniconda3/envs/beaglejre/bin")
BCF = Path("/home/ds85/miniconda3/envs/ranpanmixer/bin")
EXPERIMENT = "1"
CHROM = "21"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", type=int, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--xmx", default="32g")
    ap.add_argument("--exclude-self", action="store_true")
    ap.add_argument("--keep-vcf", action="store_true")
    # a second mirror (e.g. the exploratory high-`tau` rows) must be read from ITS
    # root: constants.py there points EXPERIMENT_PATH at that tree
    ap.add_argument("--root", type=Path, default=ROOT)
    # ⚠ exp_6 holds PANMIXER'S OWN 8 capacity releases, so the same Beagle attack can be
    # run on their output as on ours. Their rows are NOT in the main manifest (whose
    # row_ids 0..307 are ours), so `--subject` overrides the lookup; without it row 0 of
    # exp_6 would silently be scored against OUR row 0's subject.
    ap.add_argument("--experiment", default=EXPERIMENT)
    ap.add_argument("--subject", default=None,
                    help="override the manifest lookup (required for exp_6)")
    a = ap.parse_args(argv)
    exp = a.experiment
    root = a.root
    jar = root / "downloaded_tools" / "beagle.27Feb25.75f.jar"

    sys.path.insert(0, str(root))                 # the mirror: constants.py + tools
    import numpy as np
    from tools.common.convert_2_vcf import produce_vcf_from_np
    from tools.beagle.accuracy_stats import accuracy_computer

    if a.subject is None:
        man = json.loads((root / "manifest.json").read_text())
        row = next(r for r in man["rows"] if r["row_id"] == a.row)
        subject, tau = row["subject"], row["tau"]
    else:
        subject, tau = a.subject, f"exp{exp}"
    rd = root / "experiments" / f"exp_{exp}" / "data" / f"chr{CHROM}" / str(a.row)
    tag = f"row {a.row} {subject} tau={tau}" + (" [exclude-self]" if a.exclude_self else "")
    print(f"=== {tag} ===", flush=True)

    env = dict(os.environ)
    env["PATH"] = f"{BCF}:{env['PATH']}"
    env["PYTHONDONTWRITEBYTECODE"] = "1"

    # ---- 1. their converter: release array -> VCF on the pangenome template
    vcf = rd / "new_haplotypes.vcf"
    t0 = time.time()
    if not (rd / "new_haplotypes.vcf.gz").exists():
        arr = np.load(rd / "new_haplotypes.npy")
        produce_vcf_from_np(
            arr, str(root / "starting_data" / f"chr{CHROM}" / "pangenome.vcf.gz"),
            subject, str(vcf))
        subprocess.run([str(BCF / "bgzip"), "-@2", "-f", str(vcf)], check=True, env=env)
        subprocess.run([str(BCF / "bcftools"), "index", "-f", f"{vcf}.gz"],
                       check=True, env=env)
    print(f"  VCF ready ({time.time()-t0:.0f}s)", flush=True)

    # ---- 2. Beagle, the command lifted from their beagle_refinement.get_command
    out_prefix = rd / ("beagle_selfout.vcf.gz" if a.exclude_self
                       else "beagle_refined.vcf.gz")
    # ⚠⚠ THE PAPER'S PROTOCOL NEEDS A PER-TARGET PANEL, NOT `excludesamples=`. Measured
    # 2026-10-14: `excludesamples=` drops the named sample from BOTH the reference and the
    # STUDY set (3,202 -> 3,201 reference AND 44 -> 43 study samples on row 7), so the target
    # is absent from Beagle's output entirely and `accuracy_stats.py` raises at
    # `list(beagle_subjects).index(subject)`. That raise is the only reason this was caught;
    # a lookup that silently returned a different sample would have produced a complete set
    # of plausible numbers for the wrong individual.
    ref = PANEL
    if a.exclude_self:
        ref = (PANEL.parent / "per_target" / f"panel_minus_{subject}.vcf.gz")
        if not ref.exists():
            print(f"  MISSING per-target panel {ref}")
            return 1
    cmd = [str(JRE / "java"), f"-Xmx{a.xmx}", "-jar", str(jar),
           f"gt={vcf}.gz", f"ref={ref}", f"out={out_prefix}",
           f"nthreads={a.threads}"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    (rd / ("beagle_selfout.log" if a.exclude_self else "beagle_run.log")).write_text(
        r.stdout + "\n==== stderr ====\n" + r.stderr)
    if r.returncode != 0:
        print(f"  BEAGLE FAILED rc={r.returncode}\n{r.stderr[-2000:]}")
        return 1
    print(f"  Beagle done ({time.time()-t0:.0f}s)", flush=True)

    # ---- 3. THEIR accuracy_stats, unmodified
    # ⚠⚠ BOTH of their hardcoded names have to be parked in exclude-self mode, and missing
    # the SECOND one would have silently destroyed the protocol column. Their
    # `accuracy_stats.py:182` reads `beagle_refined.vcf.gz.vcf.gz` and their
    # `accuracy_stats.py:196` WRITES `accuracy_results.json`, both unconditionally. Swapping
    # only the VCF would leave their writer overwriting the protocol run's results with the
    # exclude-self numbers; this script would then read that file back and save it a second
    # time as `..._selfout.json`, so every row would end up with two files holding the SAME
    # numbers and the two protocols would appear to agree exactly. That is a defect that
    # presents as a finding. Caught 2026-10-14 by re-reading this function before launching
    # the second pass, not by the numbers looking wrong -- they would not have.
    vcf_live = rd / "beagle_refined.vcf.gz.vcf.gz"
    vcf_stash = rd / "beagle_refined.PROTOCOL.vcf.gz.vcf.gz"
    res_live = rd / "accuracy_results.json"
    res_stash = rd / "accuracy_results.PROTOCOL.json"
    if a.exclude_self:
        if vcf_live.exists():
            vcf_live.rename(vcf_stash)
        (rd / "beagle_selfout.vcf.gz.vcf.gz").rename(vcf_live)
        if res_live.exists():
            res_live.rename(res_stash)
    t0 = time.time()
    cwd = Path.cwd()
    os.chdir(rd)                      # their code savefig/writes relative in places
    try:
        accuracy_computer(exp, CHROM, a.row)
        res = json.loads(res_live.read_text())
    finally:
        os.chdir(cwd)
        if a.exclude_self:
            vcf_live.rename(rd / "beagle_selfout.vcf.gz.vcf.gz")
            if vcf_stash.exists():
                vcf_stash.rename(vcf_live)
            # their writer has just clobbered accuracy_results.json; put the protocol
            # numbers back and keep the exclude-self ones under their own name.
            res_live.unlink(missing_ok=True)
            if res_stash.exists():
                res_stash.rename(res_live)
    name = "accuracy_results_selfout.json" if a.exclude_self else "accuracy_results.json"
    (rd / name).write_text(json.dumps(res, indent=1) + "\n")
    if a.exclude_self:
        # ⚠ assert the protocol column SURVIVED this row, every row, rather than trusting
        # the rename dance. A silent loss here is unrecoverable without a 2.5-hour re-run.
        prot = rd / "accuracy_results_protocol.json"
        assert prot.exists(), f"row {a.row}: the protocol backup is gone"
        # the restored live file must be the PROTOCOL numbers, byte for byte. This is the
        # check the whole stash dance exists for, and it CAN fail: drop the res_stash
        # rename above and it fires immediately, because res_live then holds the
        # exclude-self numbers that their writer just put there.
        p = json.loads(prot.read_text())
        live = json.loads(res_live.read_text())
        assert p == live, (
            f"row {a.row}: accuracy_results.json no longer holds the protocol numbers -- "
            f"their writer clobbered it and the restore did not take. "
            f"protocol snp_only reconstruction_accuracy="
            f"{p['snp_only']['reconstruction_accuracy']}, live="
            f"{live['snp_only']['reconstruction_accuracy']}")
        print(f"  protocol column intact: snp_only reconstruction_accuracy "
              f"{p['snp_only']['reconstruction_accuracy']:.6f} (protocol) vs "
              f"{res['snp_only']['reconstruction_accuracy']:.6f} (exclude-self)")
    print(f"  accuracy_stats done ({time.time()-t0:.0f}s)")
    print(f"  reconstruction_accuracy[all] = {res['all']['reconstruction_accuracy']}")

    if not a.keep_vcf:
        for f in (f"{vcf}.gz", f"{vcf}.gz.csi"):
            Path(f).unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
