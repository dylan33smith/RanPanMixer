"""Steps 1-2: graph -> VCF -> the 44-sample cohort VCF.

Thin subprocess wrappers. These are the slow steps -- deconstruct is the long
pole -- and they have no numeric content, so they stay separate from the
converter that gets rewritten repeatedly.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import store

VG = Path("/data/ds85/RanPanMixer/external/PanMixer/downloaded_tools/vg")


def _run(cmd: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path:
        with open(stdout_path, "wb") as out:
            subprocess.run(cmd, stdout=out, check=True)
    else:
        subprocess.run(cmd, check=True)


def deconstruct(gfa: Path, out_vcf: Path, *, ref_prefix: str = "grch38",
                threads: int = 32, vg: Path = VG) -> dict:
    """`vg deconstruct` the GFA against the reference backbone path.

    Sample grouping was the risk here and it is settled: paths are PanSN
    `SAMPLE#HAP#CONTIG` and each assembly is split across many contigs, so
    producing sample columns rather than contig columns depends on vg
    understanding `#`. The `-H` flag that used to declare it is deprecated in
    vg 1.68 -- but a chrY smoke test (2:01, 34,014 records) produced 17 clean
    sample columns, so grouping is auto-detected and correct.
    """
    out_vcf = Path(out_vcf); out_vcf.parent.mkdir(parents=True, exist_ok=True)
    _run([str(vg), "deconstruct", "-P", ref_prefix, "-a", "-t", str(threads), str(gfa)],
         stdout_path=out_vcf)
    return {"vcf": str(out_vcf), "sha256": store.sha256_file(out_vcf)}


def drop_chm13(in_vcf: Path, out_vcf: Path, *, drop: str = "chm13") -> dict:
    """`bcftools view -s ^chm13` -> 44 samples = 88 haplotypes.

    CHM13 is in the graph -- measured, 46 path prefixes = 44 donors + chm13 +
    grch38. It is a reference assembly from a homozygous cell line, not a cohort
    donor whose privacy is at stake, so it is dropped rather than modelled.
    chrX removal is a no-op for a per-chromosome build.
    """
    out_vcf = Path(out_vcf); out_vcf.parent.mkdir(parents=True, exist_ok=True)
    _run(["bcftools", "view", "-s", f"^{drop}", "--force-samples",
          "-Oz", "-o", str(out_vcf), str(in_vcf)])
    _run(["bcftools", "index", "-f", str(out_vcf)])
    n = subprocess.run(["bcftools", "query", "-l", str(out_vcf)],
                       capture_output=True, text=True, check=True).stdout.split()
    return {"vcf": str(out_vcf), "n_samples": len(n),
            "sha256": store.sha256_file(out_vcf), "dropped": drop}


def slice_region(in_vcf: Path, out_vcf: Path, region: str) -> dict:
    """Cut the development fixture out of the cohort VCF.

    The default region is chr21:14-16 Mb: 16,622 records including the 257 kb
    hypervariable bubble and the 90-allele site, i.e. the worst case on the
    chromosome, and entirely above the chain_span start so every record is live.
    """
    out_vcf = Path(out_vcf); out_vcf.parent.mkdir(parents=True, exist_ok=True)
    _run(["bcftools", "view", "-r", region, "-Oz", "-o", str(out_vcf), str(in_vcf)])
    _run(["bcftools", "index", "-f", str(out_vcf)])
    return {"vcf": str(out_vcf), "region": region, "sha256": store.sha256_file(out_vcf)}
