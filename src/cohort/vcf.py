"""Steps 1-2: graph -> VCF -> the 44-sample cohort VCF.

Thin subprocess wrappers. These are the slow steps -- deconstruct is the long
pole -- and they have no numeric content, so they stay separate from the
converter that gets rewritten repeatedly.
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import store

# vg is a third-party binary (v1.68.0). It lives under OUR data tree, not inside
# PanMixer's checkout, so nothing in this pipeline reads a path belonging to their
# tool -- see docs/plan.md on independence. Overridable by env var for portability.
VG = Path(os.environ.get("RANPANMIXER_VG", "data/tools/vg"))


def _run(cmd: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path:
        with open(stdout_path, "wb") as out:
            subprocess.run(cmd, stdout=out, check=True)
    else:
        subprocess.run(cmd, check=True)


def resolve_vg(cfg: dict | None = None, vg: Path | None = None) -> Path:
    """Which vg binary will actually run. Precedence: explicit arg, $RANPANMIXER_VG,
    the manifest's tools.vg.path, then the default.

    ⚠ The manifest declared a vg path and sha256 that NO code read -- the pin was
    decoration, and an impostor on $RANPANMIXER_VG ran while the digest test passed.
    vg's version determines ALT allele ordering and therefore every allele index
    downstream, so this is the one pin that must not be theatre.
    """
    if vg is not None:
        return Path(vg)
    env = os.environ.get("RANPANMIXER_VG")
    if env:
        return Path(env)
    if cfg:
        pth = (cfg.get("tools", {}).get("vg", {}) or {}).get("path")
        if pth:
            return Path(pth)
    return VG


def verify_vg(resolved: Path, cfg: dict | None) -> dict:
    """Check the binary that will run against the manifest's declared digest."""
    info = {"path": str(resolved), "sha256": store.sha256_file(resolved)}
    try:
        info["version"] = subprocess.run([str(resolved), "version"], capture_output=True,
                                         text=True, check=True).stdout.splitlines()[0]
    except Exception:
        info["version"] = "unknown"
    declared = ((cfg or {}).get("tools", {}).get("vg", {}) or {}).get("sha256")
    if declared and declared != info["sha256"]:
        raise ValueError(
            f"vg digest mismatch. The manifest declares {declared} but the binary that "
            f"would run, {resolved}, hashes to {info['sha256']}. vg's version decides "
            f"ALT allele ordering and therefore every allele index downstream -- refuse."
        )
    info["digest_checked_against_manifest"] = bool(declared)
    return info


def deconstruct(gfa: Path, out_vcf: Path, *, ref_prefix: str = "grch38",
                threads: int = 32, vg: Path | None = None, cfg: dict | None = None) -> dict:
    """`vg deconstruct` the GFA against the reference backbone path.

    Sample grouping was the risk here and it is settled: paths are PanSN
    `SAMPLE#HAP#CONTIG` and each assembly is split across many contigs, so
    producing sample columns rather than contig columns depends on vg
    understanding `#`. The `-H` flag that used to declare it is deprecated in
    vg 1.68 -- but a chrY smoke test (2:01, 34,014 records) produced 17 clean
    sample columns, so grouping is auto-detected and correct.
    """
    gfa = Path(gfa)
    resolved = resolve_vg(cfg, vg)
    vginfo = verify_vg(resolved, cfg)
    out_vcf = Path(out_vcf); out_vcf.parent.mkdir(parents=True, exist_ok=True)

    # vg cannot read a gzipped GFA -- it fails with "invalid Graph message", which
    # is not obviously a compression error. Decompress beside the output first.
    # (This path was initially untested because the first chr21 run decompressed by
    # hand; the CLI would have failed on the manifest's .gz path.)
    if gfa.suffix == ".gz":
        plain = out_vcf.parent / gfa.with_suffix("").name
        if not plain.exists() or plain.stat().st_size == 0:
            print(f"    decompressing {gfa.name} -> {plain.name}")
            with gzip.open(gfa, "rb") as fi, open(plain, "wb") as fo:
                shutil.copyfileobj(fi, fo, length=1 << 24)
        gfa = plain

    _run([str(resolved), "deconstruct", "-P", ref_prefix, "-a", "-t", str(threads), str(gfa)],
         stdout_path=out_vcf)
    out = {"vcf": str(out_vcf), "gfa_used": str(gfa),
           "sha256": store.sha256_file(out_vcf), "vg": vginfo}
    # close the lineage: GFA + vg -> deconstruct.vcf had no sidecar at all
    store.attach_prov(out_vcf, inputs={"gfa": store.sha256_file(gfa), "vg": vginfo["sha256"]},
                      params={"ref_prefix": ref_prefix, "threads": threads,
                              "vg_path": vginfo["path"], "vg_version": vginfo["version"]},
                      note="vg deconstruct output. The vg VERSION is load-bearing: it "
                           "determines ALT allele ordering and so every allele index.")
    return out


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
    store.attach_prov(out_vcf, inputs={"deconstruct_vcf": store.sha256_file(in_vcf)},
                      params={"dropped_sample": drop, "n_samples": len(n)},
                      note="the 44-sample cohort VCF; chm13 removed")
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
