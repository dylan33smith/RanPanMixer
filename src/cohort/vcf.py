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
        decompressed_here = True
    else:
        decompressed_here = False

    _run([str(resolved), "deconstruct", "-P", ref_prefix, "-a", "-t", str(threads), str(gfa)],
         stdout_path=out_vcf)

    # ⚠⚠ THE GFA DIGEST IS TAKEN BEFORE THE FILE IS DELETED, AND THE ORDER IS THE
    # WHOLE POINT. An earlier version of this function unlinked the decompressed
    # copy first and then passed the same path to `sha256_file` three lines later,
    # which raises FileNotFoundError AFTER vg has already done ~6 minutes of wall work
    # and about 2 CPU-hours across 32 threads.
    # Six chromosomes failed that way on 2026-10-10 (chr10-14, chr17) and the
    # traceback named `store.sha256_file`, not the unlink, so the cause read as a
    # provenance bug rather than an ordering one. Digest first, delete second.
    gfa_digest = store.sha256_file(gfa)

    # ⚠ DELETE THE DECOMPRESSED COPY. It is 4-9 GB per chromosome and is
    # regenerable from the .gz in seconds, so caching it costs about 86 GB across
    # 22 autosomes to save nothing that matters. The earlier code kept it
    # deliberately ("if not plain.exists()") as a resume convenience for one
    # chromosome, which does not scale to a genome: measured 2026-10-10, six
    # finished chromosomes were holding 25 GB of graphs nothing would read again.
    # ⚠ Only removed if WE created it, so a hand-decompressed graph a caller staged
    # on purpose is left alone.
    if decompressed_here and gfa.exists():
        gfa.unlink()
        print(f"    removed the decompressed {gfa.name} "
              f"(regenerable from the .gz; keeping it costs GB per chromosome)")
    out = {"vcf": str(out_vcf), "gfa_used": str(gfa),
           "sha256": store.sha256_file(out_vcf), "vg": vginfo}
    # close the lineage: GFA + vg -> deconstruct.vcf had no sidecar at all
    store.attach_prov(out_vcf, inputs={"gfa": gfa_digest, "vg": vginfo["sha256"]},
                      params={"ref_prefix": ref_prefix, "threads": threads,
                              "vg_path": vginfo["path"], "vg_version": vginfo["version"]},
                      note="vg deconstruct output. The vg VERSION is load-bearing: it "
                           "determines ALT allele ordering and so every allele index.")
    return out


def attach_deconstruct_prov(gfa: Path, out_vcf: Path, *, ref_prefix: str = "grch38",
                            threads: int = 32, vg: Path | None = None,
                            cfg: dict | None = None) -> dict:
    """Give an ALREADY-WRITTEN deconstruct VCF the provenance sidecar it is missing.

    ⚠ WHY A SEPARATE ENTRY POINT AND NOT A FLAG ON `deconstruct`. vg writes the VCF
    and the sidecar is attached afterwards, so a crash between the two leaves a
    COMPLETE, CORRECT VCF with no provenance -- and this project treats an artifact
    without provenance as not quotable, so the VCF is unusable while being perfectly
    good. Re-running vg to recover a 1 KB JSON file costs about 6 minutes of wall time
    and 2 CPU-hours per chromosome. This recovers it in 9 s, measured on chr10.

    ⚠ IT REFUSES TO OVERWRITE an existing sidecar. A sidecar that is already there
    was written by the run that produced the VCF, and silently replacing it with one
    assembled afterwards would turn real provenance into a reconstruction while
    looking identical.
    """
    out_vcf = Path(out_vcf)
    gfa = Path(gfa)
    prov = out_vcf.with_name(out_vcf.name + store.PROV_SUFFIX)
    if not out_vcf.exists() or out_vcf.stat().st_size == 0:
        raise FileNotFoundError(
            f"{out_vcf} does not exist or is empty, so there is nothing to attach "
            f"provenance to. Run deconstruct instead.")
    if prov.exists():
        raise FileExistsError(
            f"{prov} already exists and was written by the run that produced the "
            f"VCF. Refusing to replace real provenance with a reconstruction.")
    vginfo = verify_vg(resolve_vg(cfg, vg), cfg)
    # ⚠ The digest must be of the DECOMPRESSED bytes, to match the field the
    # chromosomes built before the crash carry. Streamed from the .gz, never staged.
    gfa_digest = (store.sha256_gz_plain(gfa) if gfa.suffix == ".gz"
                  else store.sha256_file(gfa))
    store.attach_prov(out_vcf, inputs={"gfa": gfa_digest, "vg": vginfo["sha256"]},
                      params={"ref_prefix": ref_prefix, "threads": threads,
                              "vg_path": vginfo["path"], "vg_version": vginfo["version"],
                              "prov_recovered_after_crash": True},
                      note="vg deconstruct output. The vg VERSION is load-bearing: it "
                           "determines ALT allele ordering and so every allele index. "
                           "⚠ This sidecar was attached AFTER the fact, by "
                           "attach_deconstruct_prov, because the run that wrote the VCF "
                           "crashed between writing it and recording it. The VCF bytes "
                           "are vg's own output and unmodified; the GFA digest was "
                           "recomputed from the .gz, whose decompression is "
                           "deterministic.")
    return {"vcf": str(out_vcf), "sha256": store.sha256_file(out_vcf),
            "gfa": gfa_digest, "vg": vginfo, "recovered": True}


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
