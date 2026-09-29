"""Steps 1-2: graph -> VCF -> the 44-sample cohort VCF.

Thin subprocess wrappers. These are the slow steps (deconstruct is 1-4 h on
chr21) and the ones with no numeric content, so they stay separate from the
converter that will be rewritten repeatedly.
"""
from __future__ import annotations
from pathlib import Path


def deconstruct(gfa: Path, out_vcf: Path, ref_prefix: str = "grch38", threads: int = 48) -> Path:
    """`vg deconstruct` the GFA against the reference backbone path.

    ⚠ UNVERIFIED, and the reason the chrY smoke test exists: sample grouping.
    Paths in the GFA are PanSN `SAMPLE#HAP#CONTIG` and each assembly is split
    across many contigs, so producing 45 sample columns rather than hundreds of
    contig columns depends on vg understanding `#` as the separator. The `-H`
    flag that declares it is reported as deprecated in vg 1.68. If grouping
    silently changes, every downstream array is wrong -- obviously so, but only
    after hours of compute. Run chrY first.
    """
    raise NotImplementedError


def drop_chm13(in_vcf: Path, out_vcf: Path) -> Path:
    """`bcftools view -s ^chm13` -> 44 samples = 88 haplotypes.

    CHM13 is in the graph (measured: 46 path prefixes = 44 donors + chm13 +
    grch38). It is a reference assembly from a homozygous cell line, not a
    cohort donor whose privacy is at stake, so it is dropped rather than
    modelled. Note chrX removal is a no-op for a per-chromosome build.
    """
    raise NotImplementedError
