"""Step 7: a target VCF -> its path over the site axis, plus the audit.

The join is on (POS, REF, ALT) against sites.tsv. No alignment and no
`vg giraffe`: the 1000G panel is phased and GRCh38-called, so a sample column
already IS a path once joined -- which also sidesteps HPRC marking short-read
mapping "untested" for the PGGB graph.
"""
from __future__ import annotations
from pathlib import Path


def build(target_vcf: Path, sites_tsv: Path, outdir: Path) -> dict:
    """Emit path.npy (n_sites, 2) and representability.json.

    ⚠ representability is an AUDIT ARTIFACT ONLY. About 9.9% of a real target's
    chr21 non-reference calls have no record in the graph and carry 20.5% of its
    -log f information. The fix is to MEASURE that loss and report it against the
    target_fidelity ceiling -- never to add a record for a target-specific
    variant, which would make output_support depend on the target and void
    Theorem 1.
    """
    raise NotImplementedError
