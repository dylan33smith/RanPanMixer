"""Step 3: the cohort VCF -> the arrays. THE IRREVERSIBLE STEP.

This is the only preprocessing step with real design content, and the one that
will be rewritten most, which is why it is alone in its own module.

Irreversible because the VCF line is the only place four things coexist: the
genotype, the INFO tags (LV, PS, CONFLICT), the REF/ALT strings, and the record's
position in the file. A converter that keeps only position and genotype -- which
is what the upstream one does -- destroys the rest permanently. Cheap now,
unrecoverable later.

TWO PASSES, because "inherited" is not a property of the record it applies to:
  pass 1  matrix, per-record codes (not-applicable / conflict), lengths, sites,
          and a PS -> row index
  pass 2  fill code 2 where the PARENT record's call for that sample is missing

PRECEDENCE, pinned in docs/terms.md under `reason_array`:
  conflict (3) > not applicable (1) > inherited (2) > uncategorised (4)
Conflict wins because the CONFLICT tag is the caller stating outright that it
could not choose, which beats a structural inference.

ASSEMBLY GAPS GET NO CODE. They are a property of a RUN, not of a cell, and
classifying them per cell needs a length threshold that has never been
validated -- it would claim 92.06% of all missing cells at >=10 and 68.67% at
>=1000. They go to `missing_runs.tsv` as intervals so any threshold can be
applied, changed or swept downstream without re-reading the VCF.
"""
from __future__ import annotations
from pathlib import Path

CALLED         = 0
NOT_APPLICABLE = 1   # nested child whose parent allele this haplotype lacks
INHERITED      = 2   # the parent record's own call was missing
CONFLICT       = 3   # the CONFLICT tag names this sample
UNCATEGORISED  = 4   # a no-call with no structural explanation (4.57% on chr21)


def build(cohort_vcf: Path, outdir: Path) -> dict:
    """Emit haplotypes, reason, missing_runs, allele_lengths, sites, ids.

    Returns a summary dict for the manifest. Takes NO target argument, and must
    not acquire one -- see cohort/__init__.py.
    """
    raise NotImplementedError
