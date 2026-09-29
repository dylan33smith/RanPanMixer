"""Step 7: a target VCF -> its path over the site axis, plus the audit.

The join is strict on (POS, REF, full ALT string). Strictness is not fussiness:
it is what makes the allele INDICES transferable. If the ALT strings match
exactly then target allele 2 and cohort allele 2 are the same piece of DNA. Under
a relaxed (POS, REF) join they need not be, and at a record declaring 90 alleles
"the first match" is arbitrary -- which is a silent way to write the wrong allele
into the private input.

No alignment step and no `vg giraffe`: the 1000G panel is phased and
GRCh38-called, so a sample column already IS a path once joined.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np

import sites as sites_mod
import store

# What a site means when the target VCF has no record for it at all.
#   "missing"    -> -1. Honest: the panel simply does not cover that record, so
#                   the target's allele there is unknown. Costs tilt, because a
#                   position with no target allele contributes no phi_v.
#   "reference"  -> 0.  Defensible on the reading that a phased panel lists every
#                   site polymorphic in it, so absence means "nobody varies here".
#                   Cheaper on budget, but it INVENTS a call.
ABSENT_POLICIES = ("missing", "reference")


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _parse_gt(field: str) -> tuple[int, int]:
    gt = field.split(":", 1)[0]
    if "|" in gt:
        a, b = gt.split("|", 1)
    elif "/" in gt:
        a, b = gt.split("/", 1)
    else:
        a, b = gt, None
    x = -1 if a in (".", "") else int(a)
    y = -1 if b is None or b in (".", "") else int(b)
    return x, y


def build(target_vcf: Path, cohort_dir: Path, outdir: Path, *,
          sample: str, absent_policy: str = "missing") -> dict:
    """Emit path.npy (n_sites, 2) and representability.json.

    ⚠ representability is an AUDIT ARTIFACT ONLY. About 9.9% of a real target's
    chr21 non-reference calls have no record in the graph and carry 20.5% of its
    -log f information. The fix is to MEASURE that loss and report it against the
    target_fidelity ceiling -- never to add a record for a target-specific
    variant, which would make output_support depend on the target and void
    Theorem 1.
    """
    if absent_policy not in ABSENT_POLICIES:
        raise ValueError(f"absent_policy must be one of {ABSENT_POLICIES}")
    target_vcf, cohort_dir, outdir = Path(target_vcf), Path(cohort_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    s = sites_mod.load_sites(cohort_dir / "sites.tsv")
    strict: dict[tuple[int, str, str], int] = {}
    pos_seen: dict[int, int] = {}
    for i in range(s.n_sites):
        strict.setdefault((int(s.pos[i]), s.ref[i], s.alt[i]), i)
        pos_seen.setdefault(int(s.pos[i]), i)

    # The cohort's coordinate range. Representability must be reported RELATIVE
    # to the region actually modelled: on a 2 Mb development slice, 95% of a
    # whole-chromosome target VCF's calls fall outside the window, and counting
    # those as "invisible" would report a catastrophic-looking number that is
    # purely an artifact of the fixture. Both figures are emitted; the in-region
    # one is the meaningful one.
    lo, hi = int(s.pos.min()), int(s.pos.max())

    fill = -1 if absent_policy == "missing" else 0
    path = np.full((s.n_sites, 2), fill, dtype=np.int16)
    matched = np.zeros(s.n_sites, dtype=bool)

    n_target_records = 0
    n_nonref = n_nonref_in = 0
    exact_hit = pos_only_hit = invisible = 0        # over ALL target calls
    ex_in = po_in = inv_in = 0                      # restricted to [lo, hi]
    info_exact = info_pos_only = info_invisible = 0.0
    info_ex_in = info_po_in = info_inv_in = 0.0
    col = None

    with _open(target_vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                if sample not in hdr[9:]:
                    raise ValueError(f"{target_vcf} has no sample {sample!r}")
                col = 9 + hdr[9:].index(sample)
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            n_target_records += 1
            pos, ref, alt = int(f[1]), f[3], f[4]
            a, b = _parse_gt(f[col])

            row = strict.get((pos, ref, alt))
            if row is not None:
                path[row, 0], path[row, 1] = a, b
                matched[row] = True

            if (a > 0) or (b > 0):          # a non-reference call by the target
                n_nonref += 1
                # allele frequency as the panel reports it, for the -log f weighting
                fr = None
                for kv in f[7].split(";"):
                    if kv.startswith("AF="):
                        try:
                            fr = min(float(x) for x in kv[3:].split(",") if x)
                        except ValueError:
                            fr = None
                        break
                info = -math.log(fr) if (fr and 0 < fr < 1) else 0.0
                in_region = lo <= pos <= hi
                if in_region:
                    n_nonref_in += 1
                if row is not None:
                    exact_hit += 1;  info_exact += info
                    if in_region: ex_in += 1; info_ex_in += info
                elif pos in pos_seen:
                    pos_only_hit += 1; info_pos_only += info
                    if in_region: po_in += 1; info_po_in += info
                else:
                    invisible += 1;  info_invisible += info
                    if in_region: inv_in += 1; info_inv_in += info

    total_info = info_exact + info_pos_only + info_invisible
    total_in = info_ex_in + info_po_in + info_inv_in
    rep = {
        "sample": sample,
        "absent_policy": absent_policy,
        "n_sites": s.n_sites,
        "cohort_region": [lo, hi],
        "sites_matched_exactly": int(matched.sum()),
        "sites_unmatched": int((~matched).sum()),
        "sites_unmatched_frac": float((~matched).mean()),
        "target_records_read": n_target_records,
        "target_nonref_calls": n_nonref,
        "nonref_exact_record_in_G": exact_hit,
        "nonref_pos_only_different_spelling": pos_only_hit,
        "nonref_invisible_no_record": invisible,
        "nonref_exact_frac": exact_hit / n_nonref if n_nonref else None,
        "nonref_invisible_frac": invisible / n_nonref if n_nonref else None,
        "information_nats": {
            "exact": info_exact, "pos_only": info_pos_only, "invisible": info_invisible,
            "invisible_frac": (info_invisible / total_info) if total_info else None,
        },
        "in_region": {
            "note": ("THE MEANINGFUL FIGURES. Restricted to the coordinate range the "
                     "cohort actually covers, so a partial-chromosome fixture does not "
                     "report its own window as a representability failure."),
            "target_nonref_calls": n_nonref_in,
            "exact_record_in_G": ex_in,
            "pos_only_different_spelling": po_in,
            "invisible_no_record": inv_in,
            "exact_frac": ex_in / n_nonref_in if n_nonref_in else None,
            "invisible_frac": inv_in / n_nonref_in if n_nonref_in else None,
            "information_nats": {"exact": info_ex_in, "pos_only": info_po_in,
                                 "invisible": info_inv_in,
                                 "invisible_frac": (info_inv_in / total_in) if total_in else None},
        },
        "ceiling_note": (
            "invisible calls can never be retained at any tau -- they have no chain "
            "position. Report target_fidelity against this ceiling. A held-out panel "
            "sample UNDERSTATES it, since its variants are catalogued by construction."
        ),
    }

    inputs = {"target_vcf": store.sha256_file(target_vcf),
              "cohort_sites": store.sha256_file(cohort_dir / "sites.tsv")}
    params = {"sample": sample, "absent_policy": absent_policy, "join": "strict (POS,REF,ALT)"}
    store.save(outdir / "path.npy", path, inputs=inputs, params=params,
               axis_digest=s.digest, note="the PRIVATE INPUT: target allele per site, per strand")
    store.save(outdir / "representability.json", rep, inputs=inputs, params=params,
               axis_digest=s.digest, note="AUDIT ONLY -- never used to extend output_support")
    return rep
