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

# Why the TARGET needs its own reason codes, separate from the cohort's.
#
# A -1 in path.npy has two completely different meanings and phi_v cannot tell
# them apart, because both simply produce "no tilt here":
#
#   NO_RECORD  the target's callset has no record at this graph position at all.
#              A limitation of OUR JOIN, identical for every target drawn from the
#              same panel. It says nothing about the person.
#   NO_CALL    the callset HAS a record here but reported no genotype for this
#              individual. A real failure to call THIS person.
#
# They behave identically in the mechanism and differently in the reporting: a
# position the panel never covered arguably should not count against
# target_fidelity, while one it covered and could not call for this person
# arguably should. Without these codes that denominator is ambiguous, and the
# cohort's reason codes cannot express either case -- theirs are about graph
# nesting (LV/PS/CONFLICT), which has nothing to do with callset coverage.
T_CALLED    = 0
T_NO_RECORD = 1
T_NO_CALL   = 2
T_CODE_NAMES = {T_CALLED: "called", T_NO_RECORD: "no_record_in_target_callset",
                T_NO_CALL: "record_present_but_no_call"}


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
          sample: str, absent_policy: str = "missing",
          split_alt_join: bool = True) -> dict:
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
    # ⚠ check_axis was dead code until 2026-09-29 -- defined, unit-tested, and never
    # called from src/, so a cohort directory mixing two site axes was consumed
    # without error. Verify every cohort array this stage will read.
    for art in ("in_chain.npy", "positions.npy", "support.npy"):
        f = cohort_dir / art
        if f.exists():
            store.check_axis(f, s.digest)
    strict: dict[tuple[int, str, str], int] = {}
    pos_seen: dict[int, int] = {}
    # SPLIT-ALT INDEX. The 1000G panel is entirely decomposed into biallelic
    # records -- every one of its 1,002,752 chr21 records declares exactly one ALT
    # -- while 7.9% of ours are multi-allelic. A strict full-ALT-string match can
    # therefore NEVER succeed on those 27,021 records, for a purely notational
    # reason. That threw away target information we actually hold: 8,377 records
    # are recoverable by matching each of OUR alts against its own panel record and
    # mapping back to OUR allele numbering. Keyed (pos, ref, one_alt) -> (row, j)
    # where j is the 1-based index of that alt in our record.
    split_ix: dict[tuple[int, str, str], tuple[int, int]] = {}
    for i in range(s.n_sites):
        pos_i, ref_i = int(s.pos[i]), s.ref[i]
        strict.setdefault((pos_i, ref_i, s.alt[i]), i)
        pos_seen.setdefault(pos_i, i)

    # The cohort's coordinate range. Representability must be reported RELATIVE
    # to the region actually modelled: on a 2 Mb development slice, 95% of a
    # whole-chromosome target VCF's calls fall outside the window, and counting
    # those as "invisible" would report a catastrophic-looking number that is
    # purely an artifact of the fixture. Both figures are emitted; the in-region
    # one is the meaningful one.
    lo, hi = int(s.pos.min()), int(s.pos.max())

    # ⚠ FULL-COVERAGE REQUIREMENT, and why the obvious version of this is unsafe.
    #
    # The reconstruction infers "this haplotype is REFERENCE" from the absence of a
    # 1 at every panel record for our alts. That inference is only valid if the
    # panel actually HAS a record for every one of our alts. If it covers 2 of our
    # 4, a haplotype carrying alt 3 shows 0 at both covered records and would be
    # written as reference -- a FABRICATED homozygous-reference call. The first
    # version of this code did exactly that for 6,309 rows per target.
    #
    # So the recoverable set is restricted to rows the panel FULLY covers, and --
    # critically -- it is computed from the PANEL'S SITE LIST ALONE, before any
    # target is read. That keeps the set identical for every target, which is the
    # precondition for normalising beta_v over it; deciding it per-person instead
    # would make beta_v target-dependent and void Theorem 1.
    fully_covered: set[int] = set()
    if split_alt_join:
        want: dict[tuple[int, str], set[str]] = {}
        for i in range(s.n_sites):
            if s.n_alt[i] > 1:
                want.setdefault((int(s.pos[i]), s.ref[i]), set()).update(s.alt[i].split(","))
        have: dict[tuple[int, str], set[str]] = {}
        with _open(target_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 6)
                k = (int(f[1]), f[3])
                if k in want:
                    have.setdefault(k, set()).add(f[4])
        for i in range(s.n_sites):
            if s.n_alt[i] > 1:
                k = (int(s.pos[i]), s.ref[i])
                alts = set(s.alt[i].split(","))
                if alts <= have.get(k, set()):
                    fully_covered.add(i)
                    for j, a in enumerate(s.alt[i].split(","), start=1):
                        split_ix.setdefault((k[0], k[1], a), (i, j))

    # the panel's site keys, for the target-independent `readable` set below
    panel_site_keys: set = set()
    with _open(target_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 6)
            panel_site_keys.add((int(f[1]), f[3], f[4]))
    inputs_panel = {"target_vcf_sitelist": "see path.npy inputs"}

    fill = -1 if absent_policy == "missing" else 0
    path = np.full((s.n_sites, 2), fill, dtype=np.int16)
    matched = np.zeros(s.n_sites, dtype=bool)
    # every site starts as "no record"; the join downgrades to no-call or promotes
    # to called as it learns better
    treason = np.full(s.n_sites, T_NO_RECORD, dtype=np.int8)
    # row -> {our_allele_j: (gt0, gt1)} accumulated from the panel's split records
    split_hits: dict[int, dict[int, tuple[int, int]]] = {}
    # non-ref target calls whose classification must wait for the split resolution:
    # row -> [(info, in_region), ...]
    split_pending: dict[int, list] = {}

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
                treason[row] = T_NO_CALL if (a == -1 and b == -1) else T_CALLED
            elif split_ix:
                hit = split_ix.get((pos, ref, alt))
                if hit is not None:
                    r2, j = hit
                    split_hits.setdefault(r2, {})[j] = (a, b)

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
                elif (pos, ref, alt) in split_ix:
                    # defer: representable IF the reconstruction succeeds
                    split_pending.setdefault(split_ix[(pos, ref, alt)][0], []).append(
                        (info, in_region))
                elif pos in pos_seen:
                    pos_only_hit += 1; info_pos_only += info
                    if in_region: po_in += 1; info_po_in += info
                else:
                    invisible += 1;  info_invisible += info
                    if in_region: inv_in += 1; info_inv_in += info

    # resolve the split-ALT reconstruction. For each haplotype, the panel record
    # whose single ALT the haplotype carries names OUR allele index; if no panel
    # record for this site says 1, the haplotype is reference.
    n_recovered = 0
    n_contradictions = 0
    recovered_calls = [0, 0]      # [all, in-region] non-ref calls credited via split
    for r2, byj in split_hits.items():
        # only reachable for fully covered rows, so "no alt carries a 1" really does
        # mean reference rather than "we could not see the alt this haplotype has"
        out = [fill, fill]
        ok = True
        for h in (0, 1):
            carried = [j for j, gt in byj.items() if gt[h] == 1]
            if len(carried) > 1:
                ok = False
                break
            if carried:
                out[h] = carried[0]
            elif any(gt[h] == -1 for gt in byj.values()):
                out[h] = -1
            else:
                out[h] = 0          # conclusive: every alt is covered and none is carried
        if not ok:
            n_contradictions += 1
            # not representable after all -- a record exists at that POS but we
            # cannot resolve which of our alleles the haplotype carries
            for info, in_region in split_pending.pop(r2, []):
                pos_only_hit += 1; info_pos_only += info
                if in_region: po_in += 1; info_po_in += info
            continue
        path[r2, 0], path[r2, 1] = out
        matched[r2] = True
        treason[r2] = T_NO_CALL if (out[0] == -1 and out[1] == -1) else T_CALLED
        n_recovered += 1
        # the release CAN express these calls, so they belong with the exact matches
        for info, in_region in split_pending.pop(r2, []):
            exact_hit += 1; info_exact += info
            recovered_calls[0] += 1
            if in_region:
                ex_in += 1; info_ex_in += info
                recovered_calls[1] += 1

    for leftover in split_pending.values():
        for info, in_region in leftover:
            pos_only_hit += 1; info_pos_only += info
            if in_region: po_in += 1; info_po_in += info

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
        "target_reason_counts": {T_CODE_NAMES[k]: int((treason == k).sum()) for k in T_CODE_NAMES},
        "split_alt_join": {
            "enabled": split_alt_join,
            "our_multiallelic_records": int((s.n_alt > 1).sum()),
            "fully_covered_by_panel": len(fully_covered),
            "rows_recovered": n_recovered,
            "contradictions": n_contradictions,
            "nonref_calls_credited": recovered_calls[0],
            "nonref_calls_credited_in_region": recovered_calls[1],
            "contradiction_cause": ("observed at poly-A repeat sites where the panel's "
                                    "biallelic records are NOT a mutually exclusive "
                                    "decomposition -- one haplotype reads as carrying two "
                                    "different alts. The row is left unresolved rather than "
                                    "guessed, and its calls fall back to pos_only."),
            "note": ("The panel is decomposed into biallelic records while 7.9% of our "
                     "records are multi-allelic, so a strict full-ALT match cannot succeed "
                     "on those. These rows were recovered by matching each of our alts "
                     "against its own panel record. 'contradictions' counts rows where a "
                     "haplotype appeared to carry two different alts, which should not "
                     "happen at a properly decomposed site."),
        },
        "in_chain_usable": None,     # filled below
        "ceiling_note": (
            "invisible calls can never be retained at any tau -- they have no chain "
            "position. Report target_fidelity against this ceiling. A held-out panel "
            "sample UNDERSTATES it, since its variants are catalogued by construction."
        ),
    }

    # THE TARGET-INDEPENDENT SET. A row is "readable" if the target's panel has a
    # record our join can use: either a strict (POS,REF,ALT) hit or a fully covered
    # multi-allelic row. Both tests read only the panel's SITE LIST, which is the
    # same file for every target, so this array is identical for all of them --
    # verified by the audit check. That invariance is what makes it legal to
    # normalise beta_v over this set: doing it over "positions where THIS target has
    # a call" would be target-dependent and void Theorem 1, and it is NOT the same
    # set, because a contradiction at a poly-A repeat drops a row for one target and
    # not another.
    readable = np.zeros(s.n_sites, dtype=bool)
    for i in range(s.n_sites):
        if (int(s.pos[i]), s.ref[i], s.alt[i]) in panel_site_keys or i in fully_covered:
            readable[i] = True
    store.save(outdir / "readable.npy", readable, inputs=inputs_panel,
               params={"definition": "panel has a usable record for this row",
                       "target_independent": True},
               axis_digest=s.digest,
               note="TARGET-INDEPENDENT. Normalise beta_v over `in_chain & readable`, "
                    "never over positions where a particular target happens to have a "
                    "call -- that set varies between targets and would make beta_v "
                    "target-dependent.")

    # how much of the tilt budget this target can actually move
    usable = s.in_chain & (treason == T_CALLED)
    rep["in_chain_usable"] = {
        "chain_positions": int(s.in_chain.sum()),
        "with_a_target_call": int(usable.sum()),
        "readable_target_independent": int((s.in_chain & readable).sum()),
        "frac": float(usable.sum() / s.in_chain.sum()) if s.in_chain.any() else None,
        "readable_frac": float((s.in_chain & readable).sum() / s.in_chain.sum()) if s.in_chain.any() else None,
        "note": ("positions where phi_v is defined. The rest consume beta_v weight and "
                 "can never contribute to u, so max achievable u is this fraction unless "
                 "beta is normalised over this set instead -- see docs/plan.md."),
    }

    inputs = {"target_vcf": store.sha256_file(target_vcf),
              "cohort_sites": store.sha256_file(cohort_dir / "sites.tsv")}
    params = {"sample": sample, "absent_policy": absent_policy,
              "join": "strict (POS,REF,ALT)" + (" + split-ALT recovery" if split_alt_join else ""),
              "split_alt_join": split_alt_join}
    store.save(outdir / "path.npy", path, inputs=inputs, params=params,
               axis_digest=s.digest, note="the PRIVATE INPUT: target allele per site, per strand")
    store.save(outdir / "target_reason.npy", treason, inputs=inputs,
               params={**params, "codes": {str(k): v for k, v in T_CODE_NAMES.items()}},
               axis_digest=s.digest,
               note="why a path entry is -1: our join had no record (1) vs the callset "
                    "could not call this person (2). phi_v cannot distinguish them; "
                    "target_fidelity's denominator depends on which.")
    store.save(outdir / "representability.json", rep, inputs=inputs, params=params,
               axis_digest=s.digest, note="AUDIT ONLY -- never used to extend output_support")
    return rep
