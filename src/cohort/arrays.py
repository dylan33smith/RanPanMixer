"""Step 3: the cohort VCF -> the arrays. THE IRREVERSIBLE STEP.

This is the only preprocessing step with real design content, and the one that
will be rewritten most, which is why it is alone in its own module.

Irreversible because the VCF line is the only place four things coexist: the
genotype, the INFO tags (LV, PS, CONFLICT, AN), the REF/ALT strings, and the
record's position in the file. A converter that keeps only position and genotype
-- which is what the upstream one does -- destroys the rest permanently.
Cheap now, unrecoverable later.

TWO PASSES, because "inherited" is not a property of the record it applies to:
  pass 1  matrix, per-record facts, lengths, sites, and an ID -> row index
  pass 2  resolve the reason codes, which need the PARENT record's call

PRECEDENCE, pinned in docs/terms.md under `reason_array`:
  conflict (3) > not applicable (1) > inherited (2) > uncategorised (4)
Conflict wins because the CONFLICT tag is the caller stating outright that it
could not choose, which beats a structural inference.

ASSEMBLY GAPS GET NO CODE. They are a property of a RUN, not of a cell, and
classifying them per cell needs a length threshold that has never been validated
-- it would claim 92.06% of all missing cells at >=10 and 68.67% at >=1000. They
go to `missing_runs.tsv` as intervals so any threshold can be applied, changed or
swept downstream without re-reading the VCF.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np

import sites as sites_mod
import store

CALLED         = 0
NOT_APPLICABLE = 1   # nested child; this haplotype's parent allele lacks the bubble
INHERITED      = 2   # the parent record's own call was missing
CONFLICT       = 3   # the CONFLICT tag names this sample
UNCATEGORISED  = 4   # a no-call with no structural explanation

CODE_NAMES = {CALLED: "called", NOT_APPLICABLE: "not_applicable",
              INHERITED: "inherited", CONFLICT: "conflict",
              UNCATEGORISED: "uncategorised"}

_INFO_RE = re.compile(r"(^|;)(LV|PS|CONFLICT|AN)=([^;]*)")


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _parse_info(info: str) -> dict:
    out = {}
    for _, key, val in _INFO_RE.findall(info):
        out[key] = val
    return out


def _parse_gt(field: str) -> tuple[int, int]:
    """One GT -> (allele0, allele1), with -1 for missing.

    A separator-less GT is written into strand 0 with -1 in strand 1. On chr21
    every such field is an explicit '.', so this branch only ever produces
    (-1, -1) there -- but the asymmetry is preserved because it is real for any
    other input VCF, and because the two haplotype columns are then not
    exchangeable, which everything downstream implicitly assumes they are.
    """
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


def build(cohort_vcf: Path, outdir: Path, *, chrom: str | None = None) -> dict:
    """Emit haplotypes, reason, missing_runs, allele_lengths, sites, ids.

    Takes NO target argument, and must not acquire one -- see cohort/__init__.py.
    """
    cohort_vcf, outdir = Path(cohort_vcf), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- pass 0: sample names and record count -----------------------------
    samples, n_records = None, 0
    with _open(cohort_vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                n_records += 1
    if samples is None:
        raise ValueError(f"{cohort_vcf}: no #CHROM header line")
    n_samples = len(samples)
    n_hap = 2 * n_samples

    hap = np.full((n_hap, n_records), -1, dtype=np.int16)
    reason = np.full((n_hap, n_records), UNCATEGORISED, dtype=np.int8)
    positions = np.zeros(n_records, dtype=np.int64)
    an_info = np.full(n_records, -1, dtype=np.int32)

    site_rows = []
    id_to_row: dict[str, int] = {}
    ps_of_row: list[str] = []
    conflict_of_row: list[set] = []
    lengths_flat: list[int] = []
    lengths_off = np.zeros(n_records + 1, dtype=np.int64)

    sample_ix = {s: i for i, s in enumerate(samples)}

    # ---- pass 1: stream the VCF -------------------------------------------
    r = 0
    with _open(cohort_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            c, pos, vid, ref, alt, info = f[0], int(f[1]), f[2], f[3], f[4], f[7]
            alts = alt.split(",")
            inf = _parse_info(info)

            positions[r] = pos
            if vid and vid != ".":
                id_to_row.setdefault(vid, r)
            ps_of_row.append(inf.get("PS", ""))
            conflict_of_row.append(
                set(inf["CONFLICT"].split(",")) if inf.get("CONFLICT") else set()
            )
            an_info[r] = int(inf["AN"]) if inf.get("AN", "").isdigit() else -1

            lengths_off[r] = len(lengths_flat)
            lengths_flat.append(len(ref))
            lengths_flat.extend(len(a) for a in alts)

            site_rows.append((c, pos, ref, alt, len(alts),
                              int(inf.get("LV", -1)) if inf.get("LV", "").lstrip("-").isdigit() else -1,
                              inf.get("PS", "")))

            for i in range(n_samples):
                a, b = _parse_gt(f[9 + i])
                hap[2 * i, r] = a
                hap[2 * i + 1, r] = b
            r += 1
    lengths_off[n_records] = len(lengths_flat)
    assert r == n_records, f"pass-1 read {r} records, pass-0 counted {n_records}"

    # ---- pass 2: resolve reason codes -------------------------------------
    # Order matters and is the pinned precedence. Start everything CALLED,
    # then overwrite the missing cells from lowest precedence upward so the
    # highest-precedence rule is the last writer.
    reason[:] = CALLED
    missing = hap == -1
    reason[missing] = UNCATEGORISED

    parent_row = np.full(n_records, -1, dtype=np.int64)
    for i, ps in enumerate(ps_of_row):
        if ps and ps in id_to_row:
            parent_row[i] = id_to_row[ps]

    nested = np.flatnonzero(parent_row >= 0)
    for i in nested:
        p = parent_row[i]
        miss_here = missing[:, i]
        if not miss_here.any():
            continue
        parent_missing = hap[:, p] == -1
        # inherited first, then not-applicable overwrites where the parent WAS
        # called -- i.e. the DNA genuinely is not on that chromosome.
        reason[miss_here & parent_missing, i] = INHERITED
        reason[miss_here & ~parent_missing, i] = NOT_APPLICABLE

    # conflict last: highest precedence
    for i, names in enumerate(conflict_of_row):
        if not names:
            continue
        for s in names:
            j = sample_ix.get(s)
            if j is None:
                continue
            for h in (2 * j, 2 * j + 1):
                if hap[h, i] == -1:
                    reason[h, i] = CONFLICT

    assert not (reason[~missing] != CALLED).any(), "a called cell got a missing-reason code"
    assert not (reason[missing] == CALLED).any(), "a missing cell was left as CALLED"

    # ---- missing runs (intervals, no threshold) ---------------------------
    runs = []
    for h in range(n_hap):
        m = missing[h]
        if not m.any():
            continue
        d = np.diff(np.concatenate(([0], m.view(np.int8), [0])))
        starts = np.flatnonzero(d == 1)
        ends = np.flatnonzero(d == -1)
        for s, e in zip(starts, ends):
            runs.append((samples[h // 2], h % 2, int(s), int(e)))

    # ---- write -------------------------------------------------------------
    src_digest = store.sha256_file(cohort_vcf)
    inputs = {"cohort_vcf": src_digest}
    params = {"n_samples": n_samples, "n_haplotypes": n_hap, "n_sites": n_records,
              "precedence": "conflict>not_applicable>inherited>uncategorised"}

    digest = sites_mod.write_sites_tsv(outdir / "sites.tsv", site_rows)
    store.attach_prov(outdir / "sites.tsv", inputs=inputs, params=params,
                      axis_digest=digest,
                      note="THE site axis: what row i means. Every other artifact is "
                           "indexed by it and carries this digest.")

    store.save(outdir / "haplotypes.npy", hap, inputs=inputs, params=params,
               axis_digest=digest, note="(n_haplotypes, n_sites) int16, -1 = missing")
    store.save(outdir / "reason.npy", reason, inputs=inputs,
               params={**params, "codes": {str(k): v for k, v in CODE_NAMES.items()}},
               axis_digest=digest, note="why each -1 is missing; assembly gaps are NOT coded here")
    store.save(outdir / "positions.npy", positions, inputs=inputs, params=params,
               axis_digest=digest)
    store.save(outdir / "an_info.npy", an_info, inputs=inputs, params=params,
               axis_digest=digest, note="INFO/AN as written by the caller; cross-checked against support")
    store.save(outdir / "allele_lengths.npy", np.asarray(lengths_flat, dtype=np.int32),
               inputs=inputs, params=params, axis_digest=digest,
               note="ragged; record i occupies [offsets[i], offsets[i+1])")
    store.save(outdir / "allele_lengths_offsets.npy", lengths_off, inputs=inputs,
               params=params, axis_digest=digest)
    store.save(outdir / "haplotype_ids.npy",
               np.array([f"{samples[h // 2]}:{h % 2}" for h in range(n_hap)], dtype="<U32"),
               inputs=inputs, params=params, axis_digest=digest,
               note="haplotype h is (sample h//2, strand h%2)")

    runs_tsv = "sample\tstrand\tstart_row\tend_row\n" + "".join(
        f"{s}\t{st}\t{a}\t{b}\n" for s, st, a, b in runs)
    store.save(outdir / "missing_runs.tsv", runs_tsv, inputs=inputs,
               params={**params, "threshold": None},
               axis_digest=digest,
               note="every run of consecutive -1 as a half-open interval. No length "
                    "threshold is applied here, deliberately: the >=100 cut has never "
                    "been validated and claims 92.06%/68.67% of missing cells at >=10/>=1000.")

    summary = {
        "n_samples": n_samples, "n_haplotypes": n_hap, "n_sites": n_records,
        "axis_digest": digest,
        "missing_cells": int(missing.sum()),
        "missing_frac": float(missing.mean()),
        "reason_counts": {CODE_NAMES[k]: int((reason == k).sum()) for k in CODE_NAMES},
        "n_missing_runs": len(runs),
        "n_nested_records": int((parent_row >= 0).sum()),
        "n_conflict_records": int(sum(1 for c in conflict_of_row if c)),
        "n_declared_alleles": len(lengths_flat),
    }
    return summary
