# configs

Every file in this directory, verified against disk 2026-10-08.

| file | what |
|---|---|
| `inputs.chr21.yaml` | the input manifest: URLs, sizes, sha256, `chain_span`, target-selection rule, and what is deliberately NOT an input. ⚠ Its `chain_span` block is **superseded by decision C** — read the warning at the top of the block |
| `eligible_targets.chr21.txt` | the 3,085 samples eligible to stand in as a test target (Arm 2 only) |
| `clean_pool.chr21.txt` | the strict 1,399 — **the pool actually used** (Arm 2 only) |
| `targets.chr21.txt` | the 10 selected test targets, two per super-population (Arm 2 only) |
| `README.md` | this file |

⚠ **The table above used to list only two of these five**, while the table below called
them "the four sample lists" and named a fourth file that is not in this directory at all.
Both are fixed: the table above now lists every file in `configs/`, and the fourth list's
real location is stated where it appears.

⚠ **Nothing in this directory is per-chromosome except the manifest.** All three list
files carry a `.chr21` suffix and **contain nothing chr21-specific** — they are sample ids
selected on genome-wide relatedness and genome-wide pedigree, so the suffix is misleading.
Under decision D (22 autosomes) the chromosome loop **must not copy or regenerate them per
chromosome**; one copy serves every chromosome. They are left named as they are because
the manifest's `pool:`, `wide_pool:` and `selected:` keys point at these exact filenames
and renaming them is a change to code-visible paths, not to documentation.

⚠ **Arm 1 needs a fifth list that does not exist: the 44 HPRC cohort members.** Arm 1's
targets ARE the cohort, so none of the three lists above applies to it — every one of them
is an Arm 2 construct (a held-out external sample drawn from the 1000G 30x panel). The
Arm 1 target list is **PLANNED**; until it exists, no file here answers "which 44".

## Why `eligible_targets` is not just "the panel minus the cohort"

The 39 HPRC donors present in the 1000G 30x panel are the **children** of trios,
and their **78 parents sit in the unrelated-2504 set**. So a target drawn at
random from the 2504 has roughly a 3.1% chance of being a cohort member's parent,
sharing about half its genome with a donor haplotype — which would inflate
apparent fidelity and confound the privacy measurement.

3,202 − 39 (HPRC in panel) − 78 (first-degree relatives) = **3,085**.

Regenerate with the 1000G pedigree at
`.../1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt`.

## The four sample lists — three here, one in the data tree

All four are **Arm 2** artifacts: they exist to pick and protect a stand-in *external*
target. Arm 1 draws no external sample.

| file | n | what |
|---|---|---|
| `eligible_targets.chr21.txt` | 3,085 | panel minus HPRC members and their first-degree relatives |
| `clean_pool.chr21.txt` | 1,399 | **use this one** — also excludes anyone with a first-degree relative *anywhere* in the 3,202 |
| `targets.chr21.txt` | 10 | the selected test targets, two per super-population |
| `data/attack/excluded_targets.txt` | 10 | the same ten, as consumed by the attack-DB build. ⚠ **It plays a configs role but lives in the data tree**, which is why it was in no `docs/data.md` row for weeks; it is now recorded in the storage-layout block of `docs/data.md` §1. Moving it into `configs/` would be the cleaner fix and is not done here — the attack-DB build reads it from `data/attack/`, so the move is a code change |

**Why the strict pool.** 1,686 of the 3,085 — 55% — have a first-degree relative
somewhere in the panel. If a target's own parent or child sits in the attack
database, the attack can match the relative instead of the target, which makes
privacy look better than it is.

## Reproducing the selection

The pedigree is at
`.../1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt` and the
population labels at
`.../release/20130502/integrated_call_samples_v3.20130502.ALL.panel`. Selection
used `random.Random(20260929)`, so it replays exactly.
