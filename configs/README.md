# configs

| file | what |
|---|---|
| `inputs.chr21.yaml` | the input manifest: URLs, sizes, sha256, `chain_span`, target-selection rule, and what is deliberately NOT an input |
| `eligible_targets.chr21.txt` | the 3,085 samples eligible to stand in as a test target |

## Why `eligible_targets` is not just "the panel minus the cohort"

The 39 HPRC donors present in the 1000G 30x panel are the **children** of trios,
and their **78 parents sit in the unrelated-2504 set**. So a target drawn at
random from the 2504 has roughly a 3.1% chance of being a cohort member's parent,
sharing about half its genome with a donor haplotype — which would inflate
apparent fidelity and confound the privacy measurement.

3,202 − 39 (HPRC in panel) − 78 (first-degree relatives) = **3,085**.

Regenerate with the 1000G pedigree at
`.../1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt`.

## The four sample lists

| file | n | what |
|---|---|---|
| `eligible_targets.chr21.txt` | 3,085 | panel minus HPRC members and their first-degree relatives |
| `clean_pool.chr21.txt` | 1,399 | **use this one** — also excludes anyone with a first-degree relative *anywhere* in the 3,202 |
| `targets.chr21.txt` | 10 | the selected test targets, two per super-population |
| `data/attack/excluded_targets.txt` | 10 | the same ten, as consumed by the attack-DB build |

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
