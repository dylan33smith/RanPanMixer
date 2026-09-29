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
