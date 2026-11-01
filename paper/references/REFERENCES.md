# Reference tracker

Every claim in the paper that needs a citation, and where the citation goes.
One row per CITATION SLOT, not per paper: a paper cited in three places has three rows.

- **Have** — the PDF is in this folder and the slot is covered.
- **Needed** — the slot is open; the sentence does not go into prose without it.
- **Standard** — textbook or universally known; no PDF kept, cite in the usual form.

File naming in this folder: `YYYY_FirstAuthor_venue-or-id_slug.pdf`.
Ours is `2026_OURS_*`; the head-to-head method is `2026_PanMixer_NatComms.pdf`.

## ⚠⚠ THERE IS NO BIBLIOGRAPHY YET

`find . -name '*.bib'` outside `external/` returns **nothing**. Every `\cite{...}` key in
`paper/intro_draft_2026-10-31.tex` is a **PLACEHOLDER**: it will not resolve, and the
candidate papers named in the table below are **my best identification from memory and are
UNVERIFIED**. Do not treat a key as a citation until its PDF is in this folder and its row
reads `Have`.

## Citation keys used in the Introduction draft

Source: `paper/intro_draft_2026-10-31.tex`.

| key | where | the claim it supports | status | candidate (UNVERIFIED unless `Have`) |
|---|---|---|---|---|
| `schneider2017` | 1.1 s1 | a linear reference is one arbitrary mosaic from few donors | **Needed** | Schneider et al. 2017, Genome Res — GRCh38 / de novo haploid assembly evaluation |
| `ballouz2019` | 1.1 s2 | reference bias is structural, and falls on under-represented ancestries | **Needed** | Ballouz, Dobin & Gillis 2019, Genome Biol — "Is it time to change the reference genome?" |
| `eizenga2020` | 1.1 s3, 1.2 s1 | what a pangenome graph is and what it supports | **Needed** | Eizenga et al. 2020, Annu Rev Genomics Hum Genet — "Pangenome Graphs". **DYLAN'S OWN KEY**, from his draft |
| `siren2021` | 1.1 s3, s4 | read alignment on a graph; alignment bias reduced | **Needed** | Sirén et al. 2021, Science — graph alignment / giraffe |
| `ebler2022` | 1.1 s3 | genotyping on a graph | **Needed** | Ebler et al. 2022, Nat Genet — PanGenie |
| `hickey2024` | 1.1 s4, 1.2 s4 | structural variation modelled natively; one of the two builders | **Needed** | Hickey et al., Nat Biotechnol — Minigraph-Cactus |
| `liao2023` | 1.2 s3, s4 | LD and phase along a haplotype path; the first human pangenome reference | ✅ **Have** | `2023_FirstPangenome.pdf` — Liao et al., "A draft human pangenome reference", Nature 617:312 (2023), doi:10.1038/s41586-023-05896-x |
| `garrison2024` | 1.2 s4 | the second independent construction method | **Needed** | Garrison et al., PGGB — "Building pangenome graphs" |
| `wang2022` | 1.2 s4 | the consortium continues to add individuals | **Needed** | Wang et al. 2022, Nature — "The Human Pangenome Project" |
| `erlich2014` | 1.3 s3 | a path is identifying; rare alleles are near-unique | **Needed** | Erlich & Narayanan 2014, Nat Rev Genet — "Routes for breaching and protecting genetic privacy" |
| `blindenbach2026` | 1.3 s4, s6 | rarity-weighted linkage; PanMixer as the existing method | ✅ **Have** | `2026_PanMixer_NatComms.pdf` |
| `venkatesaramani2023` | 1.3 s4 | membership inference on summary releases | ✅ **Have** | `2023_Venkatesaramani_arXiv2302.01763_beacons_preprint.pdf` |
| `browning2018` | 1.3 s4 | imputation-style reconstruction from panel LD | **Needed** | Browning, Zhou & Browning 2018, AJHG — Beagle 5 |
| `li2003` | 2.4 | the Li-Stephens haplotype-copying model | **Needed** | Li & Stephens 2003, Genetics — modelling LD and identifying recombination hotspots |
| `privvg2022` | 1.4 | a graph-level DP pangenome mechanism already exists | ⚠⚠ **Needed — BLOCKING** | privvg / `odgi priv`. See the rule below |

## Introduction slots

| # | Claim needing a reference | Where | Status | Note |
|---|---|---|---|---|
| 1 | A linear reference collapses population diversity into one arbitrary mosaic | 1.1 | **Needed** | `schneider2017`. No source supports the word "arbitrary"; it is Dylan's framing |
| 2 | Reference bias is structural, and the alleles lost are those common in under-represented ancestries | 1.1 | **Needed** | `ballouz2019` |
| 3 | Pangenome gains — alignment bias, structural variation, LD and phasing | 1.1, 1.2 | **Needed** | ONE CITATION PER CLAUSE, not one for the group. Split across `siren2021`, `hickey2024`, `liao2023` |
| 4 | Downstream payoff — read alignment, genotyping, variant calling | 1.1 | **Needed** | `siren2021`, `ebler2022`, plus one for variant calling |
| 5 | The first human pangenome reference, released with two independent builders | 1.2 | ✅ **Have** | `liao2023`; the builders themselves still need `garrison2024` and `hickey2024` |
| 6 | The consortium continues to add individuals | 1.2 | **Needed** | `wang2022` |
| 7 | A haplotype path is identifying | 1.3 | **Needed** | `erlich2014` |
| 8 | Participation and consent — people may decline without protection | ⚠ NOT IN THE DRAFT | **Needed** | The sentence was drafted and then CUT, because no source we hold supports it. Dylan's decision whether the point belongs in 1.3 at all |
| 9 | Rarity-weighted linkage against an external panel | 1.3 | ✅ **Have** | `blindenbach2026` |
| 10 | Membership inference; the Beacon narrowing does not remove the exposure | 1.3 | ✅ **Have** | `venkatesaramani2023` |
| 11 | Imputation-style reconstruction from panel LD | 1.3 | **Needed** | `browning2018` |
| 12 | An existing method selects edits under a utility budget | 1.3 | ✅ **Have** | `blindenbach2026` |
| 13 | A graph-level differentially private pangenome mechanism already exists | 1.4 | ⚠⚠ **Needed — BLOCKING** | `privvg2022` |
| 26 | Non-reference sequence placed in its correct context against flanking homology | ⚠ NOT IN THE DRAFT | **Needed** | Was in the bullet scaffold, was OURS rather than Dylan's, and would be a fifth clause in an already long sentence. Available if he wants it |

## Methods, Results, Discussion slots

| # | Claim needing a reference | Where | Status | Note |
|---|---|---|---|---|
| 14 | The graph build both arms run on | 2.1 | ✅ **Have** | `liao2023` — `2023_FirstPangenome.pdf`. The build, stated once, with its counts |
| 15 | Li–Stephens haplotype copying model | 2.4 | **Needed** | `li2003`. ⚠ The TRANSITION CONSTANTS (N_e = 10,000, r = 1.26) come from PanMixer’s published Methods, not from Li & Stephens, so 2.4 needs BOTH citations |
| 16 | Forward filtering, backward sampling | 2.4 | Standard | — |
| 17 | The privacy guarantee this paper implements | 2.6 | ✅ **Have** | `2026_OURS_Private_Genome_Path_Release.pdf` — Theorem 1, Corollary 1, Definition 1 |
| 18 | Pinsker's inequality | 2.9 | Standard | — |
| 19 | The exponential mechanism family | 2.3, 4.5 | **Needed** | McSherry & Talwar |
| 20 | A bounded-rational Bayesian attacker is worse for the defender than the likelihood-ratio-optimal test | 2.12, 4.3 | ✅ **Have** | `2024_Zhang_Malin_Vorobeychik_arXiv2406.01811_game_theoretic.pdf` |
| 21 | Suppression and noise against likelihood-ratio membership inference on AF releases | 4.5 | ✅ **Have** | `2023_Venkatesaramani_arXiv2302.01763_beacons_preprint.pdf` |
| 22 | PanMixer's own evaluators, run unmodified on our releases | 2.12, 3.3 | ✅ **Have** | `2026_PanMixer_NatComms.pdf` |
| 23 | 0–1 knapsack formulation of their optimiser | 4.4 | ✅ **Have** | `2026_PanMixer_NatComms.pdf` |
| 24 | Their stated scope — "the guarantees are empirical" — and their naming of formal guarantees as future work | 4.4, 4.10 | ✅ **Have** | `2026_PanMixer_NatComms.pdf`, quoted once |
| 25 | Their DP paragraph cites a pangenome-informed language model for synthetic sequence generation, NOT a graph-noising method | 4.5 | **Needed** | the paper they cite, so we answer their sentence and not one they never made |

## Rules

- A **Needed — blocking** row stops a claim, not a draft: write the sentence, leave the slot,
  do not soften the claim to avoid the citation.
- **Slot 13 is the one that constrains a claim of ours.** Until privvg's primary source is read,
  nothing may say or imply "first with a guarantee", and "to our knowledge" may not be attached
  to the novelty sentence in 1.4. The draft obeys this.
- A reference added here gets its PDF in this folder in the same change, or the row stays
  **Needed**.
- A candidate named in this file is an identification from memory and is not evidence. Verify
  against the PDF before it is typeset.
