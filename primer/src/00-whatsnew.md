### What changed in this revision

You stopped at §A.6. This revision folds in everything the working sessions established since then, and from §A.7 onward it adds the questions you are likely to ask next.

**Behind your bookmark, and worth going back for**

- **§1.1.8 (new)** explains how a graph and a VCF really correspond: indels, bubbles nested inside bubbles and the `LV` / `PS` / `AT` tags that record them, all five reasons a haplotype can have no value at a site, cycles, and the chr21 site where every haplotype differs. Steps 9 and 12 in Part B lean on it.
- **§1.1.4** gains a Caution: LD blocks are not disjoint in the way the original text implied.
- **§2.3.4** writes out the ranking attack formally (p, q, Q_p, Q_q, y), with a worked three-candidate example. **§2.4.3** explains why reconstruction works even though the space of possible paths is astronomically large.
- **§A.3–§A.6** are revised. They now give the privacy score for a block with several variants (§A.4.1), show that the score reads anchor alleles only (§A.4.2), and retract the claim that `−log 0` zeroes the rarest alleles (§A.4.3). They also cover what a move that changes nothing costs (§A.5) and what the optimizer credits it with anyway (§A.6).

**From your bookmark on**

- **§A.8–§A.10 (new)** cover what an anchor actually is and whether a design needs one, which dataset supplies f_v on the frequency route, and what a missing genotype does to the privacy score and to the cost.
- **Part B**: Steps 9, 12, 13, 14, 17, 18 and 26 are substantially revised, and there is a new **worked trace of block 1077**, from raw records to the released VCF.
- **Departures (f)–(j) (new)** are five more findings. Some are places where the released code differs from the published Methods; others are metrics that measure something narrower than their names suggest, or consequences of the data model that the Methods never take up. A closing section then applies the same test to our own pipeline.
- **Section 4, Part C (new)** holds the pinned v1 specification (§4.16) and the four candidate missing-data policies (§4.17). The chain-length question and the choice of missing-data policy are both left open. The primer lays out what each decision turns on and does not make it.

**The violet "Question" boxes** mostly begin at §A.7. From there on, each one predicts a question you are likely to ask, based on how you have asked so far, and answers it in full. The few in §1.1.8, §2.3.4 and §2.4.3 are different: they are questions you actually asked in the working sessions, answered where the primer should have answered them.
