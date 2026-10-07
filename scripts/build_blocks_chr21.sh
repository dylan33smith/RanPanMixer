#!/usr/bin/env bash
# PLINK LD blocks for the model-B chain unit, built from OUR staged 30x panel.
#
# ⚠ THE 10 TEST TARGETS ARE EXCLUDED. PanMixer's own block set was computed on all
# 3,202 panel samples, which CONTAINS our test targets -- so its block boundaries are
# weakly a function of the target, and standing constraint 1 says the segmentation of
# the output support must be a function of G and D only. We drop the 10, leaving 3,192
# (the same set the attack database uses, for the same reason).
#
# ⚠ RESIDUAL, and it is the same shape as deferred decision R4: a first-degree relative
# of a target still carries half its haplotypes, and we do not exclude those. For a real
# external target the panel the blocks come from must be DECLARED IN ADVANCE and must not
# contain the target or its relatives. Stated, not solved.
set -euo pipefail
cd "$(dirname "$0")/.."
PLINK=data/tools/plink
PANEL=data/raw/1000g_30x_phased.vcf.gz
OUT=data/cohort/chr21/blocks
TMP=data/temp/chr21/blocks_build
mkdir -p "$TMP" "$(dirname "$OUT")"

cut -f1 configs/targets.chr21.txt | awk 'NF{print $1"\t"$1}' > "$TMP/exclude.txt"
echo "excluding $(wc -l < "$TMP/exclude.txt") test targets from the LD panel"

$PLINK --vcf "$PANEL" --double-id --allow-extra-chr --snps-only just-acgt \
       --remove "$TMP/exclude.txt" --make-bed --out "$TMP/panel" --memory 64000
$PLINK --bfile "$TMP/panel" --blocks no-pheno-req --allow-extra-chr \
       --out "$OUT" --memory 64000
echo "blocks written to ${OUT}.blocks.det"
