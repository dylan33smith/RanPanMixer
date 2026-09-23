set -euo pipefail
E=/home/ds85/miniconda3/envs/panmixer
export PATH=$E/bin:$PATH
P=/data/ds85/RanPanMixer/external/PanMixer
REFDIR=$P/starting_data/references
mkdir -p $REFDIR $P/read_fastqs
cd $REFDIR

echo "### [1/4] samtools $(date +%H:%M:%S)"
command -v samtools >/dev/null || conda install -y -q -p $E -c bioconda -c conda-forge samtools >/dev/null
samtools --version | head -1

echo "### [2/4] GRCh38 analysis-set reference (needed to DECODE the CRAMs) $(date +%H:%M:%S)"
REF=GRCh38_full_analysis_set_plus_decoy_hla.fa
if [ ! -s $REF.fai ]; then
  wget -nv -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/$REF" -O $REF
  samtools faidx $REF
fi
ls -la $REF $REF.fai

echo "### [3/4] chr21.fa for vg autoindex (sliced from the same reference) $(date +%H:%M:%S)"
if [ ! -s chr21.fa.fai ]; then samtools faidx $REF chr21 > chr21.fa; samtools faidx chr21.fa; fi
head -c 60 chr21.fa; echo; grep -c "" chr21.fa.fai

echo "### [4/4] chr21 read slices for the paper's five external donors $(date +%H:%M:%S)"
# Paper (Methods, read mapping): HG00138 EUR, HG00635 EAS, HG01112 AMR, HG02698 SAS, NA18853 AFR.
# NOTE constants.py READ_SUBJECTS substitutes HG01600 for NA18853 -- we follow the PAPER.
while read -r S U; do
  HTTPS=$(echo "$U" | sed 's#ftp://ftp.sra.ebi.ac.uk#https://ftp.sra.ebi.ac.uk#')
  OUT=$P/read_fastqs/${S}
  if [ -s ${OUT}_1.fastq.gz ]; then echo "  $S already present"; continue; fi
  echo "  --- $S  $(date +%H:%M:%S)"
  samtools view -T $REFDIR/$REF -b -o ${OUT}.chr21.bam "$HTTPS" chr21
  samtools collate -u -O ${OUT}.chr21.bam | samtools fastq -1 ${OUT}_1.fastq.gz -2 ${OUT}_2.fastq.gz -0 /dev/null -s /dev/null -n
  rm -f ${OUT}.chr21.bam
  ls -la ${OUT}_1.fastq.gz ${OUT}_2.fastq.gz
done <<'SAMPLES'
HG00138 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240144/HG00138.final.cram
HG00635 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3988806/HG00635.final.cram
HG01112 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241828/HG01112.final.cram
HG02698 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989056/HG02698.final.cram
NA18853 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239382/NA18853.final.cram
SAMPLES
echo "### DONE $(date +%H:%M:%S)"; ls -la $P/read_fastqs/
