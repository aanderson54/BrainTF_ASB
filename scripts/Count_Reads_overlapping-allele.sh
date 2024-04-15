#!/bin/bash

sam=$1
pre=$(echo $sam | cut -d "/" -f 9)
tis=$(echo $pre | cut -d "_" -f 1)
tf=$(echo $pre | cut -d "_" -f 2)
rep=$(echo $pre | cut -d "_" -f 3)

EXP=${tis}_${tf}_${rep}

if [ $rep == "1224" ]; then
  ID="5397-JL-0001"
else
  ID="5397-JL-0002"
fi

cat HEADER.tmp $sam > /cluster/home/aanderson/allele_bias/Counts_WASP_noMAPQ0/${EXP}.sam
samtools view -bS ${EXP}.sam > ${EXP}.tmp.bam
samtools sort ${EXP}.tmp.bam > ${EXP}.srt.bam
samtools index ${EXP}.srt.bam

rm ${EXP}.sam
rm ${EXP}.tmp.bam

python ~/software/WASP/CHT/bam2h5.py --chrom ../chromInfo.hg38.txt \
      --snp_index /cluster/home/aanderson/allele_bias/Mapping/vg/${rep}/het_${rep}_snp_index.h5 \
      --snp_tab /cluster/home/aanderson/allele_bias/Mapping/vg/${rep}/het_${rep}_snp_tab.h5 \
      --haplotype /cluster/home/aanderson/allele_bias/Mapping/vg/${rep}/het_${rep}_haplotypes.h5 \
      --individual $ID \
      --ref_as_counts H5/${EXP}_ref.h5 \
      --alt_as_counts H5/${EXP}_alt.h5 \
      --other_as_counts H5/${EXP}_other.h5 \
      --read_counts H5/${EXP}_reads.h5 \
      --txt_counts ${EXP}.txt \
     ${EXP}.srt.bam

