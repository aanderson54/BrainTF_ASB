#!/bin/bash

module load cluster/R/4.2.1

file=$1

ID=$(echo $file | cut -d "/" -f 6)
if [ $ID == "5397-JL-0001" ]; then
   rep="1224"
else
   rep="1230"
fi


tmp=$(echo $file | cut -d "/" -f 8)
EXP=$(echo $tmp | cut -d "." -f 1)



# START

## step3: Intersect SNPs

mkdir -p ${rep}/intersect/

python ~/software/WASP/mapping/find_intersecting_snps.py \
      --is_sorted \
      --output_dir ${rep}/intersect/${EXP} \
      --snp_tab ${rep}/het_${rep}_snp_tab.h5 \
      --snp_index ${rep}/het_${rep}_snp_index.h5 \
      --haplotype ${rep}/het_${rep}_haplotypes.h5 \
      --samples $ID \
      --max_snps 2 \
      ${file}

## step3.5: remove space in read name to filter later
gunzip ${rep}/intersect/${EXP}/${EXP}.filt.nodup.srt.remap.fq.gz
sed 's@ @@g' ${rep}/intersect/${EXP}/${EXP}.filt.nodup.srt.remap.fq | gzip > ${rep}/intersect/${EXP}/${EXP}.filt.nodup.srt.remap.fq.gz

## step4: Reallign
###### Must vg command
mkdir -p ${rep}/map2

ALN_GAM=${rep}/map2/${EXP}.raw.gam
ALN_BAM=${rep}/map2/${EXP}.raw.bam
vg map -t 1 -A -K -M 3 -f ${rep}/intersect/${EXP}/${EXP}.filt.nodup.srt.remap.fq.gz -x VariantGenome_Graph/hg38_with_variants.xg -g VariantGenome_Graph/hg38_with_variants.gcsa -1 VariantGenome_Graph/hg38_with_variants.gbwt > $ALN_GAM

vg surject -t 1 -b -x VariantGenome_Graph/hg38_with_variants.xg -b $ALN_GAM > $ALN_BAM


samtools sort -o ${rep}/map2/${EXP}.sort.bam $ALN_BAM 
samtools index ${rep}/map2/${EXP}.sort.bam

# write to sam to determine mapping position

samtools view  ${rep}/map2/${EXP}.sort.bam >  ${rep}/trial/${EXP}_remapped.txt
samtools view ${rep}/intersect/${EXP}/${EXP}.filt.nodup.srt.to.remap.bam >  ${rep}/trial/${EXP}_original.txt
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $10}' ${rep}/trial/${EXP}_original.txt > tmp_${EXP}
mv tmp_${EXP} ${rep}/trial/${EXP}_original.txt 

#get het reads
awk '{print $1}' ${ID}/MappedReads_vg/${EXP}.filt.nodup.srt.hetsOnly.sam > HET_reads/${EXP}_het.txt

# Only output het reads with no mapping issues
Rscript ./intersect_v2.R ${rep} ${EXP}

samtools view ${file} | grep -f HET_reads/${EXP}_het_filtered.txt > HET_reads/${EXP}_het_filtered.sam



