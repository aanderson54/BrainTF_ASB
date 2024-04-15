#Building_vg_genomes_2022May01.sh



module load g/vg/1.20.0



#/gpfs/gpfs2/software/tabix/bgzip /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants_PASSED.vcf
#/gpfs/gpfs2/software/tabix/tabix -p vcf /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants_PASSED.vcf.gz


vg construct --threads 8 -a -r /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/ReferenceGenome/genome.fa -v /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants_PASSED.vcf.gz > /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.vg


vg prune -P /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.vg > /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.pruned.vg
vg index -x /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.xg /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.vg

bsub -n 1 -R rusage[mem=60000] -We 24:00 -q c7normal -o /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/gbwtIndexing.bsubLog.txt -J vgGBWT "module load g/vg/1.20.0; vg index --threads 8 -G /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.gbwt -v /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants.vcf.gz /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.vg"


bsub -n 1 -R rusage[mem=300000] -We 24:00 -q c7highmem -o /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/gscaIndexing.bsubLog.txt -J vgIndex "module load g/vg/1.20.0; mkdir ${HOME}/vgtemp && export TMPDIR=${HOME}/vgtemp; vg index -p -t 8 -x /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.xg -v /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants_PASSED.vcf.gz -G /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.gbwt -g /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.gcsa -k 11 /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.vg"








#/gpfs/gpfs2/software/tabix/bgzip -c /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf > /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.bgz
#/gpfs/gpfs2/software/tabix/tabix -p vcf /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.bgz


vg construct -r /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/ReferenceGenome/genome.fa -v /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants_PASSED.vcf.gz > /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.vg


vg prune -P /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.vg > /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.pruned.vg
vg index -x /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.xg /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.vg


bsub -n 1 -R rusage[mem=60000] -We 24:00 -q c7normal -o /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/gbwtIndexing.bsubLog.txt -J vgGBWT "module load g/vg/1.20.0; vg index --threads 8 -G /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.gbwt -v /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.gz /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.vg"


bsub -n 1 -R rusage[mem=300000] -We 24:00 -q c7highmem -o /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/gscaIndexing.bsubLog.txt -J vgIndex "module load g/vg/1.20.0; mkdir ${HOME}/vgtemp && export TMPDIR=${HOME}/vgtemp; vg index -X 3 -Z 4000 -p -t 8 -x /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.xg -v /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants_PASSED.vcf.gz -G /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.gbwt -g /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.gcsa -k 11 /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.vg"



#vg index -X 3 -Z 4000 -p -t 8 -x /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.xg -v /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.gz -G /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.gbwt -g /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.gcsa -k 11 /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/VariantGenome_Graph/hg38_with_variants.vg
