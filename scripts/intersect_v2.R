#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)

ID=args[1]
EXP=args[2]

print(ID)

.libPaths("~/miniconda3/envs/R/lib/R/library")
library(rtracklayer)
library(stringr)


remap<-read.table(paste0("/cluster/home/aanderson/allele_bias/Mapping/vg/",ID,"/trial/",EXP,"_remapped.txt"), sep="\t")
orig<-read.table(paste0("/cluster/home/aanderson/allele_bias/Mapping/vg/",ID,"/trial/",EXP,"_original.txt"), sep="\t")

remap$p1<-sapply(strsplit(remap$V1, ".", fixed=T),`[`,2)
remap$p2<-sapply(strsplit(remap$V1, ".", fixed=T),`[`,3) #how many snps were flipped 
remap$p3<-sapply(strsplit(remap$V1, ".", fixed=T),`[`,4) #potential snps in read
remap<-remap[which(remap$p3<3 & remap$p2==remap$p3),] #p3 %in% c(3,4) is for variants that weren't phased so does iteration
remap<-remap[which(remap$V5>30),] #remove zero quality before proceding

remap$V1<-sapply(strsplit(remap$V1, ".", fixed=T),`[`,1)

mer<-merge( orig, remap,by="V1", all.x=T)  #keep all of orig in case it didn't align anywhere in map2
mer$agree<- ifelse(is.na(mer$V3.y)==T,F,ifelse(mer$V3.x==mer$V3.y & mer$V4.x==mer$V4.y, T,F))  #did it map to same position? Or did read not map at all?
mer2<-mer
mer2$seqnames<-mer2$V3.x
mer2$start<-mer2$V4.x
mer2$end<-mer2$V4.x
mer2.g<-GRanges(mer2)
mer2.g<-mer2.g[which(mer2.g$agree==F),]
red<-reduce(mer2.g, min.gapwidth=100)


export(resize(red, fix="center",width=100), paste0("~/allele_bias/Mapping/vg_noMAPQ0/DROP/",EXP,".bed"), format="bed")

write.table(unique(mer2.g$V1), paste0("~/allele_bias/Mapping/vg_noMAPQ0/DROP/",EXP,"_reads.txt"), row.names=F, col.names=F, quote=F)
print( paste0("~/allele_bias/Mapping/vg/DROP/",ID,"_",EXP,"_reads.txt"))

#remove reads that mapped incorrectly
hets<-read.table(paste0("~/allele_bias/Mapping/vg/HET_reads/",EXP,"_het.txt"))
hets2<-hets[!(hets$V1 %in% mer2.g$V1),]
write.table(hets2, paste0("~/allele_bias/Mapping/vg_noMAPQ0/HET_reads/",EXP,"_het_filtered.txt"), row.names=F, col.names=F, quote=F )

