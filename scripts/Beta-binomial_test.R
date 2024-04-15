#!/usr/bin/env R

#.libPaths("/cluster/home/bmoyers/R/x86_64-pc-linux-gnu-library/4.1")


library(GenomicRanges)
library(VGAM)

# read in all counts tables
files=as.list(readLines("/cluster/home/aanderson/allele_bias/Counts_WASP_noMAPQ0/path_counts.txt"))
n=sapply(strsplit(unlist(files),"[/]"),"[[",7)
n<-gsub(".txt","", n)

counts<-lapply(files, function(x){
  tab<-read.table(x,header=F)
})

# define NBT and BT
nbt <- function(a, b, p = 0.5) { 
  minor<-min(b, a)
  major<-max(b, a)
  format(pnbinom(minor, major, 0.5, lower.tail=T), scientific=F) }
bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
bbt<-function(a,b,rho,p=0.5){pbetabinom(a,b,0.5,rho)}

# loop through and get binom and negbinom p-values for every experiment
## only include snps that had at least 6 reads overlapping
Binom<-matrix(nrow=nrow(counts[[1]]), ncol=length(counts))
NegBinom<-matrix(nrow=nrow(counts[[1]]), ncol=length(counts))
BetaBinom<-matrix(nrow=nrow(counts[[1]]), ncol=length(counts))
rhoL<-c()
for(i in 1:length(counts)){
  tab<-counts[[i]]
  tab$V8<-tab$V6+tab$V7
  keep<-tab$V8 > 6
  #estimate rho
    tab2<-tab[keep,]
  tab2$y<-pmin(tab2$V6,tab2$V7)
  tab2$N<-tab2$V8
  tab2$mu<-0.5
  tab2$rho<-0.5
  fit <- vglm(cbind(y, N-y) ~ 1, betabinomial, data = tab2, trace = TRUE)
  rhoL<-c(rhoL,Coef(fit)[2])
  rho<-Coef(fit)[2]
  #perform tests
  tab$betabinom_p <- 1
  tab$minor<-pmin(tab$V6,tab$V7)
  tab$betabinom_p[keep] <- as.numeric(mapply(bbt, tab$minor[keep], tab$V8[keep], rho))
  BetaBinom[,i]<-tab$betabinom_p
  counts[[i]]<-tab
  print(paste0("Finished testing for:",n[i]))
}

saveRDS(counts, "~/allele_bias/Counts_WASP_noMAPQ0/counts_tmp.rds")
saveRDS(BetaBinom, "~/allele_bias/Counts_WASP_noMAPQ0/Beta-binom_tmp.rds")
df.rho<-data.frame(rho=rhoL, names=n)
write.csv(as.data.frame(rhoL), "~/allele_bias/Counts_WASP_noMAPQ0/rho-estimates.csv")




# only keep rows where at least 1 TF was tested
pp<-apply(BetaBinom, 1, min)
BetaBinom<-BetaBinom[pp<1,]

i=1
bb.gr<-GRanges(paste0(counts[[i]][which(pp<1),1],":",counts[[i]][which(pp<1),2]))

colnames(BetaBinom)<-n
mcols(bb.gr)<-BetaBinom

#write results
saveRDS(bb.gr, "~/allele_bias/beta-binomial_pvalues_grange.rds")

