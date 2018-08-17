
##################################################
  ## Script purpose: Make a distance matrix from an Identity By State analysis of multiple VCFs. This allows for generation of a network graph comparing sample to sample based on reference genomw alignment.
  ## Date: January 2018
  ## Author: Tim Regan
  ##################################################
#Firstly, the VCF needs to be reduced to chromosomes (autosomames in this case) which need to be named by number.
#We called this "variants.vqsr.chr_renamed.vcf"
#The resulting matrix can be used to construct a network graph.
#setwd() to working directory
#install.packages("SNPRelate")
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")

library(SNPRelate)

vcf.fn<- readline("What and where is your variants.vqsr.chr_renamed.vcf please?")

snpgdsVCF2GDS(vcf.fn, "ccm.gds", method="biallelic.only")

genofile <- snpgdsOpen("ccm.gds")
#ccm_pca<-snpgdsPCA(genofile)

samples_used <- readline("Which samples are you interested in? Use quotes for each and separate with a comma.")

ibs_objects <- snpgdsIBS(genofile, sample.id=c(samples_used), snp.id = NULL, autosome.only = TRUE, remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, num.thread = 1, verbose = TRUE)

ibs<-ibs_objects$ibs
names<-bee_matrix$sample.id

matrixfile <- readline("What would you like your matrix file to be called?")

write.table(newfile, file = matrixfile, sep = "\t", row.names=names, col.names=NA)