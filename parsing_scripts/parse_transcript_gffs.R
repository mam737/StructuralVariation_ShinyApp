### This RScript Simply Parses The Gene gff  outputs found on maizegdb.org
### It also leverages the canonical gene list for each line found there as well.
### We do this to make reading this output into our shiny app easier
### Last Updated - 05/12/22; Written by Andy Read


library("optparse")
library("tidyverse")
library("stringr")
library("data.table")
library("GenomicFeatures")
library("bedtoolsr")

setwd("/Users/read0094/Desktop/Maize/FunTest")

##### Step1_create_primary_transcript_gff #####

## User needs to supply file1 and file2 manually
## We've demonstrated what that looks like below
file1="Zm-B97-REFERENCE-NAM-1.0_Zm00018ab.1.gff3"
file2="Zm-B97-REFERENCE-NAM-1.0_Zm00018ab.1.canonical_transcripts"

CanonicalGFF <- function(file1,file2){
  GFF = fread(file1, header=F, sep="\t",skip=1)
  
  Canonical = fread(file2, header=F, sep="\t")
  
  Canonical$V2 = "canonical"
  
  GFF$V10 = GFF$V9
  
  GFF1 = gsub(';.*', '', GFF$V10)
  GFF1 = as.data.frame(GFF1)
  
  GFF1 = gsub('.*=', '', GFF1$GFF1)
  GFF1 = as.data.frame(GFF1)
  
  GFF$V10 = GFF1$GFF1
  
  Canonical2 = dplyr::rename(Canonical, V10=V1)
  
  join = left_join(GFF,Canonical2, by="V10")
  
  join1 = subset(join, V3=="chromosome" | V3=="gene" | V2.y=="canonical")
  

  setwd("/Users/read0094/Desktop/Maize/FunTest")
  out_file=paste("canonical_",file1,sep="")
  write_tsv(join1, out_file, col_names=F)
  
  #make a gene orientation file
  orientation=dplyr::select(join1, 10,7)
  out_file2=paste("orientation_",file1,sep="")
  write_tsv(orientation, out_file2, col_names=F)
}