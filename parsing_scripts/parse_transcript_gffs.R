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


## For ease, we have outputted the relevant sessionInfo
## used when running this script 

#> sessionInfo()
#R version 4.2.0 (2022-04-22)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Catalina 10.15.7
#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#attached base packages:
#[1] stats4    stats     graphics  grDevices utils     datasets  methods  
#[8] base     
#other attached packages:
# [1] GenomicFeatures_1.48.0 AnnotationDbi_1.58.0   Biobase_2.56.0        
# [4] GenomicRanges_1.48.0   GenomeInfoDb_1.32.1    IRanges_2.30.0        
# [7] S4Vectors_0.34.0       BiocGenerics_0.42.0    forcats_0.5.1         
#[10] stringr_1.4.0          purrr_0.3.4            readr_2.1.2           
#[13] tidyr_1.2.0            tibble_3.1.7           ggplot2_3.3.6         
#[16] tidyverse_1.3.1        optparse_1.7.1         dplyr_1.0.9           
#[19] data.table_1.14.2     
#loaded via a namespace (and not attached):
#  [1] colorspace_2.0-3            rjson_0.2.21               
#  [3] ellipsis_0.3.2              rprojroot_2.0.3            
#  [5] XVector_0.36.0              fs_1.5.2                   
#  [7] rstudioapi_0.13             remotes_2.4.2              
#  [9] getopt_1.20.3               bit64_4.0.5                
# [11] fansi_1.0.3                 lubridate_1.8.0            
# [13] xml2_1.3.3                  cachem_1.0.6               
# [15] pkgload_1.2.4               jsonlite_1.8.0             
# [17] Rsamtools_2.12.0            broom_0.8.0                
# [19] dbplyr_2.1.1                png_0.1-7                  
# [21] BiocManager_1.30.17         compiler_4.2.0             
# [23] httr_1.4.3                  backports_1.4.1            
# [25] assertthat_0.2.1            Matrix_1.4-1               
# [27] fastmap_1.1.0               cli_3.3.0                  
# [29] bedtoolsr_2.30.0-4          prettyunits_1.1.1          
# [31] tools_4.2.0                 gtable_0.3.0               
# [33] glue_1.6.2                  GenomeInfoDbData_1.2.8     
# [35] rappdirs_0.3.3              Rcpp_1.0.8.3               
# [37] cellranger_1.1.0            vctrs_0.4.1                
# [39] Biostrings_2.64.0           rtracklayer_1.56.0         
# [41] ps_1.7.0                    brio_1.1.3                 
# [43] testthat_3.1.4              rvest_1.0.2                
# [45] lifecycle_1.0.1             restfulr_0.0.13            
# [47] devtools_2.4.3              XML_3.99-0.9               
# [49] zlibbioc_1.42.0             scales_1.2.0               
# [51] vroom_1.5.7                 hms_1.1.1                  
# [53] MatrixGenerics_1.8.0        parallel_4.2.0             
# [55] SummarizedExperiment_1.26.1 yaml_2.3.5                 
# [57] curl_4.3.2                  memoise_2.0.1              
# [59] biomaRt_2.52.0              stringi_1.7.6              
# [61] RSQLite_2.2.14              BiocIO_1.6.0               
# [63] desc_1.4.1                  filelock_1.0.2             
# [65] pkgbuild_1.3.1              BiocParallel_1.30.0        
# [67] rlang_1.0.2                 pkgconfig_2.0.3            
# [69] bitops_1.0-7                matrixStats_0.62.0         
# [71] lattice_0.20-45             GenomicAlignments_1.32.0   
# [73] bit_4.0.4                   tidyselect_1.1.2           
# [75] processx_3.5.3              magrittr_2.0.3             
# [77] R6_2.5.1                    generics_0.1.2             
# [79] DelayedArray_0.22.0         DBI_1.1.2                  
# [81] pillar_1.7.0                haven_2.5.0                
# [83] withr_2.5.0                 KEGGREST_1.36.0            
# [85] RCurl_1.98-1.6              modelr_0.1.8               
# [87] crayon_1.5.1                utf8_1.2.2                 
# [89] BiocFileCache_2.4.0         tzdb_0.3.0                 
# [91] progress_1.2.2              usethis_2.1.5              
# [93] grid_4.2.0                  readxl_1.4.0               
# [95] blob_1.2.3                  callr_3.7.0                
# [97] reprex_2.0.1                digest_0.6.29              
# [99] munsell_0.5.0               sessioninfo_1.2.2 

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