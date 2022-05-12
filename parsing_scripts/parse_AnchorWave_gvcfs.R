### This RScript Simply Parses The GVCF  output from AnchorWave
### We do this to make reading this output into our shiny app easier
### Last Updated - 05/12/22; Written by Manisha Munasinghe

.libPaths('/home/brandvai/mmunasin/Rlibs4')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)

## For ease, we have outputted the relevant sessionInfo
## used when running this script 

#> sessionInfo()
#R version 4.0.4 (2021-02-15)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)#

#Matrix products: default
#BLAS:   /panfs/roc/msisoft/R/4.0.4/lib64/R/lib/libRblas.so
#LAPACK: /panfs/roc/msisoft/R/4.0.4/lib64/R/lib/libRlapack.so#

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C#

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base#

#other attached packages:
#[1] tidyr_1.2.0       data.table_1.14.2 stringr_1.4.0     dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] fansi_1.0.3      assertthat_0.2.1 utf8_1.2.2       crayon_1.5.1
# [5] R6_2.5.1         DBI_1.1.2        lifecycle_1.0.1  magrittr_2.0.3
# [9] pillar_1.7.0     stringi_1.7.6    rlang_1.0.2      cli_3.3.0
#[13] vctrs_0.4.1      generics_0.1.2   ellipsis_0.3.2   tools_4.0.4
#[17] glue_1.6.2       purrr_0.3.4      compiler_4.0.4   pkgconfig_2.0.3
#[21] tidyselect_1.1.2 tibble_3.1.6

## This function simply removes extraneous text from columns
remove_extra_txt <- function(col_string) {
	rel_text <- unlist(str_split(col_string,pattern='='))[2]
	return(rel_text)
}

## This function looks at variant regions in the gvcf and classifies them as either a 
## SNP (a single base pair differences)
## InDel (a structural variant less than 100 bp long)
## Structural Variant (a variant larger than 100 bp long)
determine_variant_type <- function(row.df) {

	ID_size <- (as.numeric(row.df[3]) - as.numeric(row.df[2]))+1
	ASM_size <- (as.numeric(row.df[6]) - as.numeric(row.df[5])) +1

	size_dif <- abs(ID_size - ASM_size)

	if (size_dif == 0 ) {
		return('SNP')
	}
	if (size_dif > 0 & size_dif < 100) {
		return('InDel')
	} 
	if (size_dif >= 100) {
		return('structural_variant')
	}

}

## We supply the gvcf file as an argument to the RScript
## Prior to the Rscript we remove the header from the gvcf
args <- commandArgs(trailingOnly = TRUE)
gvcf_filename <- args[1]

## For our supplied AnchorWave gvcfs, we know that the ID lineage is always B73
## The ASM lineage varies depending on the file but can be identified in the filename
ID_lineage <- 'B73'
ASM_lineage <- unlist(str_split(unlist(str_split(gvcf_filename,pattern='/'))[8],pattern='_'))[1]

## Read in the gvcf and adjust any colnames as necessary
gvcf <- fread(gvcf_filename)
names(gvcf)[1] <- 'CHROM'

## We start by first pulling in all nonvariant regions in the gvcf
## We identify these by the presence of the END tag in the Info column
nonvariant <- gvcf[grep("END",gvcf$INFO),]

## We separate the INFO column as it contains relevant information
nonvariant <- nonvariant %>% tidyr::separate(col=INFO,into=c('ASM_CHR','ASM_END','ASM_Start','ASM_Strand', 'END'),sep=';')

## We remove any extraneous text from these columns
nonvariant$ASM_CHR <- unlist(lapply(nonvariant$ASM_CHR,remove_extra_txt))
nonvariant$ASM_END <- unlist(lapply(nonvariant$ASM_END,remove_extra_txt))
nonvariant$ASM_Start <- unlist(lapply(nonvariant$ASM_Start,remove_extra_txt))
nonvariant$ASM_Strand <- unlist(lapply(nonvariant$ASM_Strand,remove_extra_txt))
nonvariant$END <- unlist(lapply(nonvariant$END,remove_extra_txt))

## Bind relevant colummns for nonvariant regions
nv_regions <- data.frame(nonvariant$CHROM,nonvariant$POS,nonvariant$END,nonvariant$ASM_CHR,nonvariant$ASM_Start,nonvariant$ASM_END)

## We then update all of the column names in a standard order which makes parsing later easier
col1 <- paste(ID_lineage,'Chr',sep='_')
col2 <- paste(ID_lineage,'StartPos',sep='_')
col3 <- paste(ID_lineage,'EndPos',sep='_')
col4 <- paste(ASM_lineage,'Chr',sep='_')
col5 <- paste(ASM_lineage,'StartPos',sep='_')
col6 <- paste(ASM_lineage,'EndPos',sep='_')


colnames(nv_regions) <- c(col1,col2,col3,col4,col5,col6)
nv_regions$Type <- 'nonvariant_region'

## We now do the same but for all identified variant regions in the gvcf
## We identify these by the absence of the END tag in the Info column
## We then need to parse and identify what kind of variant it is
variant <- gvcf[!grep("END",gvcf$INFO),]
variant$ref_bp_size <- unlist(lapply(variant$REF,nchar))
variant$REF_END <- (variant$POS + variant$ref_bp_size)-1

## We separate the INFO column as it contains relevant information
variant <- variant %>% tidyr::separate(col=INFO,into=c('ASM_CHR','ASM_END','ASM_Start','ASM_Strand'),sep=';')

## We remove any extraneous text from these columns
variant$ASM_CHR <- unlist(lapply(variant$ASM_CHR,remove_extra_txt))
variant$ASM_END <- unlist(lapply(variant$ASM_END,remove_extra_txt))
variant$ASM_Start <- unlist(lapply(variant$ASM_Start,remove_extra_txt))
variant$ASM_Strand <- unlist(lapply(variant$ASM_Strand,remove_extra_txt))

## Bind relevant colummns for variant regions
v_regions <- data.frame(variant$CHROM,variant$POS,variant$REF_END, variant$ASM_CHR,variant$ASM_Start, variant$ASM_END)

## We then update all of the column names in a standard order which makes parsing later easier
colnames(v_regions) <- c(col1,col2,col3,col4,col5,col6)

## Here, we identify what kind of variant type it is and bind that
v_regions$Type <- apply(v_regions,1,determine_variant_type)

## We then combine our nonvariant and variant regions
## and reorder so that it is in the proper position order
full_region <- rbind(nv_regions,v_regions) %>% arrange(.[[1]],.[[2]])

## We then create a directory that stores all of these parsed gvcf outputs
## and split the whole gvcf into chromosomal gvcfs to minimize the amount
## of data we have to read into the app
output_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/AnchorWave_Regions/'
output_subdir <- paste(ID_lineage,ASM_lineage,'regions',sep='_')

final_output_dir <- paste(output_dir,output_subdir,sep='')

dir.create(final_output_dir)

for (i in unique(full_region[[1]])) {
	sub_full_region <- full_region %>% subset(.[[1]]==i)
	chr_header <- paste('chr',i,sep='')
	filename <- paste(chr_header,ID_lineage,ASM_lineage,'regions.tsv',sep='_')

	output_path <- paste(final_output_dir,filename,sep='/')
	write.table(sub_full_region,output_path,quote=F,col.names=T,row.names=F)
	
}

