### This R Shiny script generates is used to deploy the app
### Last Updated - 05/20/22; Written by Manisha Munasinghe
### Note: The gggenomes package works best with the newest version of R


### This script is a copy of the one actually used
### Relevant data is stored on an AWS S3 server
### and must be queried using my credentials, 
### we have edited that out of this publicly available script

library(shiny)
library(gggenomes)
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(ggplot2)
library(DT)
library(IRanges)
library(aws.s3)
library(shinydashboard)
library(shinyFeedback)


## For ease, we have outputted the relevant sessionInfo
## used when running this script 

#> sessionInfo()
#R version 4.2.0 (2022-04-22)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 11.4#

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib#

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8#

#attached base packages:
#[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     #

#other attached packages:
# [1] shinyFeedback_0.4.0  IRanges_2.30.0       S4Vectors_0.34.0     BiocGenerics_0.42.0 
# [5] shiny_1.7.1          aws.s3_0.3.21        shinydashboard_0.7.2 DT_0.22             
# [9] data.table_1.14.2    gggenomes_0.9.5.9000 snakecase_0.11.0     jsonlite_1.8.0      
#[13] tibble_3.1.7         thacklr_0.0.0.9000   tidyr_1.2.0          stringr_1.4.0       
#[17] readr_2.1.2          purrr_0.3.4          gggenes_0.4.1        ggplot2_3.3.6       
#[21] dplyr_1.0.9          rsconnect_0.8.25     BiocManager_1.30.17 #

#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8.3        assertthat_0.2.1    digest_0.6.29       packrat_0.7.0       utf8_1.2.2         
# [6] aws.signature_0.6.0 mime_0.12           R6_2.5.1            httr_1.4.3          pillar_1.7.0       
#[11] rlang_1.0.2         curl_4.3.2          fontawesome_0.2.2   rstudioapi_0.13     jquerylib_0.1.4    
#[16] labeling_0.4.2      htmlwidgets_1.5.4   munsell_0.5.0       compiler_4.2.0      httpuv_1.6.5       
#[21] pkgconfig_2.0.3     askpass_1.1         base64enc_0.1-3     htmltools_0.5.2     openssl_2.0.0      
#[26] ggfittext_0.9.1     tidyselect_1.1.2    fansi_1.0.3         crayon_1.5.1        tzdb_0.3.0         
#[31] withr_2.5.0         later_1.3.0         grid_4.2.0          xtable_1.8-4        gtable_0.3.0       
#[36] lifecycle_1.0.1     DBI_1.1.2           magrittr_2.0.3      scales_1.2.0        cachem_1.0.6       
#[41] cli_3.3.0           stringi_1.7.6       farver_2.1.0        promises_1.2.0.1    bslib_0.3.1        
#[46] xml2_1.3.3          ellipsis_0.3.2      generics_0.1.2      vctrs_0.4.1         tools_4.2.0        
#[51] glue_1.6.2          crosstalk_1.2.0     hms_1.1.1           yaml_2.3.5          fastmap_1.1.0      
#[56] colorspace_2.0-3    sass_0.4.1  
#Need to run before deploying App
#library(BiocManager)
#library(rsconnect)
#options(repos = BiocManager::repositories())

NAM_lines <- c('B73','B97','Ky21','M162W','Ms71','Oh43','Oh7B','M37W','Mo18W','Tx303','HP301',
               'P39','Il14H','CML52','CML69','CML103','CML228','CML247','CML277','CML322','CML333',
               'Ki3','Ki11','NC350','NC358','Tzi8')

chr_positions <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")


## AnchorWave outputs do not cover the entirety of the chromosome
## notably telomeric regions at the end of the chromosome
## This function is used to warn users that the supplied coordinate inputs
## are outside of AW's information and consequently cannot be plotted
validate_coord_inputs <- function(chr_query_data,query_line,start_pos,end_pos) {
  query_sub <- chr_query_data %>% select(contains(query_line))
  max_end_point <- as.numeric(query_sub[nrow(query_sub),3])
  
  if (start_pos >= max_end_point | end_pos >= max_end_point) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

## Given user supplied endpoints, we often have to extend them so that the
## full region it falls within is visualized
## This can extend the visualized region substantially, but we felt
## it was more appropriate then truncating ends
obtain_alt_endpoints <- function(chr_query_data,query_line,start_pos,end_pos) {
    
    query_sub <- chr_query_data %>% select(contains(query_line))
    
    query_sub_start <- query_sub %>% subset(.[[2]] <= start_pos) %>%
        subset(start_pos <= .[[3]])

    query_start_coord <- as.numeric(query_sub_start[1,2])
    
    start_row <- chr_query_data[which(query_sub[,2]==query_start_coord),]
    
    query_sub_end <- query_sub %>% subset(.[[2]] <= end_pos) %>%
        subset(end_pos <= .[[3]])
    query_end_coord <- as.numeric(query_sub_end[1,2])
    end_row <- chr_query_data[which(query_sub[,2]==query_end_coord),]
    
    return(rbind(start_row,end_row))
    
}

ui <- dashboardPage(
    # Application title
  dashboardHeader(
    title="AnchorWave Identified Structural Variation in NAM",
    titleWidth = 450
  ),

  dashboardSidebar(
    selectInput('query','Reference:',choices=NAM_lines,selected='B73',multiple=F),
    textInput('coord','Coordinates:',value='chr7:151367073..151379756'),
    selectInput('comp','Comparison:',choices=NAM_lines,selected='B97',multiple=F),
    actionButton("go","Go")
  ),
  
  dashboardBody(
    useShinyFeedback(),
    fluidRow(
      shinydashboard::box(title='Visualize Nonvariant Regions',plotOutput("gggenomes_output",width='100%',height='200px'),width=12)
    ),
    fluidRow(
      shinydashboard::box(title='Summarized AnchorWave Output',DT::dataTableOutput(outputId = "table"),width=12)
    ),
    fluidRow(
      shinydashboard::box(title='Gene Tracks',DT::dataTableOutput(outputId = "gene_track"),width=12)
    ),
    fluidRow(
      shinydashboard::box(title='Present TEs',DT::dataTableOutput(outputId = "TE_track"),width=12)
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    ## Get summarized AnchorWave output given starting coordinates
    anchorwave_tsv <- eventReactive(input$go,{
          
        # Parse supplied inputs for all relevant details
        query_lineage <- input$query
        lineage_inputs <- c(input$query,input$comp)
        
        # Right now, this results in B73 always being the ID lineage
        # But in the future, if we have other direct comparisons
        # We would use this to find the relevant files
        ID_lineage <- sort(lineage_inputs)[1]
        ASM_lineage <- sort(lineage_inputs)[2]
        
        coord_inputs <- input$coord
        chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
        range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
        start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
        end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])

        # Ensure Start Coord < End Coord
        if (start_input > end_input) {
          showFeedbackWarning(
            inputId='coord',
            text='Start Coord < End Coord'
          )
          req(start_input < end_input)
        }
        
        # Ensure Supplied Chr is Valid
        if ( !(chr_input %in% chr_positions) ) {
          showFeedbackWarning(
            inputId='coord',
            text='Not Accepted Chromosome'
          )
          req(chr_input %in% chr_positions)
        }

        # Data is stored on an AWS S3 Server Through UMN - access it this way
        AnchorWave_dir <- 's3://shiny-namsv-data-storage/AnchorWave_Regions/'
        sub_dir <- paste(ID_lineage,ASM_lineage,'regions',sep='_')
        final_dir <- paste(AnchorWave_dir,sub_dir,'/',sep='')
        
        filename <- paste(chr_input,ID_lineage,ASM_lineage,'regions.tsv',sep='_')
        
        full_path <- paste(final_dir,filename,sep='')
        
        # Save file and read in data, and then delete file to save space
        chr_data <- save_object(full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
        unlink(filename)
        
        # Check Coordinates Within AW Range
        coord_within_AWrange <- validate_coord_inputs(chr_data,query_lineage,as.numeric(start_input),as.numeric(end_input))        
        if (!coord_within_AWrange) {
          showFeedbackWarning(
            inputId='coord',
            text='Coords Outside AW Alignment'
          )
          req(coord_within_AWrange)
        }

        # Find external end points for visualized region
        updated_end_points <- obtain_alt_endpoints(chr_data,query_lineage,as.numeric(start_input),as.numeric(end_input))
        

        # Add check if region results in one end point block
        # Update visualization endpoints
        if (nrow(updated_end_points)==1) {
          ID_start_pos <- as.numeric(updated_end_points[1,2])
          ID_end_pos <- as.numeric(updated_end_points[1,3])
          ASM_start_pos <- as.numeric(updated_end_points[1,5])
          ASM_end_pos <- as.numeric(updated_end_points[1,6]) 
        } else {
          ID_start_pos <- as.numeric(updated_end_points[1,2])
          ID_end_pos <- as.numeric(updated_end_points[2,3])
          ASM_start_pos <- as.numeric(updated_end_points[1,5])
          ASM_end_pos <- as.numeric(updated_end_points[2,6])  
        }
        
        # Subset AnchorWave output to region of interest
        selected_region <- chr_data %>% subset(.[[2]] >= ID_start_pos) %>% subset(.[[3]] <= ID_end_pos)
        
        # Pull out Structural Variant Regions
        # Everything else in that region is considered nonvariant
        # We fill in the dataframe accordingly 
        sv_selected_region <- selected_region %>% subset(Type=='structural_variant')

        if (nrow(sv_selected_region)==0) {
            final_region <- data.frame(chr_input,
                                       ID_start_pos,ID_end_pos,
                                       chr_input,
                                       ASM_start_pos,ASM_end_pos,
                                       'nonvariant_region',
                                       ID_start_pos-ID_start_pos,ID_end_pos-ID_start_pos,
                                       ASM_start_pos-ASM_start_pos,ASM_end_pos-ASM_start_pos)
            colnames(final_region) <- c(colnames(sv_selected_region),
                                        paste('ADJ',ID_lineage,'StartPos',sep='_'), 
                                        paste('ADJ',ID_lineage,'EndPos',sep='_'), 
                                        paste('ADJ',ASM_lineage,'StartPos',sep='_'),
                                        paste('ADJ',ASM_lineage,'EndPos',sep='_'))
            
        } else {
            if (ID_start_pos != as.numeric(sv_selected_region[1,2])) {
                
                id_chr <- as.numeric(updated_end_points[2,1])
                id_row_start <- ID_start_pos
                id_row_end <- (as.numeric(sv_selected_region[1,2])) - 1
                asm_chr <- as.numeric(updated_end_points[2,4])
                asm_row_start <- as.numeric(updated_end_points[1,5])
                asm_row_end <- (as.numeric(sv_selected_region[1,5]))-1
                
                sv_selected_region <- sv_selected_region %>% 
                    add_row(!!!setNames(list(id_chr,id_row_start,id_row_end,asm_chr,asm_row_start,asm_row_end,'nonvariant_region'), names(.))) %>% 
                    arrange(.[[1]],.[[2]])
                
            }
            
            if (ID_end_pos != as.numeric(sv_selected_region[nrow(sv_selected_region),3])) {

                
                id_chr <- as.numeric(updated_end_points[2,1])
                id_row_start <- (as.numeric(sv_selected_region[nrow(sv_selected_region),3])) + 1
                id_row_end <- ID_end_pos
                asm_chr <- as.numeric(updated_end_points[2,4])
                asm_row_start <- (as.numeric(sv_selected_region[nrow(sv_selected_region),6])) + 1
                asm_row_end <- ASM_end_pos	

                sv_selected_region <- sv_selected_region %>% 
                    add_row(!!!setNames(list(as.numeric(id_chr),id_row_start,id_row_end,as.numeric(asm_chr),asm_row_start,asm_row_end,'nonvariant_region'), names(.))) %>% 
                    arrange(.[[1]],.[[2]])

            }
            
            missing_rows <- data.frame(ID_Chr=as.numeric(),ID_StartPos=as.numeric(),
                                       ID_EndPos=as.numeric(),ASM_Chr=as.numeric(),ASM_StartPos=as.numeric(),
                                       ASM_EndPos=as.numeric(),Type=as.character())	
            
            for (i in seq(1,nrow(sv_selected_region)-1)) {
                junct_end <- as.numeric(sv_selected_region[i,3]) +1
                junct_start <- as.numeric(sv_selected_region[i+1,2])
                
                if (junct_end!=junct_start) {
                    id_chr <- as.numeric(updated_end_points[2,1])
                    id_row_start <- as.numeric(sv_selected_region[i,3]) + 1
                    id_row_end <- as.numeric(sv_selected_region[i+1,2]) -1
                    asm_chr <- as.numeric(updated_end_points[2,4])
                    asm_row_start <- as.numeric(sv_selected_region[i,6]) + 1
                    asm_row_end <- as.numeric(sv_selected_region[i+1,5]) -1
                    
                    missing_rows <- missing_rows %>% 
                        add_row(!!!setNames(list(id_chr,id_row_start,id_row_end,asm_chr,asm_row_start,asm_row_end,'nonvariant_region'), names(.))) %>% 
                        arrange(.[[1]],.[[2]])
                }
            }    
            
            colnames(missing_rows) <- colnames(sv_selected_region)
            
            final_region <- rbind(sv_selected_region,missing_rows) %>% arrange(.[[1]],.[[2]])
            
            final_region$ADJ_ID_StartPos <- final_region[,2] - ID_start_pos
            final_region$ADJ_ID_EndPos <- final_region[,3] - ID_start_pos
            final_region$ADJ_ASM_StartPos <- final_region[,5] - ASM_start_pos
            final_region$ADJ_ASM_EndPos <- final_region[,6] - ASM_start_pos
            
            colnames(final_region)[8] <- paste('ADJ',ID_lineage,'StartPos',sep='_')
            colnames(final_region)[9] <- paste('ADJ',ID_lineage,'EndPos',sep='_')
            colnames(final_region)[10] <- paste('ADJ',ASM_lineage,'StartPos',sep='_')
            colnames(final_region)[11] <- paste('ADJ',ASM_lineage,'EndPos',sep='_')
            
        }
        
        return(final_region)
    })
    
    #Obtain Seq Track From Full Region
    # For visualization, this is simply how long the region is for each line
    seq_track <- eventReactive(input$go,{
    
        ID_region_length <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3]) - as.numeric(anchorwave_tsv()[1,2])
        ASM_region_length <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6]) -as.numeric(anchorwave_tsv()[1,5])
        
        lineage_inputs <- c(input$query,input$comp)
        
        ID_lineage <- sort(lineage_inputs)[1]
        ASM_lineage <- sort(lineage_inputs)[2]
        
        seq_track <- tibble(
            seq_id= c(ID_lineage,ASM_lineage),
            length = c(ID_region_length,ASM_region_length)
        )
        
        return(seq_track)
        
    })
    
    # Obtain Link Track From Full Region
    # This is functionally just the nonvariant regions as they link the two lineages together
    link_track <- eventReactive(input$go,{
        nv_tracks <- anchorwave_tsv() %>% subset(Type=='nonvariant_region')

        # Add check if entire region is a structural variant
        if (nrow(nv_tracks)!=0) {
          lineage_inputs <- c(input$query,input$comp)
          
          ID_lineage <- sort(lineage_inputs)[1]
          ASM_lineage <- sort(lineage_inputs)[2]
          
          ID_adj <- nv_tracks %>% select(contains(ID_lineage)) %>% select(contains("ADJ"))
          ID_nv_tracks <- cbind(seq_id=ID_lineage,ID_adj)
          colnames(ID_nv_tracks) <- c("seq_id",'start','end')
          
          ASM_adj <- nv_tracks %>% select(contains(ASM_lineage)) %>% select(contains("ADJ"))
          ASM_nv_tracks <- cbind(seq_id=ASM_lineage,ASM_adj)
          colnames(ASM_nv_tracks) <- c("seq_id",'start','end')
          
          link_track <- cbind(ID_nv_tracks,ASM_nv_tracks)
          colnames(link_track) <- c('seq_id','start','end','seq_id2','start2','end2')
        } else {
          link_track <- NULL
        }

        return(link_track)
    })
    
    #Obtain Gene Track From Full Region
    # We now parse the gffs for relevant exon/intron structure
    gene_track <- eventReactive(input$go,{
        
        parent_dir <- 's3://shiny-namsv-data-storage/gff_gene_info/'
        
        lineage_inputs <- c(input$query,input$comp)
        
        ID_lineage <- sort(lineage_inputs)[1]
        ASM_lineage <- sort(lineage_inputs)[2]
        
        ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
        ID_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3])
        ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
        ASM_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6])
        
        coord_inputs <- input$coord
        chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
        range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
        start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
        end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
        
        ID_subdir <- paste(parent_dir,ID_lineage,'/',sep='')
        ID_file <- paste(paste(chr_input,ID_lineage,sep='_'),'.gff',sep='')
        ID_full_path <- paste(ID_subdir,ID_file,sep='')
        
        ID_gff_exons <- save_object(ID_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
        unlink(ID_file)
        
        sub_ID_gff_exons <- ID_gff_exons %>% subset(Start >=ID_start_pos) %>% subset(End <= ID_end_pos)
        
        
        ASM_subdir <- paste(parent_dir,ASM_lineage,'/',sep='')
        ASM_file <- paste(paste(chr_input,ASM_lineage,sep='_'),'.gff',sep='')
        ASM_full_path <- paste(ASM_subdir,ASM_file,sep='')        
        
        ASM_gff_exons <- save_object(ASM_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
        unlink(ASM_file)
        
        sub_ASM_gff_exons <- ASM_gff_exons %>% subset(Start >= ASM_start_pos) %>% subset(End <= ASM_end_pos)

        
        # We need to supply specific intron positions to the visualization package
        # So we essentially reverse engineer this from the exons in the gff
        if (nrow(sub_ID_gff_exons !=0 & nrow(sub_ASM_gff_exons) !=0)) {
            gene_track <- data.frame(
                ID_lineage,sub_ID_gff_exons$Start-ID_start_pos,sub_ID_gff_exons$End-ID_start_pos,
                sub_ID_gff_exons$Strand,'exon',sub_ID_gff_exons$Canonical_Tx
            )
        
        
            colnames(gene_track) <- c('seq_id','start','end','strand','type','name')
        
            gene_track <- gene_track %>% add_row(seq_id=ASM_lineage,
                                                 start=sub_ASM_gff_exons$Start-ASM_start_pos,
                                                 end=sub_ASM_gff_exons$End-ASM_start_pos,
                                                 strand=sub_ASM_gff_exons$Strand,
                                                 type='exon',
                                                 name=sub_ASM_gff_exons$Canonical_Tx)
        
        
            intron_gene_track <- tibble(
                seq_id=character(),
                start=numeric(),
                end=numeric(),
                strand=character(),
                introns=list()
            )
        
            for (gene in unique(gene_track$name)) {
                gene_sub <- gene_track %>% subset(name==gene)
                
                intron_positions <- c()
                
                if (nrow(gene_sub) ==1) {
                  intron_positions <- c()
                } else {
                  for (i in seq(1,nrow(gene_sub)-1)) {
                    intron_start <- (gene_sub[i,3]+1) - (as.numeric(gene_sub[1,2]))
                    intron_end  <- (gene_sub[i+1,2]-1) - (as.numeric(gene_sub[1,2]))
                    intron_positions <- c(intron_positions,intron_start,intron_end)
                  }                  
                }

                upd_df <- tibble(
                    seq_id=unique(gene_sub$seq_id),
                    start=as.numeric(gene_sub[1,2]),
                    end=as.numeric(gene_sub[nrow(gene_sub),3]),
                    strand=as.character(gene_sub[1,4]),
                    introns=list(intron_positions)
                )
                
                intron_gene_track <- rbind(intron_gene_track,upd_df)
            }
        } else {
            intron_gene_track <- NULL
        }
        
        return(intron_gene_track)
        
    })
    
    ## This function does the same as above, just stops before identifying introns for visualization
    ## We just output a cleaned up version of the dataframe
    output_gene_track <- eventReactive(input$go,{
        
      parent_dir <- 's3://shiny-namsv-data-storage/gff_gene_info/'
      
      lineage_inputs <- c(input$query,input$comp)
      
      ID_lineage <- sort(lineage_inputs)[1]
      ASM_lineage <- sort(lineage_inputs)[2]
      
      ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
      ID_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3])
      ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
      ASM_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6])

      coord_inputs <- input$coord
      chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
      range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
      start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
      end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
      
      ID_subdir <- paste(parent_dir,ID_lineage,'/',sep='')
      ID_file <- paste(paste(chr_input,ID_lineage,sep='_'),'.gff',sep='')
      ID_full_path <- paste(ID_subdir,ID_file,sep='')
      

      ID_gff_exons <- save_object(ID_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ID_file)
      
      sub_ID_gff_exons <- ID_gff_exons %>% subset(Start >=ID_start_pos) %>% subset(End <= ID_end_pos)
      
      
      ASM_subdir <- paste(parent_dir,ASM_lineage,'/',sep='')
      ASM_file <- paste(paste(chr_input,ASM_lineage,sep='_'),'.gff',sep='')
      ASM_full_path <- paste(ASM_subdir,ASM_file,sep='')        
      
      ASM_gff_exons <- save_object(ASM_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_file)
      
      sub_ASM_gff_exons <- ASM_gff_exons %>% subset(Start >= ASM_start_pos) %>% subset(End <= ASM_end_pos)
      
        
        if (nrow(sub_ID_gff_exons !=0 & nrow(sub_ASM_gff_exons) !=0)) {
            gene_track <- data.frame(
                ID_lineage,sub_ID_gff_exons$Start,sub_ID_gff_exons$End,
                sub_ID_gff_exons$Strand,'exon',sub_ID_gff_exons$Canonical_Tx
            )
            
            
            colnames(gene_track) <- c('seq_id','start','end','strand','type','name')
            
            gene_track <- gene_track %>% add_row(seq_id=ASM_lineage,
                                                 start=sub_ASM_gff_exons$Start,
                                                 end=sub_ASM_gff_exons$End,
                                                 strand=sub_ASM_gff_exons$Strand,
                                                 type='exon',
                                                 name=sub_ASM_gff_exons$Canonical_Tx)
        }
        
        colnames(gene_track) <- c("Lineage","Start",'Stop','Type','Transcript')
        return(gene_track)
    })
    
    ## This section parses the EDTA TE annotations to find
    ## fully encapsulated TEs in the region and then stores their coordinates for visualization
    TE_track <- eventReactive(input$go,{
      TE_anno_dir <- 's3://shiny-namsv-data-storage/TE_anno/'
      
      lineage_inputs <- c(input$query,input$comp)
      
      ID_lineage <- sort(lineage_inputs)[1]
      ASM_lineage <- sort(lineage_inputs)[2]
      
      ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
      ID_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3])
      ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
      ASM_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6])
      
      coord_inputs <- input$coord
      chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
      range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
      start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
      end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
      
      ID_TEAnno_subdir <- paste(TE_anno_dir,ID_lineage,'/',sep='')
      ID_TEAnno_file <- paste(paste(ID_lineage,chr_input,sep='_'),'_EDTA.TEanno.bed',sep='')
      ID_TEAnno_full_path <- paste(ID_TEAnno_subdir,ID_TEAnno_file,sep='')
      
      ID_TEAnno <- save_object(ID_TEAnno_full_path,base_url = 's3.msi.umn.edu', region = '',key='RU54HOMAC1II5NJ835UM',secret='FN8olrOP0DW51XqeKslASkLJpCHYWy3u4I8SJyNG')  %>% data.table::fread()
      unlink(ID_TEAnno_file)
      
      ID_TEAnno_subregion <- ID_TEAnno %>% subset(start >= ID_start_pos) %>%
        subset(end <= ID_end_pos)
      
      ID_present_TEs <- data.frame(
        seq_id=ID_lineage,
        start=ID_TEAnno_subregion$start - ID_start_pos,
        end=ID_TEAnno_subregion$end - ID_start_pos,
        type='transposable_element'
      )
      
      ASM_TEAnno_subdir <- paste(TE_anno_dir,ASM_lineage,'/',sep='')
      ASM_TEAnno_file <- paste(paste(ASM_lineage,chr_input,sep='_'),'_EDTA.TEanno.bed',sep='')
      ASM_TEAnno_full_path <- paste(ASM_TEAnno_subdir,ASM_TEAnno_file,sep='')
      
      ASM_TEAnno <- save_object(ASM_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_TEAnno_file)
      
      ASM_TEAnno_subregion <- ASM_TEAnno %>% subset(start >= ASM_start_pos) %>%
        subset(end <= ASM_end_pos)
      
      ASM_present_TEs <- data.frame(
        seq_id=ASM_lineage,
        start=ASM_TEAnno_subregion$start - ASM_start_pos,
        end=ASM_TEAnno_subregion$end - ASM_start_pos,
        type='transposable_element'
      )
      
      present_TEs <- rbind(ID_present_TEs,ASM_present_TEs)
      
      return(present_TEs)
      
    })
    
    # This function does the same as the above but outputs
    # a smaller, cleaner version of the data frame for outputting
    output_TE_track <- eventReactive(input$go,{
      TE_anno_dir <- 's3://shiny-namsv-data-storage/TE_anno/'
      
      
      lineage_inputs <- c(input$query,input$comp)
      
      ID_lineage <- sort(lineage_inputs)[1]
      ASM_lineage <- sort(lineage_inputs)[2]
      
      ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
      ID_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3])
      ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
      ASM_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6])

      coord_inputs <- input$coord
      chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
      range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
      start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
      end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
      
      ID_TEAnno_subdir <- paste(TE_anno_dir,ID_lineage,'/',sep='')
      ID_TEAnno_file <- paste(paste(ID_lineage,chr_input,sep='_'),'_EDTA.TEanno.bed',sep='')
      ID_TEAnno_full_path <- paste(ID_TEAnno_subdir,ID_TEAnno_file,sep='')
      
      ID_TEAnno <- save_object(ID_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ID_TEAnno_file)
      
      ID_TEAnno_subregion <- ID_TEAnno %>% subset(start >= ID_start_pos) %>%
        subset(end <= ID_end_pos)
      
      output_ID_TEAnno <- ID_TEAnno_subregion %>% 
        select(chr,start,end,strand,type,raw_superfamily,raw_family,method) %>% 
        separate(chr,c("lineage",'chr'))
      
      colnames(output_ID_TEAnno) <- c('Lineage','Chr','Start','End','Strand','Type','Superfamily','Family','Method')
      
      ASM_TEAnno_subdir <- paste(TE_anno_dir,ASM_lineage,'/',sep='')
      ASM_TEAnno_file <- paste(paste(ASM_lineage,chr_input,sep='_'),'_EDTA.TEanno.bed',sep='')
      ASM_TEAnno_full_path <- paste(ASM_TEAnno_subdir,ASM_TEAnno_file,sep='')
      
      ASM_TEAnno <- save_object(ASM_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_TEAnno_file)
      
      ASM_TEAnno_subregion <- ASM_TEAnno %>% subset(start >= ASM_start_pos) %>%
        subset(end <= ASM_end_pos)

      output_ASM_TEAnno <- ASM_TEAnno_subregion %>% 
        select(chr,start,end,strand,type,raw_superfamily,raw_family,method) %>% 
        separate(chr,c("lineage",'chr'))
      
      colnames(output_ASM_TEAnno) <- c('Lineage','Chr','Start','End','Strand','Type','Superfamily','Family','Method')
      
      output_TE_track <- rbind(output_ID_TEAnno, output_ASM_TEAnno)
      return(output_TE_track)
      
    })
    
    # This renders the actual visualization
    output$gggenomes_output <- renderPlot({
        
        lineage_inputs <- c(input$query,input$comp)
        
        ID_lineage <- sort(lineage_inputs)[1]
        ASM_lineage <- sort(lineage_inputs)[2]
        
        coord_inputs <- input$coord
        chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
        range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
        start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
        end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
        
        ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
        ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
        
        if (ID_lineage == input$query) {
            coord_start <- as.numeric(start_input) - ID_start_pos
            coord_end <- as.numeric(end_input) - ID_start_pos
        }
        if (ID_lineage == input$comp) {
            coord_start <- as.numeric(start_input) - ASM_start_pos
            coord_end <- as.numeric(end_input) - ASM_start_pos            
        }
        
        
        output_TE <- TE_track()

        if (is.null(gene_track()) & is.null(link_track())) {
          gggenomes(seqs=seq_track(),feats=list(output_TE)) +
            geom_seq(color='#742615',size=2) + 
            geom_feat(data=feats(output_TE),color='#F6C564') +
            geom_bin_label()+
            geom_vline(xintercept=coord_start,color='gray68',linetype='longdash') +
            geom_vline(xintercept=coord_end,color='gray68',linetype='longdash')
        } else if (is.null(gene_track())) {
          gggenomes(seqs=seq_track(),links=link_track(),feats=list(output_TE)) +
            geom_seq(color='#742615',size=2) + 
            geom_feat(data=feats(output_TE),color='#F6C564') +
            geom_bin_label() +
            geom_link() + 
            geom_vline(xintercept=coord_start,color='gray68',linetype='longdash') +
            geom_vline(xintercept=coord_end,color='gray68',linetype='longdash')
          
        } else if (is.null(link_track())) {
          gggenomes(seqs=seq_track(),genes=gene_track(),feats=list(output_TE)) +
            geom_seq(color='#742615',size=2) + 
            geom_feat(data=feats(output_TE),color='#F6C564') +
            geom_bin_label() +
            geom_gene(intron_shape=4.0,size=5,color='#06485E',fill='#06485E') +
            geom_vline(xintercept=coord_start,color='gray68',linetype='longdash') +
            geom_vline(xintercept=coord_end,color='gray68',linetype='longdash')
        } else {
          gggenomes(genes=gene_track(),seqs=seq_track(),links=link_track(),feats=list(output_TE)) +
            geom_seq(color='#742615',size=2) + 
            geom_feat(data=feats(output_TE),color='#F6C564') +
            geom_bin_label() + 
            geom_gene(intron_shape=4.0,size=5,color='#06485E',fill='#06485E') + 
            geom_link(color='#89B790',fill='#89B790') + 
            geom_vline(xintercept=coord_start,color='gray68',linetype='longdash') +
            geom_vline(xintercept=coord_end,color='gray68',linetype='longdash')
        }

    })
    
    # This outputs the summarized AnchorWave subregions (i.e., where is nonvariant/variant)
    output$table <- DT::renderDT(
        anchorwave_tsv()[,c(1,2,3,4,5,6,7)]
    )
    
    # This outputs the summarized gff exon structure for genes in the selected region
    output$gene_track <- DT::renderDT(
        
      if (!is.null(gene_track())) {
        output_gene_track()  
      }
      
    )
    
    # This outputs the summarized EDTA TE Annotation for TEs in the selected region 
    output$TE_track <- DT::renderDT(
      if (!is.null(output_TE_track())) {
        output_TE_track()
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
