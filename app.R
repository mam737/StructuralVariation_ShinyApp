### This R Shiny script generates is used to deploy the app
### Last Updated - 10/12/22; Written by Manisha Munasinghe
### Note: The gggenomes package works best with R 4.2


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



## This function is used to warn users if the supplied coordinate inputs
## are outside of AW's outputted information and consequently cannot be plotted
validate_coord_inputs <- function(chr_query_data,query_line,start_pos,end_pos) {
  query_sub <- chr_query_data %>% select(contains(query_line))
  max_end_point <- as.numeric(query_sub[nrow(query_sub),3])
  
  if (start_pos >= max_end_point | end_pos >= max_end_point) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#Check if Supplied Endpoints Fall In Unalignable Region
#This is required to crop visualizations at user supplied endpoints
no_unalignable_endpoints <- function(end_points.df) {
  if ('unalignable' %in% end_points.df$Reformat_Type) {
    if (end_points.df[[1,7]] =='unalignable') {
      return('Start Unalignable')
    } else {
      return('End Unalignable')
    }
  } else {
    return(FALSE)
  }
}

#Provide alternative coordinates if supplied endpoints
#are in unalignable regions
provide_alignable_endpoints <- function(all_data,end_points.df,type) {
  if (type=='Start Unalignable') {
    ind <- which(all_data$AW_Block== end_points.df[[1,8]])
    potential <- all_data[[ind-1,2]]
  } else {
    ind <- which(all_data$AW_Block== end_points.df[[nrow(end_points.df),8]])
    potential <- all_data[[ind+1,3]]
  }
  return(potential)
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

#This function adjusts the relevant AnchorWave region to adjust the returned
#starting row that contains the user supplied start coordinate
#with the supplied start coordinate
return_original_start_endpoints <- function(end_points.df,query_line,start_input) {
  if (query_line == 'B73') {
    # Check to see if start_input equals the start of the
    # starting AW Block -- if not, we need to adjust
    if (start_input != end_points.df[[1,2]]) {
      # Pull out the starting row
      start_row <- end_points.df %>% dplyr::slice(1)
      # If the starting row is alignable
      if (start_row[[1,7]] == "alignable_region") {
        # determine how offset the supplied starting point
        # is from the start of the AW starting block start coord
        start_offset <- start_input - end_points.df[[1,2]]
        # Update AW starting block start coord to start input
        start_row[[1,2]] <- start_input
        # Evaluate what the offset start coord in comparison line
        # will be
        asm_new_start <- start_row[[1,5]]+start_offset
        # If comparison offset start coord is smaller than the 
        # comparison end coord, change comparison start to comparison
        # offset start
        if (asm_new_start < start_row[[1,6]]) {
          start_row[[1,5]] <- asm_new_start
        } else { # If not, take comparison end coord and subtract 10 bp
          start_row[[1,5]] <- start_row[[1,6]] - 10
        } # This can be caused by differences in small InDels
      } else if (start_row[[1,7]] == 'structural_insertion_inB73') {
        # determine how offset the supplied starting point
        # is from the start of the AW starting block start coord
        start_offset <- start_input - end_points.df[[1,2]]
        start_row[[1,2]] <- start_input
        # We do not have to updated the asm start
        # because it is a single base pair
      } else {
        # if the start coord != block position
        # it cannot be a structural insertion in NAM, as that should
        # only be a single bp, that means only remaining option is unalignable
        start_row <- end_points.df %>% dplyr::slice(1) 
        # We are not making any adjustments, since
        # we will flag this later to not visualize
      }
    } else {
      start_row <- end_points.df %>% dplyr::slice(1) 
    } 
  } else {
    # Check to see if start_input equals the start of the
    # starting AW Block -- if not, we need to adjust
    if (start_input != end_points.df[[1,5]]) {
      # Pull out the starting row
      start_row <- end_points.df %>% dplyr::slice(1)
      # If the starting row is alignable
      if (start_row[[1,7]] == "alignable_region") {
        # determine how offset the supplied starting point
        # is from the start of the AW starting block start coord
        start_offset <- start_input - end_points.df[[1,5]]
        # Update AW starting block start coord to start input
        start_row[[1,5]] <- start_input
        # Evaluate what the offset start coord in comparison line
        # will be
        id_new_start <- start_row[[1,2]]+start_offset
        # If comparison offset start coord is smaller than the 
        # comparison end coord, change comparison start to comparison
        # offset start  
        if (asm_new_start < start_row[[1,3]]) {
          start_row[[1,2]] <- asm_new_start
        } else { # If not, take comparison end coord and subtract 10 bp
          start_row[[1,2]] <- start_row[[1,3]] - 10
        } # This can be caused by differences in small InDels
      } else if (!(start_row[[1,7]] %in% c("structural_insertion_inB73","unalignable","alignable_region"))) {
        # determine how offset the supplied starting point
        # is from the start of the AW starting block start coord
        start_offset <- start_input - end_points.df[[1,5]]
        start_row[[1,5]] <- start_input
        # We do not have to updated the asm start
        # because it is a single base pair
      } else {
        # if the start coord != block position
        # it cannot be a structural insertion in B73, as that should
        # only be a single bp, that means only remaining option is unalignable
        start_row <- end_points.df %>% dplyr::slice(1) 
        # We are not making any adjustments, since
        # we will flag this later to not visualize
      }
    } else {
      start_row <- end_points.df %>% dplyr::slice(1) 
    }
  }
  return(start_row)
}

#This function adjusts the relevant AnchorWave region to adjust the returned
#ending row that contains the user supplied end coordinate
#with the supplied end coordinate
return_original_end_endpoints <- function(end_points.df,query_line,end_input) {
  if (query_line == 'B73') {
    # Check to see if end_input equals the end of the
    # end AW Block -- if not, we need to adjust
    if (end_input != end_points.df[[nrow(end_points.df),3]]) {
      # Pull out the end row
      end_row <- end_points.df %>% dplyr::slice(n())
      # If the ending row is alignable
      if (end_row[[1,7]] == "alignable_region") {
        # determine how offset the supplied end point
        # is from the end of the AW ending block start coord
        end_offset <- end_points.df[[nrow(end_points.df),3]] - end_input
        # Update AW end block end coord to end input
        end_row[[1,3]] <- end_input
        # Evaluate what the offset end coord in comparison line
        # will be
        asm_new_end <- end_row[[1,6]]-end_offset
        # If comparison offset end coord is larger than the 
        # comparison start coord, change comparison end to comparison
        # offset end
        if (asm_new_end > end_row[[1,5]]) {
          end_row[[1,6]] <- asm_new_end
        } else { # If not, take comparison end coord and add 10 bp
          end_row[[1,6]] <- end_row[[1,5]] + 10
        } # This can be caused by differences in small InDels
      } else if (end_row[[1,7]] == 'structural_insertion_inB73') {
        # determine how offset the supplied starting point
        # is from the start of the AW starting block start coord
        end_offset <- end_points.df[[nrow(end_points.df),3]] - end_input
        end_row[[1,3]] <- end_input
        # We do not have to updated the asm start
        # because it is a single base pair
      } else {
        # if the start coord != block position
        # it cannot be a structural insertion in NAM, as that should
        # only be a single bp, that means only remaining option is unalignable
        end_row <- end_points.df %>% dplyr::slice(n()) 
        # We are not making any adjustments, since
        # we will flag this later to not visualize
      }
    } else {
      end_row <- end_points.df %>% dplyr::slice(n()) 
    }   
  } else {
    # Check to see if end_input equals the end of the
    # ending AW Block -- if not, we need to adjust
    if (end_input != end_points.df[[1,6]]) {
      # Pull out the starting row
      end_row <- end_points.df %>% dplyr::slice(n())
      # If the ending row is alignable
      if (end_row[[1,7]] == "alignable_region") {
        # determine how offset the supplied end point
        # is from the start of the AW ending block end coord
        end_offset <- end_points.df[[nrow(end_points.df),6]] - end_input
        # Update AW ending block start coord to start input
        end_row[[1,6]] <- end_input
        # Evaluate what the offset start coord in comparison line
        # will be
        id_new_end <- end_row[[1,3]]-end_offset
        # If comparison offset start coord is smaller than the 
        # comparison end coord, change comparison start to comparison
        # offset start  
        if (id_new_end > end_row[[1,2]]) {
          end_row[[1,3]] <- id_new_end
        } else { # If not, take comparison end coord and subtract 10 bp
          end_row[[1,3]] <- end_row[[1,2]] + 10
        } # This can be caused by differences in small InDels
      } else if (!(end_row[[1,7]] %in% c("structural_insertion_inB73","unalignable","alignable_region"))) {
        # determine how offset the supplied starting point
        # is from the start of the AW starting block start coord
        end_offset <- end_points.df[[nrow(end_points.df),6]] - end_input
        end_row[[1,6]] <- end_input
        # We do not have to updated the asm start
        # because it is a single base pair
      } else {
        # if the start coord != block position
        # it cannot be a structural insertion in B73, as that should
        # only be a single bp, that means only remaining option is unalignable
        end_row <- end_points.df %>% dplyr::slice(n()) 
        # We are not making any adjustments, since
        # we will flag this later to not visualize
      }
    } else {
      end_row <- end_points.df %>% dplyr::slice(n()) 
    }
  }
  return(end_row)
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
      shinydashboard::box(title='Warning Messages',verbatimTextOutput("feature_warning_message"),width=12)
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
        
        filename <- paste(chr_input,ID_lineage,ASM_lineage,'summarised_regions.tsv',sep='_')
        
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

        # Check and Output Error Message if Supplied Endpoints
        # Fall in an Unalignable Rgion
        unalignable_endpoints <- no_unalignable_endpoints(updated_end_points)
        if (unalignable_endpoints != 'FALSE') {
          if (unalignable_endpoints=='Start Unalignable') {
            alt_start_coord <- provide_alignable_endpoints(chr_data,updated_end_points,unalignable_endpoints)
            out_text <- paste('Start Coordinate Falls In Unalignable Region. Start Coordinate of Upstream Region - ',alt_start_coord,sep='')
          } else {
            alt_start_coord <- provide_alignable_endpoints(chr_data,updated_end_points,unalignable_endpoints)
            out_text <- paste('End Coordinate Falls In Unalignable Region. End Coordinate of Downstream Region - ',alt_start_coord,sep='')
          }
          
          showFeedbackWarning(
            inputId='coord',
            text=out_text
          )
          req(unalignable_endpoints=='FALSE')     
        }
        
        # Pull out Start and End Region where supplied end points fall
        updated_end_points_start_id <- updated_end_points[[1,8]]
        updated_end_points_end_id <- updated_end_points[[nrow(updated_end_points),8]]
        
        # Pull out all AW Regions Within Start and End Region
        final_region_raw <- chr_data %>% dplyr::slice(which(chr_data[,8]==updated_end_points_start_id):which(chr_data[,8]==updated_end_points_end_id))
        
        # Extract Start and End Region + Adjust to Use Supplied End Points
        final_region <- final_region_raw %>% dplyr::slice(2:(n()-1 ))
        original_start_row <- return_original_start_endpoints(updated_end_points,query_lineage,as.numeric(start_input))
        original_end_row <- return_original_end_endpoints(updated_end_points,query_lineage,as.numeric(end_input))
        
        final_region <- rbind(final_region,original_start_row,original_end_row) %>% arrange(B73_StartPos)
        colnames(final_region) <- c('ID_Chr','ID_StartPos','ID_EndPos','ASM_Chr','ASM_StartPos','ASM_EndPos','Reformat_Type','AW_Block')
        
        # Format AW Regions for Visualization
        final_region$ADJ_ID_StartPos <- final_region[,2] - as.numeric(final_region[[1,2]])
        final_region$ADJ_ID_EndPos <- final_region[,3] - as.numeric(final_region[[1,2]])
        final_region$ADJ_ASM_StartPos <- final_region[,5] - as.numeric(final_region[[1,5]])
        final_region$ADJ_ASM_EndPos <- final_region[,6] -as.numeric(final_region[[1,5]])
        
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
      nv_tracks <- anchorwave_tsv() %>% subset(Reformat_Type=='alignable_region')
      
      if (nrow(nv_tracks)!=0) {
        lineage_inputs <- c(input$query,input$comp)
        
        ID_lineage <- sort(lineage_inputs)[1]
        ASM_lineage <- sort(lineage_inputs)[2]
        
        ID_links <- nv_tracks %>% select(contains('ID'))
        ID_nv_tracks <- cbind(seq_id=ID_lineage,ID_links) %>% select(1,5,6)
        colnames(ID_nv_tracks) <- c("seq_id",'start','end')
        
        ASM_links <- nv_tracks %>% select(contains('ASM'))
        ASM_nv_tracks <- cbind(seq_id=ASM_lineage,ASM_links) %>% select(1,5,6)
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
          sub_ID_gff_exons$Strand,sub_ID_gff_exons$Canonical_Tx
        )
        
        
        colnames(gene_track) <- c('seq_id','start','end','strand','name')
        
        gene_track <- gene_track %>% add_row(seq_id=ASM_lineage,
                                             start=sub_ASM_gff_exons$Start,
                                             end=sub_ASM_gff_exons$End,
                                             strand=sub_ASM_gff_exons$Strand,
                                             name=sub_ASM_gff_exons$Canonical_Tx)
      }
      
      colnames(gene_track) <- c("Lineage","Start",'Stop','Strand','Transcript')
      gene_track <- gene_track %>% distinct(.keep_all=TRUE)
      return(gene_track)
    })
    
    ## This section parses the EDTA TE annotations to find
    ## fully encapsulated TEs in the region and then stores their coordinates for visualization
    TE_track <- eventReactive(input$go,{
      #TE_anno_dir <- './TE_anno/'
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
      ID_TEAnno_file <- paste(paste(ID_lineage,chr_input,sep='_'),'_EDTA_filtered.bed',sep='')
      ID_TEAnno_full_path <- paste(ID_TEAnno_subdir,ID_TEAnno_file,sep='')
      
      ID_TEAnno <- save_object(ID_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ID_TEAnno_file)
      
      ID_TEAnno_subregion <- ID_TEAnno %>% subset(start >= ID_start_pos) %>%
        subset(end <= ID_end_pos)
      
      if (nrow(ID_TEAnno_subregion!=0)) {
        ID_present_TEs <- data.frame(
          seq_id=ID_lineage,
          start=ID_TEAnno_subregion$start - ID_start_pos,
          end=ID_TEAnno_subregion$end - ID_start_pos,
          type='transposable_element'
        ) 
      } else {
        ID_present_TEs <- data.frame(seq_id=as.character(),start=as.numeric(),end=as.numeric(),type=as.character())
      }
      
      
      ASM_TEAnno_subdir <- paste(TE_anno_dir,ASM_lineage,'/',sep='')
      ASM_TEAnno_file <- paste(paste(ASM_lineage,chr_input,sep='_'),'_EDTA_filtered.bed',sep='')
      ASM_TEAnno_full_path <- paste(ASM_TEAnno_subdir,ASM_TEAnno_file,sep='')
      
      ASM_TEAnno <- save_object(ASM_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_TEAnno_file)
      
      ASM_TEAnno_subregion <- ASM_TEAnno %>% subset(start >= ASM_start_pos) %>%
        subset(end <= ASM_end_pos)
      
      if (nrow(ASM_TEAnno_subregion!=0)) {
        ASM_present_TEs <- data.frame(
          seq_id=ASM_lineage,
          start=ASM_TEAnno_subregion$start - ASM_start_pos,
          end=ASM_TEAnno_subregion$end - ASM_start_pos,
          type='transposable_element'
        )  
      } else {
        ASM_present_TEs <- data.frame(seq_id=as.character(),start=as.numeric(),end=as.numeric(),type=as.character())
      }
      
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
      ID_TEAnno_file <- paste(paste(ID_lineage,chr_input,sep='_'),'_EDTA_filtered.bed',sep='')
      ID_TEAnno_full_path <- paste(ID_TEAnno_subdir,ID_TEAnno_file,sep='')
      
      ID_TEAnno <- save_object(ID_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ID_TEAnno_file)
      
      ID_TEAnno_subregion <- ID_TEAnno %>% subset(start >= ID_start_pos) %>%
        subset(end <= ID_end_pos)
      
      output_ID_TEAnno <- ID_TEAnno_subregion %>% 
        select(lineage,chr,start,end,strand,name,type,condense_superfamily,raw_family,method)
      
      colnames(output_ID_TEAnno) <- c('Lineage','Chr','Start','End','Strand','EDTA_Name','Type','Superfamily','Family','Method')
      
      ASM_TEAnno_subdir <- paste(TE_anno_dir,ASM_lineage,'/',sep='')
      ASM_TEAnno_file <- paste(paste(ASM_lineage,chr_input,sep='_'),'_EDTA_filtered.bed',sep='')
      ASM_TEAnno_full_path <- paste(ASM_TEAnno_subdir,ASM_TEAnno_file,sep='')
      
      ASM_TEAnno <- save_object(ASM_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_TEAnno_file)
      
      ASM_TEAnno_subregion <- ASM_TEAnno %>% subset(start >= ASM_start_pos) %>%
        subset(end <= ASM_end_pos)
      
      output_ASM_TEAnno <- ASM_TEAnno_subregion %>% 
        select(lineage,chr,start,end,strand,name,type,condense_superfamily,raw_family,method)
      
      colnames(output_ASM_TEAnno) <- c('Lineage','Chr','Start','End','Strand','EDTA_Name','Type','Superfamily','Family','Method')
      
      output_TE_track <- rbind(output_ID_TEAnno, output_ASM_TEAnno)
      output_TE_track <- output_TE_track %>% distinct(.keep_all=TRUE)
      return(output_TE_track)
      
    })
    
    #Check if Unalignable Region 
    check_unalignableregion_warning_message <- eventReactive(input$go,{
      AW_df <- anchorwave_tsv()[,c(1,2,3,4,5,6,7)] %>% distinct(.keep_all=TRUE)
      if ('unalignable' %in% AW_df$Reformat_Type) {
        warning_message <- 'WARNING: Unalignable Region in Selected Region. Refer to Summarised AnchorWave Output.'
      } else {
        warning_message <- ''
      }
      return(warning_message)
    })
 
    #Check if ID Start/End Falls Within TE
    check_IDTE_feature_warning_message <- eventReactive(input$go,{
      TE_anno_dir <- 's3://shiny-namsv-data-storage/TE_anno/'
      
      lineage_inputs <- c(input$query,input$comp)
      
      ID_lineage <- sort(lineage_inputs)[1]
      
      ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
      ID_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3])
      
      coord_inputs <- input$coord
      chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
      range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
      start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
      end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
      
      ID_TEAnno_subdir <- paste(TE_anno_dir,ID_lineage,'/',sep='')
      ID_TEAnno_file <- paste(paste(ID_lineage,chr_input,sep='_'),'_EDTA_filtered.bed',sep='')
      ID_TEAnno_full_path <- paste(ID_TEAnno_subdir,ID_TEAnno_file,sep='')
      
      ID_TEAnno <- save_object(ID_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ID_TEAnno_file)
      
      start_betweenTE <- ID_TEAnno %>% subset(between(ID_start_pos,start,end))
      
      if (nrow(start_betweenTE) >0) {
        betweenTE_name <- start_betweenTE %>% subset(start==min(start)) %>% pull(name)
        betweenTE_name <- betweenTE_name[1] # just in case two TEs with same start
        betweenTE_start <- start_betweenTE %>% subset(start==min(start)) %>% pull(start)
        betweenTE_start <-betweenTE_start[1]
        betweenTE_start <- betweenTE_start-1
        start_warning <- paste("Start Position in ",ID_lineage,' falls within a TE (',betweenTE_name,'): Consider Adjusting ',ID_lineage, ' Start Position to ',betweenTE_start,sep='')
      } else {
        start_warning <- ''
      }
      
      end_betweenTE <- ID_TEAnno %>% subset(between(ID_end_pos,start,end))
      if (nrow(end_betweenTE) >0) {
        betweenTE_name <- end_betweenTE %>% subset(end==max(end)) %>% pull(name)
        betweenTE_name <- betweenTE_name[1] # just in case two TEs with same start
        betweenTE_end <- end_betweenTE %>% subset(end==max(end)) %>% pull(end)
        betweenTE_end <-betweenTE_end[1]
        betweenTE_end <- betweenTE_end + 1
        end_warning <- paste("End Position in ",ID_lineage,' falls within a TE (',betweenTE_name,'): Consider Adjusting ',ID_lineage, ' End Position to ',betweenTE_end,sep='')
      } else {
        end_warning <- ''
      }      
      
      ID_TE_warning <- c(start_warning,end_warning)
      return(ID_TE_warning)
    })
 
    #Check if ID Start/End Falls Within a Gene
    check_IDGene_feature_warning_message <- eventReactive(input$go,{
      parent_dir <- 's3://shiny-namsv-data-storage/gff_gene_info/'
      
      lineage_inputs <- c(input$query,input$comp)
      
      ID_lineage <- sort(lineage_inputs)[1]
      
      ID_start_pos <- as.numeric(anchorwave_tsv()[1,2])
      ID_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),3])
      
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
      
      gene_coords <- ID_gff_exons  %>% group_by(Canonical_Tx) %>% arrange(Start) %>% filter(row_number()==1 | row_number()==n()) %>% select(Chr,Start,End,Canonical_Tx) %>% ungroup() %>% group_by(Chr,Canonical_Tx) %>% dplyr::summarise(Start=min(Start),End=max(End))
      
      start_betweenGene <- gene_coords %>% subset(between(ID_start_pos,Start,End))
      
      if (nrow(start_betweenGene) >0) {
        betweenGene_name <- start_betweenGene %>% subset(Start==min(Start)) %>% pull(Canonical_Tx)
        betweenGene_name <- betweenGene_name[1] # just in case two TEs with same start
        betweenGene_start <- start_betweenGene %>% subset(Start==min(Start)) %>% pull(Start)
        betweenGene_start <-betweenGene_start[1]
        betweenGene_start <- betweenGene_start -1
        start_warning <- paste("Start Position in ",ID_lineage,' falls within a Gene (',betweenGene_name,'): Consider Adjusting ',ID_lineage, ' Start Position to ',betweenGene_start,sep='')
      } else {
        start_warning <- ''
      }
      
      end_betweenGene <- gene_coords %>% subset(between(ID_end_pos,Start,End))
      if (nrow(end_betweenGene) >0) {
        betweenGene_name <- end_betweenGene %>% subset(end==max(end)) %>% pull(name)
        betweenGene_name <- betweenGene_name[1] # just in case two Genes with same start
        betweenGene_end <- end_betweenGene %>% subset(end==max(end)) %>% pull(end)
        betweenGene_end <-betweenGene_end[1]
        betweenGene_end <- betweenGene_end+1
        end_warning <- paste("End Position in ",ID_lineage,' falls within a Gene (',betweenGene_name,'): Consider Adjusting ',ID_lineage, ' End Position to ',betweenGene_end,sep='')
      } else {
        end_warning <- ''
      }
      
      ID_Gene_warning <- c(start_warning,end_warning)
      return(ID_Gene_warning)
    })    
    
    #Check if ASM Start/End Falls Within TE
    check_ASMTE_feature_warning_message <- eventReactive(input$go,{
      TE_anno_dir <- 's3://shiny-namsv-data-storage/TE_anno/'
      
      lineage_inputs <- c(input$query,input$comp)
      
      ASM_lineage <- sort(lineage_inputs)[2]
      
      ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
      ASM_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6])
      
      coord_inputs <- input$coord
      chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
      range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
      start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
      end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
      
      ASM_TEAnno_subdir <- paste(TE_anno_dir,ASM_lineage,'/',sep='')
      ASM_TEAnno_file <- paste(paste(ASM_lineage,chr_input,sep='_'),'_EDTA_filtered.bed',sep='')
      ASM_TEAnno_full_path <- paste(ASM_TEAnno_subdir,ASM_TEAnno_file,sep='')
      
      ASM_TEAnno <- save_object(ASM_TEAnno_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_TEAnno_file)
      
      start_betweenTE <- ASM_TEAnno %>% subset(between(ASM_start_pos,start,end))
      
      if (nrow(start_betweenTE) >0) {
        betweenTE_name <- start_betweenTE %>% subset(start==min(start)) %>% pull(name)
        betweenTE_name <- betweenTE_name[1] # just in case two TEs with same start
        betweenTE_start <- start_betweenTE %>% subset(start==min(start)) %>% pull(start)
        betweenTE_start <-betweenTE_start[1]
        betweenTE_start <- betweenTE_start - 1
        start_warning <- paste("Start Position in ",ASM_lineage,' falls within a TE (',betweenTE_name,'): Consider Adjusting ',ASM_lineage, ' Start Position to ',betweenTE_start,sep='')
      } else {
        start_warning <- ''
      }
      
      end_betweenTE <- ASM_TEAnno %>% subset(between(ASM_end_pos,start,end))
      if (nrow(end_betweenTE) >0) {
        betweenTE_name <- end_betweenTE %>% subset(end==max(end)) %>% pull(name)
        betweenTE_name <- betweenTE_name[1] # just in case two TEs with same start
        betweenTE_end <- end_betweenTE %>% subset(end==max(end)) %>% pull(end)
        betweenTE_end <-betweenTE_end[1]
        betweenTE_end <- betweenTE_end + 1
        end_warning <- paste("End Position in ",ASM_lineage,' falls within a TE (',betweenTE_name,'): Consider Adjusting ',ASM_lineage, ' End Position to ',betweenTE_end,sep='')
      } else {
        end_warning <- ''
      }      
      
      ASM_TE_warning <- c(start_warning,end_warning)
      return(ASM_TE_warning)
    })        

    #Check if ID Start/End Falls Within a Gene
    check_ASMGene_feature_warning_message <- eventReactive(input$go,{
      parent_dir <- 's3://shiny-namsv-data-storage/gff_gene_info/'
      
      lineage_inputs <- c(input$query,input$comp)
      
      ASM_lineage <- sort(lineage_inputs)[2]
      
      ASM_start_pos <- as.numeric(anchorwave_tsv()[1,5])
      ASM_end_pos <- as.numeric(anchorwave_tsv()[nrow(anchorwave_tsv()),6])
      
      coord_inputs <- input$coord
      chr_input <- unlist(str_split(tolower(coord_inputs),pattern =':'))[1]
      range_input <- gsub(',','',unlist(str_split(coord_inputs,pattern =':'))[2])
      start_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[1])
      end_input <- as.numeric(unlist(str_split(range_input,pattern='[..]'))[3])
      
      ASM_subdir <- paste(parent_dir,ASM_lineage,'/',sep='')
      ASM_file <- paste(paste(chr_input,ASM_lineage,sep='_'),'.gff',sep='')
      ASM_full_path <- paste(ASM_subdir,ASM_file,sep='')
      
      ASM_gff_exons <- save_object(ASM_full_path,base_url = 'host_url', region = 'host_region',key='access_ID',secret='access_key')  %>% data.table::fread()
      unlink(ASM_file)
      
      gene_coords <- ASM_gff_exons  %>% group_by(Canonical_Tx) %>% arrange(Start) %>% filter(row_number()==1 | row_number()==n()) %>% select(Chr,Start,End,Canonical_Tx) %>% ungroup() %>% group_by(Chr,Canonical_Tx) %>% dplyr::summarise(Start=min(Start),End=max(End))
      
      start_betweenGene <- gene_coords %>% subset(between(ASM_start_pos,Start,End))
      
      if (nrow(start_betweenGene) >0) {
        betweenGene_name <- start_betweenGene %>% subset(Start==min(Start)) %>% pull(Canonical_Tx)
        betweenGene_name <- betweenGene_name[1] # just in case two TEs with same start
        betweenGene_start <- start_betweenGene %>% subset(Start==min(Start)) %>% pull(Start)
        betweenGene_start <-betweenGene_start[1]
        betweenGene_start <- betweenGene_start -1
        start_warning <- paste("Start Position in ",ASM_lineage,' falls within a Gene (',betweenGene_name,'): Consider Adjusting ',ASM_lineage, ' Start Position to ',betweenGene_start,sep='')
      } else {
        start_warning <- ''
      }
      
      end_betweenGene <- gene_coords %>% subset(between(ASM_end_pos,Start,End))
      if (nrow(end_betweenGene) >0) {
        betweenGene_name <- end_betweenGene %>% subset(end==max(end)) %>% pull(name)
        betweenGene_name <- betweenGene_name[1] # just in case two Genes with same start
        betweenGene_end <- end_betweenGene %>% subset(end==max(end)) %>% pull(end)
        betweenGene_end <-betweenGene_end[1]
        betweenGene_end <- betweenGene_end + 1
        end_warning <- paste("End Position in ",ASM_lineage,' falls within a Gene (',betweenGene_name,'): Consider Adjusting ',ASM_lineage, ' End Position to ',betweenGene_end,sep='')
      } else {
        end_warning <- ''
      }
      
      ASM_Gene_warning <- c(start_warning,end_warning)
      return(ASM_Gene_warning)
    })    

    #Output Warning Messages
    output$feature_warning_message <- renderText({
      
      unalignable_warning_message <- check_unalignableregion_warning_message()
      ID_TEwarning_messages <- check_IDTE_feature_warning_message()
      ASM_TEwarning_messages <- check_ASMTE_feature_warning_message()
      ID_Genewarning_messages <- check_IDGene_feature_warning_message()
      ASM_Genewarning_messages <- check_ASMGene_feature_warning_message()
      warning_messages <- c(unalignable_warning_message,ID_TEwarning_messages,ASM_TEwarning_messages,ID_Genewarning_messages,ASM_Genewarning_messages)
      warning_messages <- warning_messages[!warning_messages%in% c('')]
      warning_messages <- paste(warning_messages,collapse='\n')
      print(warning_messages)
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
      anchorwave_tsv()[,c(1,2,3,4,5,6,7)] %>% distinct(.keep_all=TRUE)
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
