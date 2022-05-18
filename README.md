# StructuralVariation_ShinyApp

## Last Updated - 05/18/22

In order to easily visualize structural variation between the 26 lines of the maize [Nested Association Mapping](https://www.science.org/doi/10.1126/science.abg5289) (NAM) population, I constructed a R Shiny App which reads in parsed and filtered gvcf outputs from [AnchorWave](https://github.com/baoxingsong/AnchorWave), gff gene annotations, and [EDTA](https://github.com/oushujun/EDTA). It then leverages the R package [gggenomes](https://github.com/thackl/gggenomes) to visualize nonvariant and structural variants ( > 75 bp) between two maize lines for a user supplied region.

As of 05/18/22, this app supports comparisons between B73 and any of the other 25 NAM lines. It currently only supports two line comparisons, but we are working on extensions that allow users to compare differences across several lines. The reference and comparison inputs do allow the user to pick any of the 26 NAM lines, but the app will break if one of those is not B73. Also, unfortunately, Oh43 is not available at this time, but it will be soon. I'll update this accordingly. The app only visualizes one region at a time, inputed as either chr#:1234..5678 or Chr#:1,234..5,678 allowing for easy copying of coordinates from either JBrowse or GBrowse in maizeGDB. 

The app produces a single visualizing showing nonvariant (connected by aquamarine blocks) and variant regions. We condense true nonvariant regions  with SNPs and small InDels (<75 bp) into one broader nonvariant region. It then overlays exon/intron structure (colored in dark blue) for the canonical transcripts of any genes in the user-selected region. Finally, it also overlays tracks showing annotated transposable elements (colored in yellow) fully encapsualted in the region as well. Further information on all aspects of the visualizations is outputted underneath in a series of tables. 

This Shiny app is very much still in development, and I anticipate that there may be comparisons made that break the app. If you identify any trends in regions that break the app or would like me to investigate a specific region in particular, please contact me at mmunasin@umn.edu. 

## Find R Shiny App AnchorWave Identified Structural Variation in NAM [Here](https://mmunasin.shinyapps.io/nam_sv/)

## Updates

05/18/22 - Accepts two forms of coordinate inputs to match those generated by JBrowse and GBrowse from maizeGDB

05/18/22 - Adjusted variant classification to call variants larger than 75 bp (formerly 100 bp) as structural variants

05/18/22 - Fixed bug that caused the app to break when visualizing regions with no genes

05/18/22 - Fixed bug that caused Oh7B to fail

05/18/22 - Updated app aesthetic using shinydashboard