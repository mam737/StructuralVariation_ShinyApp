# StructuralVariation_ShinyApp


In order to easily visualize structural variation between the 26 lines of the maize [Nested Association Mapping](https://www.science.org/doi/10.1126/science.abg5289) (NAM) population, I constructed a R Shiny App which reads in parsed and filtered gvcf outputs from [AnchorWave](https://github.com/baoxingsong/AnchorWave), gff gene annotations, and [EDTA](https://github.com/oushujun/EDTA). It then leverages the R package [gggenomes](https://github.com/thackl/gggenomes) to visualize nonvariant and structural variants ( > 100 bp) between two maize lines for a user supplied region.

As of 05/12/22, this app supports comparisons between B73 and any of the other 25 NAM lines. It currently only supports two line comparisons, but we are working on extensions that allow users to compare differences across several lines. The reference and comparison inputs do allow the user to pick any of the 26 NAM lines, but the app will break if one of those is not B73. The app only visualizes one region at a time, inputed as chr#:start..end.

The app produces a single visualizing showing nonvariant (connected by aquamarine blocks) and variant regions. We condense true nonvariant regions  with SNPs and small InDels (<100 bp) into one broader nonvariant region. It then overlays exon/intron structure (colored in dark blue) for the canonical transcripts of any genes in the user-selected region. Finally, it also overlays tracks showing annotated transposable elements (colored in yellow) fully encapsualted in the region as well. Further information on all aspects of the visualizations is outputted underneath in a series of tables. 

This Shiny app is very much still in development. I anticipate that there may be comparisons made that break the app. If you identify any trends in regions that break the app or would like me to investigate a specific region in particular, please contact me at mmunasin@umn.edu. 

## Find R Shiny App AnchorWave Identified Structural Variation in NAM [Here](https://mmunasin.shinyapps.io/nam_sv/)
