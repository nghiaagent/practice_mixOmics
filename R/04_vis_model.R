here::i_am("R/04_vis_model.R")

# Import packages
library(mixOmics)
library(tidyverse)

# Load data
breast_tcga_splsda <- readRDS(
  here::here("output/breast_tcga_splsda.rds")
)

# Sample plots
## Check correlations between component 1 of each dataset.
breast_tcga_splsda %>%
  plotDiablo(ncomp = 1)

## Projection of each sample to the first two components in each omics level.
breast_tcga_splsda %>%
  plotIndiv(
    ind.names = FALSE,
    legend = TRUE,
    title = "TCGA, block sPLS-DA comps 1 - 2"
  )

## Agreement between datasets at the sample level
## (centroid distance vs actual omics position)
breast_tcga_splsda %>%
  plotArrow(
    ind.names = FALSE,
    legend = TRUE
  )

# Variable plots
## Contribution of each variable to components 1 and 2 (similar to PCA loadings)
breast_tcga_splsda %>%
  plotVar(
    var.names = FALSE,
    legend = TRUE,
    title = "TCGA, block sPLS-DA comps 1 - 2"
  )

## Contributions of each variable to each component
## separated by omics level (loadings)
breast_tcga_splsda %>%
  plotLoadings(
    comp = 1,
    contrib = "max",
    method = "median"
  )

## Correlation between variables in different omics in circos plot
breast_tcga_splsda %>%
  circosPlot(
    cutoff = 0.7,
    line = TRUE,
    color.cor = c("red", "blue")
  )

## Correlation between variables in different omics in network plot
breast_tcga_splsda %>%
  network(
    blocks = c(1, 2, 3),
    cutoff = 0.4
  )

## Heatmap of multiomics signature expression for each sample.
breast_tcga_splsda %>%
  cimDiablo(
    comp = 1,
    margin = c(8, 20),
    legend.position = "right"
  )
