here::i_am("R/03_build_model.R")

# Import packages
library(mixOmics)
library(tidyverse)

# Load data
breast_tcga_train <- readRDS(
  here::here("output/breast_tcga_train.rds")
)

breast_tcga_train_outcome <- readRDS(
  here::here("output/breast_tcga_train_outcome.rds")
)

breast_tcga_design <- readRDS(
  here::here("output/breast_tcga_design.rds")
)
breast_tcga_ncomp <- readRDS(
  here::here("output/breast_tcga_ncomp.rds")
)
breast_tcga_keep_vars <- readRDS(
  here::here("output/breast_tcga_keep_vars.rds")
)

# Build final block sparse PLS-DA model
breast_tcga_splsda <- block.splsda(
  X = breast_tcga_train,
  Y = breast_tcga_train_outcome,
  design = breast_tcga_design,
  ncomp = breast_tcga_ncomp,
  keepX = breast_tcga_keep_vars
)

# Save data
saveRDS(
  breast_tcga_splsda,
  here::here("output/breast_tcga_splsda.rds")
)
