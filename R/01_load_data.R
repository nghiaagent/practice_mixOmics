here::i_am("R/01_load_data.R")

# Import packages
library(mixOmics)
library(tidyverse)

# This projects demos block sPLS-DA for creating a classifier
# to discriminate Basal, Her2, and LumA breast cancer subtypes
# using miRNA, mRNA, proteomics data.
# Load mixOmics example dataset
data("breast.TCGA")

# Place training data in new list
breast_tcga_train <- list(
  mrna = breast.TCGA$data.train$mrna,
  mirna = breast.TCGA$data.train$mirna,
  protein = breast.TCGA$data.train$protein
)

breast_tcga_train_outcome <- breast.TCGA$data.train$subtype

# Save data
saveRDS(
  breast_tcga_train,
  here::here("output/breast_tcga_train.rds")
)

saveRDS(
  breast_tcga_train_outcome,
  here::here("output/breast_tcga_train_outcome.rds")
)
