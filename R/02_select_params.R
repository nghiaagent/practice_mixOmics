here::i_am("R/02_select_params.R")

# Import packages
library(BiocParallel)
library(mixOmics)
library(tidyverse)

# Load data
breast_tcga_train <- readRDS(
  here::here("output/breast_tcga_train.rds")
)

breast_tcga_train_outcome <- readRDS(
  here::here("output/breast_tcga_train_outcome.rds")
)

# Select weights for sPLS-DA
# The current dataset contains paired miRNA, mRNA, and protein data
# The weight could be based on the pairwise correlation between them.
## Calculate pairwise correlations using PLS
correlations <- list(
  # Calculate PLS between pairs of omics
  mrna_protein = pls(
    breast_tcga_train$mrna,
    breast_tcga_train$protein,
    ncomp = 1
  ),
  mrna_mirna = pls(
    breast_tcga_train$mrna,
    breast_tcga_train$mirna,
    ncomp = 1
  ),
  mirna_protein = pls(
    breast_tcga_train$mirna,
    breast_tcga_train$protein,
    ncomp = 1
  )
) %>%
  # Calculate correlation within each PLS
  map(\(pls) cor(pls$variates$X, pls$variates$Y))

## As the correlation is 0.8 - 0.9, we can use a weight in this range.
## Define a design matrix with 0.85 as the weight between omics.
breast_tcga_design <- matrix(
  0.85,
  nrow = length(breast_tcga_train),
  ncol = length(breast_tcga_train),
  dimnames = list(names(breast_tcga_train), names(breast_tcga_train))
) %>%
  `diag<-`(0)

# Select number of components for sPLS-DA
## First, fit a test block PLS-DA (non-sparse)
## to choose the number of components
breast_tcga_plsda_test <- block.plsda(
  X = breast_tcga_train,
  Y = breast_tcga_train_outcome,
  design = breast_tcga_design,
  ncomp = 5
)

## Test and plot model statistics
## Use 10-fold cross-validation with 10 repeats
breast_tcga_plsda_test_perf <- breast_tcga_plsda_test %>%
  perf(validation = "Mfold", folds = 10, nrepeat = 10)

### Plot statistics
plot(breast_tcga_plsda_test_perf)

### Show the recommended number of components
breast_tcga_plsda_test_perf$choice.ncomp$WeightedVote %>%
  print()

### Select the highest recommended number of components (4)
breast_tcga_ncomp <- breast_tcga_plsda_test_perf$choice.ncomp$WeightedVote %>%
  max()

# Select variables for sPLS-DA
## Define tuning parameters; start from 5 variables then gradually increase
breast_tcga_tune_nvars <- list(
  mrna = c(seq(5, 9, 2), seq(10, 25, 5)),
  mirna = c(seq(5, 9, 2), seq(10, 25, 5)),
  protein = c(seq(5, 9, 2), seq(10, 25, 5))
)

breast_tcga_tune <- tune.block.splsda(
  X = breast_tcga_train,
  Y = breast_tcga_train_outcome,
  ncomp = 2,
  test.keepX = breast_tcga_tune_nvars,
  design = breast_tcga_design,
  validation = "Mfold",
  folds = 5,
  nrepeat = 2,
  dist = "centroids.dist",
  BPPARAM =
)

## Get selected variables
breast_tcga_keep_vars <- breast_tcga_tune$choice.keepX

# Save data
## Save the design matrix
saveRDS(
  breast_tcga_design,
  here::here("output/breast_tcga_design.rds")
)

## Save the number of components
saveRDS(
  breast_tcga_ncomp,
  here::here("output/breast_tcga_ncomp.rds")
)

## Save the selected variables
saveRDS(
  breast_tcga_keep_vars,
  here::here("output/breast_tcga_keep_vars.rds")
)
