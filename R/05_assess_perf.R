here::i_am("R/05_assess_perf.R")

# Import packages
library(mixOmics)
library(patchwork)
library(tidyverse)

# Load data
breast_tcga_splsda <- readRDS(
  here::here("output/breast_tcga_splsda.rds")
)

# Assess performance in training dataset
# Use 10-fold cross-validation with 10 repeats against the training data itself
breast_tcga_splsda_perf <- breast_tcga_splsda %>%
  perf(
    validation = "Mfold",
    folds = 10,
    nrepeat = 10,
    dist = "centroids.dist"
  )

## Show model performance
### Per-omics, per-cell type error rate
breast_tcga_splsda_perf$error.rate.per.class %>%
  print()

### Classification performance via majority vote
breast_tcga_splsda_perf$MajorityVote.error.rate %>%
  print()

### Classification performance via weighted vote
breast_tcga_splsda_perf$WeightedVote.error.rate %>%
  print()

## AUC plot per block
breast_tcga_splsda_auc_plots <- map(
  names(breast_tcga_splsda$X),
  \(omics_level) {
    ### Build list of AUC calculations
    object <- auroc(
      breast_tcga_splsda,
      roc.block = omics_level,
      roc.comp = 2,
      print = FALSE
    )

    ### Extract and return just the plot
    plot <- object[[str_c("graph.", omics_level)]][["comp2"]]

    return(plot)
  }
) %>%
  wrap_plots(nrow = 3)

plot(breast_tcga_splsda_auc_plots)

# Assess performance in test dataset
## Prepare test set, missing protein data
data(breast.TCGA)
breast_tcga_test <- list(
  mrna = breast.TCGA$data.test$mrna,
  mirna = breast.TCGA$data.test$mirna
)

## Perform prediction
breast_tcga_test_predict <- predict(
  breast_tcga_splsda,
  newdata = breast_tcga_test
)

## Check prediction performance using confusion matrix.
## Use 2 components in model.
breast_tcga_test_confmatrix <- get.confusion_matrix(
  truth = breast.TCGA$data.test$subtype,
  predicted = breast_tcga_test_predict$WeightedVote$centroids.dist[, 2]
)

## Show balanced error rate for the above
breast_tcga_test_confmatrix %>% get.BER()
