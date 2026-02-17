source("renv/activate.R")

# Set up parallel computing specs
# Install BiocParallel via pak
if (!requireNamespace("BiocParallel")) {
  if (!requireNamespace("pak")) {
    install.packages("pak")
  }
  pak::pak("BiocParallel")
}

current_bpparam <- BiocParallel::SnowParam(workers = 4)