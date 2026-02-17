# TCGA breast cancer multi-omics analysis with mixOmics


# Introduction

This repository contains the source code, data processing scripts, and
analysis results for my mixOmics demo project.

See a compiled version of the report here:
<https://nghiaagent.github.io/practice_mixOmics/>

# Prerequisites

Make sure the following are installed on your system: - Quarto - R 4.5.1

# Reproduciblity

## Report

To reproduce the report, clone this repo then run the following in a
terminal:

    #| eval: False
    quarto render report.qmd

## Full project

In the interest of speed, some parts of the project, namely variable
selection, are excluded from the report. To run the full project, run
the following in an R terminal:

``` r
renv::restore()
source(here::here("R/99_first_run.R"))
```
