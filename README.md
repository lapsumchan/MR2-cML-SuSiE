# MR2-cML-SuSiE

### Overview
This repository provides a demonstration on how to use the MR2-cML-SuSiE `R` source code. Currently, the code is dependent on the OpenGWAS database<sup>[1]</sup>.

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using the MR2-cML-SuSiE source code, users should have `R` version 4.3.0 or higher, and several packages installed.

### Installation  

First, we need to install a few dependencies [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/), [`MRcML`](https://github.com/xue-hr/MRcML) and [`MVMRcML`](https://github.com/ZhaotongL/MVMR-cML), [`susieR`](https://github.com/stephenslab/susieR) and [`MESuSiE`](https://github.com/borangao/MESuSiE):  

    install.packages("remotes")
    remotes::install_github("MRCIEU/TwoSampleMR")

    install.packages("susieR")

    install.packages("devtools")
    devtools::install_github("xue-hr/MRcML")
    devtools::install_github("ZhaotongL/MVMR-cML")
    devtools::install_github("borangao/MESuSiE")

# Demo

We first load all the source code dependencies:

```
library(TwoSampleMR)
library(susieR)
library(MRcML)
library(MVMRcML)
library(MESuSiE)
```

and the source code containing all the main functions:

```
source("main.R")
```

We will illustrate our function via the UK Biobank (UKB) metabolite example used in our manuscript, and two outcomes of interest are Alzheimer's disease (AD) and coronary heart disease (CAD). The summary statistics of the 249 UKB metabolites by Borges et al.<sup>[2]</sup> are available from the OpenGWAS database<sup>[1]</sup> with `met-d` prefix:

```
ao <- available_outcomes()

# Use grep to find ids that start with "met-d"
metd.idx <- grep("^met-d", ao$id)
exposure.ids <- ao$id[metd.idx]

# Make the ids in alphabetical order (to match the ordering of summary statistics given by TwoSampleMR package, i.e., mvdat in step 2 below)
exposure.ids <- sort(exposure.ids)
```

The AD GWAS summary statistics comes from the largest AD cohort by Bellenguez et al.<sup>[3]</sup>, which is available in OpenGWAS with the corresponding ID:
```
outcome.id1 <- "ebi-a-GCST90027158"
```
while the CAD GWAS summary statistics comes from van der Harst et al.<sup>[4]</sup>, which is also available in OpenGWAS:
```
outcome.id2 <- "ebi-a-GCST005195"
```

Notice that we need the minimum sample size amongst GWASs (exposure + outcomes) for cML. Thus, we need to prepare a vector of sample sizes corresponding to each exposures in `sample.sizes`. This has been prepared in the file `metdn.RDS` so we just need to load it:
```
sample.sizes <- readRDS("metdn.RDS")
```

With the above four `R` objects (`exposure.ids`, `outcome.id1`, `outcome.id2` and `sample.sizes`), we are ready to run step 1 of MR2-cML-SuSiE. The function is designed such that it can be run by trait, and each run would subsequently provides a vector of *p*-values for all 249 metabolite exposures for a given trait:

```
step1.res1 <- mv.cml.susie.step1(exposure.ids, outcome.id1, sample.sizes)
step1.res2 <- mv.cml.susie.step1(exposure.ids, outcome.id2, sample.sizes)
```

which, upon finishing looks like this:
```
head(step1.res1)

head(step1.res2)
```
