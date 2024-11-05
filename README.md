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

We will illustrate our function via the UK Biobank (UKB) metabolite example used in our manuscript, and two outcomes of interest are Alzheimer's disease (AD) and hypertension (HTN). The summary statistics of the 249 UKB metabolites by Borges et al.<sup>[2]</sup> are available from the OpenGWAS database<sup>[1]</sup> with `met-d` prefix:

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
while the HTN GWAS summary statistics comes from Dönertaş et al.<sup>[4]</sup>, which is also available in OpenGWAS:
```
outcome.id2 <- "ebi-a-GCST90038604"
```

Notice that we need the minimum sample size amongst GWASs (exposure + outcomes) for cML. Thus, we need to prepare a vector of sample sizes corresponding to each exposures in `sample.sizes`. This has been prepared in the file `metdn.RDS` so we just need to load it:
```
sample.sizes <- readRDS("metdn.RDS")
```

With the above four `R` objects (`exposure.ids`, `outcome.id1`, `outcome.id2` and `sample.sizes`), we are ready to run step 1 of MR2-cML-SuSiE. The function is designed such that it can be run by each trait, and each run would subsequently provide a vector of *p*-values for all 249 metabolite exposures for a given trait:

```
step1.res1 <- mv.cml.susie.step1(exposure.ids, outcome.id1, sample.sizes)
step1.res2 <- mv.cml.susie.step1(exposure.ids, outcome.id2, sample.sizes)
```

which, upon finishing looks like this:
```
head(step1.res1)
[1] 0.401999338 0.172905493 0.963813629 0.026442776 0.013637944 0.001762121

head(step1.res2)
[1] 1.613439e-01 5.540860e-01 4.398653e-01 5.708226e-07 1.430179e-01 6.191650e-02
```

Since it can take a while to run UVMR-cML on 249 metabolites for both traits, the end results of this step are provided in this Github for convenience and can be loaded using:
```
step1.res1 <- readRDS("step1res1.RDS")
step1.res2 <- readRDS("step1res2.RDS")
```

Next, we create a list of p-values according to the number of traits analyzed:
```
step1.res.list <- vector("list", length = 2)
step1.res.list[[1]] <- step1.res1
step1.res.list[[2]] <- step1.res2
```
and use the function `identify.exposure.subset.idx` to identify the indices corresponding to the exposures that are worth further investigation:
```
subset.idx <- identify.exposure.subset.idx(step1.res.list)
```
which should suggest
```
length(subset.idx)
[1] 157
```
This means that 157 exposures will be jointly analyzed in the subsequent steps. Notice that the default used for screening is Bonferroni correction (0.05 divided by the product of number of traits and the number of exposures). Change the `cutoff` argument if you want to use something different. In this case, `cutoff` is 0.05 / (2 x 249) ~ 1e-4.

With this, we can extract the exposures which we wish to further investigate:
```
exposure.ids.subset <- exposure.ids[subset.idx]
```
and correspondingly, we can also obtain the sample sizes for the 157 exposures using:
```
sample.sizes.subset <- sample.sizes[subset.idx]
```

### References

[1] Elsworth, Ben, et al. "The MRC IEU OpenGWAS data infrastructure." BioRxiv (2020): 2020-08.

[2] Borges, Maria Carolina, et al. "Role of circulating polyunsaturated fatty acids on cardiovascular diseases risk: analysis using Mendelian randomization and fatty acid genetic association data from over 114,000 UK Biobank participants." BMC medicine 20.1 (2022): 1-14.

[3] Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias." Nature genetics 54.4 (2022): 412-436.

[4] Dönertaş, Handan Melike, et al. "Common genetic associations between age-related diseases." Nature aging 1.4 (2021): 400-412.
