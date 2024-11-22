# MR2-cML-SuSiE

### Overview
This repository provides a demonstration on how to use the MR2-cML-SuSiE `R` source code. The default option of the code is dependent on the OpenGWAS database<sup>[1]</sup>, but users can also provide their own harmonized Mendelian randomization (MR) data (see TLDR towards the bottom of this README).

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
step1.res1 <- mr2.cml.susie.step1(exposure.ids, outcome.id1, sample.sizes)
step1.res2 <- mr2.cml.susie.step1(exposure.ids, outcome.id2, sample.sizes)
```

which, upon finishing looks like this:
```
head(step1.res1)
[1] 0.401999338 0.172905493 0.963813629 0.026442776 0.013637944 0.001762121

head(step1.res2)
[1] 1.613439e-01 5.540860e-01 4.398653e-01 5.708226e-07 1.430179e-01 6.191650e-02
```

Since it can take a while to run UVMR-cML<sup>[5]</sup> on 249 metabolites for both traits, the end results of this step are provided in this Github for convenience and can be loaded using:
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
In addition, we need to prepare the list of outcome IDs in a list format:
```
outcome.id.list <- vector("list", length = 2)
outcome.id.list[[1]] <- outcome.id1
outcome.id.list[[2]] <- outcome.id2
```
Now, we can use the `harmonize.mr2.data` function created for obtaining harmonized data appropriate for multi-response Mendelian randomization (MR2) analysis from OpenGWAS:
```
mr2dat <- harmonize.mr2.data(exposure.ids.subset, outcome.id.list, sample.sizes.subset)
```
which in this case, provides two harmonized dataset `mr2dat$mvdat.list[[1]]` (for AD) and mr2dat$mvdat.list[[2]] (for HTN), as well as the sample sizes for both exposures (`mr2dat$exposure.sample.sizes`) and outcomes (`mr2dat$outcome.sample.sizes`). Notice that the correct AD sample size is 487511 (OpenGWAS apparently stores the wrong sample size 85934 for this particular dataset). Thus, we need to manually modify the harmonized data:
```
mr2dat$outcome.sample.sizes[1] <- 487511
```
Again, the harmonized MR2 data (together with correct AD sample size information) has been provided for convenience:
```
mr2dat <- readRDS("mr2dat.RDS")
```
With this, we can run step 2 of MR2-cML-SuSiE:
```
step2.res1 <- mr2.cml.susie.step2(mr2dat, 1)
step2.res2 <- mr2.cml.susie.step2(mr2dat, 2)
```
which is also designed to be run on a per trait basis. Upon finish running, this provides `invalid.idx1` (from `step2.res1`) and `invalid.idx2` (from `step2.res2`), i.e., the instrumental variables (IVs) deemed invalid using UVMR-cML for the corresponding outcomes. We can obtain an overall set of invalid IVs by taking the union:
```
invalid.idx <- sort(unique(c(invalid.idx1, invalid.idx2)))
```
Using this set of invalid IVs, we can obtain better initial value estimates for the IMS step (step 4) via MVMR-cML-SuSiE in step 3. Notice that we also need `rho.mat` (the genetic correlation matrix between exposures and outcomes) in step 3 (more details can be found [here](https://github.com/lapsumchan/MVMR-cML-SuSiE), and this has been provided:
```
rho.mat <- matrix(0, 250, 250)
rho.mat[1:249,1:249] <- readRDS("metdrho.RDS")
rho.mat[250,250] <- 1

rho.mat <- rho.mat[c(subset.idx,250),c(subset.idx,250)]

invalid.idx <- sort(unique(c(step2.res1$invalid.idx, step2.res2$invalid.idx)))

step3.res1 <- mr2.cml.susie.step3(step2.res1$mvdat, invalid.idx, step2.res1$theta.vec, rho.mat)
step3.res2 <- mr2.cml.susie.step3(step2.res2$mvdat, invalid.idx, step2.res2$theta.vec, rho.mat)
```

Finally, the IMS step can be run using the initial value provided:
```
theta.vec1 <- unname(coef(step3.res1))[-1]
theta.vec2 <- unname(coef(step3.res2))[-1]
```

We need to provide `mvdat.list` (a length `Q` list of harmonized MR datasets), and `theta.vec.list` (a length `Q` list of initial value estimates for each of the corresponding outcomes):
```
mvdat.list <- vector("list", length = 2)
mvdat.list[[1]] <- step2.res1$mvdat
mvdat.list[[2]] <- step2.res2$mvdat

theta.vec.list <- vector("list", length = 2)
theta.vec.list[[1]] <- theta.vec1
theta.vec.list[[2]] <- theta.vec2

res <- mr2.cml.susie.step4(mvdat.list, invalid.idx, theta.vec.list)
```

With that, we can identify the significant exposures checking if the row sum of the PIP matrix (row corresponding to an exposure) is greater than 1/157:
```
idx <- which(rowSums(res) > 1/157)
```
and the significant exposure set is given by `res[idx,]` which looks like:
```
> head(res[idx,])
                          AD         HTN      AD_HTN
met-d-IDL_TG    0.0001776506 0.002272474 0.007515827
met-d-L_VLDL_C  0.0001741017 0.003004226 0.009922033
met-d-L_VLDL_FC 0.0001464217 0.004703768 0.013175401
met-d-L_VLDL_L  0.0001352467 0.005774971 0.014877285
met-d-L_VLDL_P  0.0001412142 0.005913658 0.016024885
met-d-L_VLDL_PL 0.0001411858 0.005783453 0.015657510
```
Note that by checking column sum of this submatrix, we also can learn which one is the most likely configuration:
```
> colSums(res[idx,])
         AD         HTN      AD_HTN 
0.005055941 0.250410900 0.535365874
```
This indicates that this set of 52 metabolite exposures are more likely to be causal to both AD and HTN, than HTN alone. Depending on the actual IMS results, we can re-run UVMR-cML<sup>[5]</sup> or MVMR-cML<sup>[6]</sup> for each of the plausible models and quantify the statistical evidence for the exposure jointly contributing to both diseases. We use one of the 52 metabolites `S_VLDL_P` as an example. Since the 52 metabolites only come from a single cluster, we run UVMR-cML:

```
exposure.dat.sub <- extract_instruments("met-d-S_VLDL_P")
outcome.dat.sub1 <- extract_outcome_data(snps = exposure.dat.sub$SNP, outcomes = outcome.id1)
outcome.dat.sub2 <- extract_outcome_data(snps = exposure.dat.sub$SNP, outcomes = outcome.id2)
dat1 <- harmonise_data(exposure.dat.sub, outcome.dat.sub1)
dat2 <- harmonise_data(exposure.dat.sub, outcome.dat.sub2)

# Perform UVMR-cML
cML.result1 <- mr_cML(dat1$beta.exposure,
                      dat1$beta.outcome,
                      dat1$se.exposure,
                      dat1$se.outcome,
                      n = n,
                      random_start = 100,
                      random_seed = 1)

cML.result2 <- mr_cML(dat2$beta.exposure,
                      dat2$beta.outcome,
                      dat2$se.exposure,
                      dat2$se.outcome,
                      n = n,
                      random_start = 100,
                      random_seed = 1)
```
and we can combine the BIC-based *p*-values for each outcome using the Cauchy combination test:
```
> cauchy.pvalue(c(cML.result1$BIC_p, cML.result2$BIC_p))
[1] 7.847029e-22
```

# TLDR
Users can provide their own version of harmonized data. For step 1, we require `Q` length `L` lists of summary statistics coefficients (beta) and standard errors (se) for both the exposures and outcomes, where `Q` is the number of outcomes analyzed. This is basically providing the univariable MR (UVMR) harmonized data for each exposure (and the outcomes summary statistics corresponding to the IVs used). Notice that the set of IVs (per exposure-outcome pair) should be independent (can be achieved by LD clumping), as this is a requirement for the cML framework: UVMR-cML<sup>[5]</sup>, MVMR-cML<sup>[6]</sup>, MVMR-cML-SuSiE<sup>[7]</sup> as well as our method (which builds upon on these former methods). In our case, `L = 249`. In addition, we also require `Q` vectors of length `L + 1` containing the sample sizes for each of the `L` exposures and the outcome GWAS (last element of each vector corresponds always correspond to an outcome). The `metdn.RDS` file contains sample sizes for the 249 UKB exposures, while 487511 and 484598 are the sample sizes for the AD and HTN GWAS, respectively. Below shows all 10 objects (5 per outcome) required for step 1 if the users were to provide their own data:

```
sample.sizes <- readRDS("metdn.RDS")
sample.sizes1 <- c(sample.sizes, 487511)
sample.sizes2 <- c(sample.sizes, 484598)

beta.exposure.ls1 <- readRDS("beta.exposure.ls1.RDS")
se.exposure.ls1 <- readRDS("se.exposure.ls1.RDS")
beta.outcome.ls1 <- readRDS("beta.outcome.ls1.RDS")
se.outcome.ls1 <- readRDS("se.outcome.ls1.RDS")

beta.exposure.ls2 <- readRDS("beta.exposure.ls2.RDS")
se.exposure.ls2 <- readRDS("se.exposure.ls2.RDS")
beta.outcome.ls2 <- readRDS("beta.outcome.ls2.RDS")
se.outcome.ls2 <- readRDS("se.outcome.ls2.RDS")
```

After running
```
step1.res1 <- mr2.cml.susie.step1(sample.sizes = sample.sizes1, beta.exposure.ls = beta.exposure.ls1, se.exposure.ls = se.exposure.ls1, beta.outcome.ls = beta.outcome.ls1, se.outcome.ls = se.outcome.ls1, use.openGWAS = FALSE)
step1.res2 <- mr2.cml.susie.step2(sample.sizes = sample.sizes2, beta.exposure.ls = beta.exposure.ls2, se.exposure.ls = se.exposure.ls2, beta.outcome.ls = beta.outcome.ls2, se.outcome.ls = se.outcome.ls2, use.openGWAS = FALSE)
```
with the `use.openGWAS` option specified as `FALSE`, it should same results as the OpenGWAS dependent version in the README.

Based on these step 1 results, it should suggest a subset of `L.star` exposures (in this case, 157) that warrants further investigation:
```
step1.res.list <- vector("list", length = 2)
step1.res.list[[1]] <- step1.res1
step1.res.list[[2]] <- step1.res2
subset.idx <- identify.exposure.subset.idx(step1.res.list)
sample.sizes.subset <- sample.sizes[subset.idx]
```

The users will then need to provide a harmonized MR2 dataset, while is a list of 3 objects: `mvdat.list`, `exposure.sample.sizes` and `outcome.sample.sizes`: `mvdat.list` is a length `Q` list of harmonized summary statistics, each `m.star` x `L.star` (i.e., same set of `m.star` IVs for all the outcomes), each containing a list of `exposure_beta`, `exposure_se`, `exposure_pval`, `outcome_beta` and `outcome_se`; `exposure.sample.sizes` is a length `L.star` vector coming from `sample.sizes.subset`; and `outcome.sample.sizes`, a length `Q` vector for the sample sizes of the outcomes (in the same order as `mvdat.list`), `c(487511, 484598)` in this case. This has been provided in the file `mr2dat.RDS` as an example and step 2 can be run using the `mr2.cml.susie.step2` function. Step 3 and 4 does not depend on OpenGWAS, so please refer back to the above README. 


### References

[1] Elsworth, Ben, et al. "The MRC IEU OpenGWAS data infrastructure." BioRxiv (2020): 2020-08.

[2] Borges, Maria Carolina, et al. "Role of circulating polyunsaturated fatty acids on cardiovascular diseases risk: analysis using Mendelian randomization and fatty acid genetic association data from over 114,000 UK Biobank participants." BMC medicine 20.1 (2022): 1-14.

[3] Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias." Nature genetics 54.4 (2022): 412-436.

[4] Dönertaş, Handan Melike, et al. "Common genetic associations between age-related diseases." Nature aging 1.4 (2021): 400-412.

[5] Xue, Haoran, Xiaotong Shen, and Wei Pan. "Constrained maximum likelihood-based Mendelian randomization robust to both correlated and uncorrelated pleiotropic effects." The American Journal of Human Genetics 108.7 (2021): 1251-1269.

[6] Lin, Zhaotong, Haoran Xue, and Wei Pan. "Robust multivariable Mendelian randomization based on constrained maximum likelihood." The American Journal of Human Genetics 110.4 (2023): 592-605.

[7] Chan, Lap Sum, Mykhaylo M. Malakhov, and Wei Pan. "A novel multivariable Mendelian randomization framework to disentangle highly correlated exposures with application to metabolomics." The American Journal of Human Genetics 111.9 (2024): 1834-1847.

