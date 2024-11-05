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
