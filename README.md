# Running MesKit on acral melanoma samples from Brazil

## Overview

This repository contains the code that was used to compare the patient and PDX-derived samples in terms of their variants and copy number alteration (CNA) events.

## Required software and environment

The following software is required to run the analysis:

- R (version 4.2.2)
- R packages:
    - `MesKit`
    - `here`
    - `ggpubr`
    - `dplyr`
    - `tidyr`
    - `stringr`
    - `logger`
    - `readr`
    - `fs`

The required packages can be installed using:
```R
install.packages(c("MesKit", "here", "ggpubr", "dplyr", "tidyr", "stringr", "logger", "readr", "fs"))
```

If you are interested in reproducing the R environment for the code used in R 4.2.2, an `renv.lock` file is available in this repository. To rebuild the environment, run the following inside R:
```R
# Install renv if needed
if (!require("renv")) install.packages("renv")

# Restore the project environment
renv::restore()
```

## Contact

If you have any questions or comments about this repository, please contact:

- Jacqueline Boccacino (<jb62@sanger.ac.uk>)
- Martin del Castillo-Velasco Herrera (<mdc1@sanger.ac.uk>)