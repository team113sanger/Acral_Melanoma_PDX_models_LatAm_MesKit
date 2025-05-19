# Running MesKit on acral melanoma samples from Brazil

## Required software and environment

The following software is required to run the analysis:

- R (version 4.2.2)
- R packages:
    - `MesKit`
    - `here`

The required packages can be installed using:
```R
install.packages(c("MesKit", "here"))
```

If you are interested in reproducing the R environment for the code used in R 4.2.2, an `renv.lock` file is available in this repository. To rebuild the environment, run the following inside R:
```R
# Install renv if needed
if (!require("renv")) install.packages("renv")

# Restore the project environment
renv::restore()
```
