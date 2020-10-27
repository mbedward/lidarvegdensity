# lidarvegdensity
Provides functions to estimate vegetation density profiles from terrestrial LiDAR based on
the methods described in:

Ashcroft MB, Gollan JR, Ramp D (2014) *Creating vegetation density profiles for a diverse range of 
  ecological habitats using terrestrial laser scanning.*
  Methods in Ecology and Evolution 5, 263–272. doi:10.1111/2041-210X.12157.

The original paper included R code (as an online appendix) for two functions: the first to
derive a 'null' or 'reference' model from a scan performed in a closed room; and the second 
to process a vegetation site scan with the null model used to correct for non-returned pulses
when calculating vegetation density within voxels.

This package provides functions for the same steps, but with sections of the code written in 
Rcpp to greatly reduce the time taken to process scans. It also provides an option to calculate 
bounds, as a [Jeffrey's interval](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Jeffreys_interval),
around the proportion of reflected rays in each voxel.

## Installation

It is best to use R version 4.0.0 or newer.

If you are a Windows user, make sure that you have Rtools installed on your system to build the C++ code included in the package. 
See https://cran.r-project.org/bin/windows/Rtools/ for more information.

To download and build the package, run the following code from the R or RStudio console...

```
# Uncomment and run the next line if you don't have the devtools package installed
# install.packages("devtools")

devtools::install_github("mbedward/lidarvegdensity")

```
