
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MIWilson

<!-- badges: start -->

[![R-CMD-check](https://github.com/hungf8342/MIWilson/workflows/R-CMD-check/badge.svg)](https://github.com/hungf8342/MIWilson/actions)
<!-- badges: end -->

The goal of MIWilson is to implement the Wilson confidence interval for
binomial proportions given imputed missing data.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hungf8342/MIWilson")
```

## Basic Usage

The functions which calculate confidence intervals (mi\_wilson and
mi\_wald) take in a mids object (produced by the mice package), the
binary response variable name (the response must be 0-1 valued), and an
optional confidence level (default is 0.95).

As an example, let’s work with the hypertension (hyp) variable in the
nhanes toy dataset. The hyp variable refers to whether an individual has
hypertension. It originally has no and yes encoded as 1 and 2, so we
need to first re-value the hyp variable. Calculating Wilson and Wald
intervals for hyp is straightfoward afterwards.

``` r
library(MIWilson)

## setting up correct response variable values and mice object
nhanes = mice::nhanes %>%
  dplyr::mutate(hyp = hyp-1)
imp = mice::mice(nhanes)
#> 
#>  iter imp variable
#>   1   1  bmi  hyp  chl
#>   1   2  bmi  hyp  chl
#>   1   3  bmi  hyp  chl
#>   1   4  bmi  hyp  chl
#>   1   5  bmi  hyp  chl
#>   2   1  bmi  hyp  chl
#>   2   2  bmi  hyp  chl
#>   2   3  bmi  hyp  chl
#>   2   4  bmi  hyp  chl
#>   2   5  bmi  hyp  chl
#>   3   1  bmi  hyp  chl
#>   3   2  bmi  hyp  chl
#>   3   3  bmi  hyp  chl
#>   3   4  bmi  hyp  chl
#>   3   5  bmi  hyp  chl
#>   4   1  bmi  hyp  chl
#>   4   2  bmi  hyp  chl
#>   4   3  bmi  hyp  chl
#>   4   4  bmi  hyp  chl
#>   4   5  bmi  hyp  chl
#>   5   1  bmi  hyp  chl
#>   5   2  bmi  hyp  chl
#>   5   3  bmi  hyp  chl
#>   5   4  bmi  hyp  chl
#>   5   5  bmi  hyp  chl

## MI-Wilson and MI-Wald 99% CIs of the proportion of patients with hypertension 
mi_wilson(imp,"hyp", 0.99)
#> [1] "Qbar:  0.224"
#> [1] "RM:  0.362453531598513"
#> [1] 0.07594168 0.50344910
mi_wald(imp, "hyp", 0.99)
#> [1] -0.007904497  0.455904497
```

Helper functions (other than Qhats) do not take in mids objects directly
since they can be calculated using Qhats (the mean proportions of each
imputed dataset).

``` r
## Qhats calculates the mean proportion of hyp for each imputed dataset
qhats = Qhats(imp,"hyp")

## Keeping track of number of imputed datasets (m) and number of observations in a dataset
m = imp$m
nrow = imp$data %>% nrow()

## Qbar (mean of Qhats) and Ubar (average response variance over imputed datasets)
Qbar(qhats)
#> [1] 0.224
Ubar(qhats, m, nrow)
#> [1] 0.0068864
```
