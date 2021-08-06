
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MIWilson

<!-- badges: start -->

[![R-CMD-check](https://github.com/hungf8342/MIWilson/workflows/R-CMD-check/badge.svg)](https://github.com/hungf8342/MIWilson/actions)
<!-- badges: end -->

The goal of MIWilson is to implement the Wilson confidence interval for
binomial proportions given multiple imputations of missing data.

TODO: if phats are 0, then output warning; is it the same as if imputed
values are the same? TODO: put together vignette (pdf thingy) TODO:
publish to CRAN

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hungf8342/MIWilson")
```

## Basic Usage: MIDs Argument

The functions which calculate confidence intervals (`mi_wilson` and
`mi_wald`) take in a mids object (produced by the mice package), the
binary response variable name (the response must be 0-1 valued), a
summaries print option (default is TRUE), and a confidence level
(default is 0.95).

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
#> [1] "Qbar:  0.232"
#> [1] "Rm:  0.602189781021898"
#> [1] "dof:  28.3152984389348"
#> [1] 0.07134426 0.54292320
mi_wald(imp, "hyp", 0.99)
#> [1] "Qbar:  0.232"
#> [1] "Tm:  0.0112384"
#> [1] "dof:  28.3152984389348"
#> [1] -0.02936932  0.49336932
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
#> [1] 0.232
Ubar(qhats, m, nrow)
#> [1] 0.0070144
```

## Basic Usage: P-hats Argument

If the user provides their imputed datasets (rather than relying on the
`mice` package), MI-Wilson allows them to input the corresponding vector
of observed binomial proportions instead of a mids object.
`mi_wilson_phat` and `mi_wald_phat` take in a vector of observed
binomial proportions (one proportion for each of \(m\) imputed
datasets), the number of total observations (must be constant across
datasets), a summaries print option (default is TRUE), and a confidence
level (default is 0.95).

``` r
mi_wilson(phats = c(0.2,0.23,0.22), n = 10)
#> [1] "Qbar:  0.216666666666667"
#> [1] "Rm:  0.0183474215320097"
#> [1] "dof:  6161.29288265307"
#> [1] 0.07684304 0.47892199
mi_wilson_phat(c(0,0,0),n=10)
#> Warning in mi_wilson_phat(c(0, 0, 0), n = 10): Imputed binomial proportions are
#> identical; degrees of freedom set to infinity.
#> [1] "Qbar:  0"
#> [1] "Rm:  0"
#> [1] "dof:  Inf"
#> [1] 1.387779e-17 2.129420e-01
```
