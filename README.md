
<!-- README.md is generated from README.Rmd. Please edit that file -->

# funnelplot

<!-- badges: start --> [![Travis build
status](https://travis-ci.org/oizin/funnelplot.svg?branch=master)](https://travis-ci.org/oizin/funnelplot)
[![codecov](https://codecov.io/github/oizin/funnelplot/branch/master/graphs/badge.svg)](https://codecov.io/github/oizin/funnelplot)
<!-- badges: end -->

The goal of funnelplot is to simplify risk adjusted comparisons of
institutional performance for the purpose of outlier discovery.

## Installation

You can install the development version of funnelplot from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("oizin/funnelplot")
```

## Example

This is a basic example which shows you how to produce a funnel plot:

``` r
library(funnelplot)
data("example_data")

# outcome ~ covariates | cluster ID
f1 <- funnel(test ~ 1 | hosp_id, data = example_data)
plot(f1)
```

<img src="man/figures/README-example-1.png" width="100%" />

See the vignettes:

``` r
# vignettes
vignette("funnelplot-intro")
```

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. Journal of
the royal statistical society. Series B (Methodological), 289-300.

Jones, H. E., Ohlssen, D. I., & Spiegelhalter, D. J. (2008). Use of the
false discovery rate when comparing multiple health care providers.
Journal of clinical epidemiology, 61(3), 232-240.

Spiegelhalter, D. J. (2005). Funnel plots for comparing institutional
performance. Statistics in medicine, 24(8), 1185-1202.
