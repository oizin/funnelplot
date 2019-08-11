## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 5,
  fig.width = 8
)

## ----setup---------------------------------------------------------------
library(funnelplot)
library(ggplot2)

## ------------------------------------------------------------------------
data("example_data", package = "funnelplot")
head(example_data)

## ------------------------------------------------------------------------
# outcome ~ covariates | cluster ID
f1 <- funnel(test ~ 1 | hosp_id, data = example_data)
plot(f1)

## ------------------------------------------------------------------------
pp <- plot(f1)
class(pp)
pp +
  theme_bw()

## ------------------------------------------------------------------------
outliers(f1)

## ------------------------------------------------------------------------
f2 <- funnel(test ~ var | hosp_id, data = example_data)
plot(f2)

## ------------------------------------------------------------------------
f3 <- funnel(test ~ var | hosp_id, data = example_data)
plot(f3)

## ------------------------------------------------------------------------
f4 <- funnel(test ~ var | hosp_id, data = example_data)
plot(f4)

## ------------------------------------------------------------------------
plot(f4,label="outliers")

## ------------------------------------------------------------------------
plot(f4,identify=13,label=13)

