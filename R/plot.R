
#' helper function for plotting the funnel plots
#'
#' Determines the position of the control limits at N values.
#'
#' @param funnelRes funnel plot object
#' @param lengthOut 500
#'
plotLimits <- function(funnelRes,lengthOut) {

  # limits
  mn <- min(funnelRes$results$n)
  mx <- max(funnelRes$results$n)
  rng <- floor(seq(from = mn,to = mx,length.out=lengthOut))
  limits <- funnelRes$control$limits

  # limits
  if (funnelRes$control$method == "point") {
    # are we controlling for over-dispersion?
    if(funnelRes$control$crtlOverDisp == TRUE){
      inflationFactor <- funnelRes$dispersion["inflationFactor"]
    } else {
      inflationFactor <- 1
    }

    rngLimits <- pointLimits(target = funnelRes$target,
                              n = rng,
                              N = nrow(funnelRes$results),
                              inflationFactor = inflationFactor,
                              control = funnelRes$control,
                              pval = funnelRes$results$pval)

    lower <- rngLimits$lower
    upper <- rngLimits$upper

  } else if (funnelRes$control$method == "distribution") {
    effectVar <- funnelRes$dispersion["effectVar"]
    rngLimits <- randomEffectLimits(target=funnelRes$target,
                                    n=rng,
                                    effectVar=effectVar,
                                    control=funnelRes$control)
    lower <- rngLimits$lower
    upper <- rngLimits$upper
  }

  tmp <- as.character(100 * (1-(limits)))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper_",tmp)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower_",tmp)

  out <- data.frame(n = rng, upper, lower)
  out
}


#' Graph funnel plot
#'
#' Produces a funnel plot using ggplot2.
#'
#' @param x funnel plot object
#' @param identify show on plot. default = all
#' @param label label on plot. default = none
#' @param lengthOut resolution of control limits. Number of points sampled to construct limit lines. default = 500.
#' @param ... Other arguments to plot.funnelRes.
#'
#' @export
plot.funnelRes <- function(x,identify="all",label="none",lengthOut=500,...) {
  ## checks on input args
  assertthat::assert_that(all(class(x) %in% c("list","funnelRes")))
  assertthat::assert_that(identify %in% c("all","outliers",x$results$id))
  assertthat::assert_that(label %in% c("outliers","none",x$results$id))
  assertthat::assert_that(lengthOut > 0)

  # calculate control limits for plotting
  rngLimits  <- plotLimits(x,lengthOut=lengthOut)

  # long form
  rngLimits <- data.table::melt(data.table::setDT(rngLimits),
                                  id.vars = c("n"), variable.name = "limit")
  rngLimits$limit_id <- paste0("ctrl_",sub("[a-z]*_","",rngLimits$limit))

  if(identify[1] == "outliers") {
    outlierVars <- names(x$results)[grep(pattern = "inside",x = names(x$results))]
    outlierRows <- rowSums(!x$results[outlierVars,drop=FALSE])
    x$results <- x$results[outlierRows,]
  } else if(identify[1] != "all") {
    stopifnot(identify %in% x$results$id)
    x$results <- x$results[x$results$id %in% identify,]
  }
  # the plot
  if (x$control$standardised == FALSE) {
    pp <- ggplot2::ggplot(data=x$results,ggplot2::aes_string(x="n",y="prop_adj")) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = x$target) +
      ggplot2::geom_line(data=rngLimits,
                         ggplot2::aes_string(x="n",y="value",group="limit",col="limit_id"))
    if(label[1] == "outliers") {
      outlierVars <- names(x$results)[grep(pattern = "inside",x = names(x$results))]
      outlierRows <- (rowSums(!x$results[,outlierVars,drop=FALSE]) >= 1)
      labelData <- x$results[outlierRows,]
      pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y="prop_adj",label="id"),size=4)
    } else if(label[1] != "none") {
      labelData <- x$results[x$results$id %in% label,]
      pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y="prop_adj",label="id"),size=4)
    }
  } else if (x$control$standardised == TRUE) {
    pp <- ggplot2::ggplot(data=x$results,ggplot2::aes_string(x="n",y="observed_expected")) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 1) +
      ggplot2::geom_line(data=rngLimits,
                         ggplot2::aes_string(x="n",y="value",group="limit",col="limit_id"))
    if(label[1] == "outliers") {
      outlierVars <- names(x$results)[grep(pattern = "inside",x = names(x$results))]
      outlierRows <- (rowSums(!x$results[,outlierVars,drop=FALSE]) >= 1)
      labelData <- x$results[outlierRows,]
      pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y="observed_expected",label="id"),size=4)
    } else if(label[1] != "none") {
      labelData <- x$results[x$results$id %in% label,]
      pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y="observed_expected",label="id"),size=4)
    }
  }
  pp
}

