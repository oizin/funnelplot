#' helper function for plotting the funnel plots
#'
#' @param funnelRes funnel plot object
#'
#'
plotLimits <- function(funnelRes,...) UseMethod("plotLimits")

#' helper function for plotting the funnel plots
#'
#' @param funnelRes funnel plot object
#' @param lengthOut 500
#'
plotLimits.funnelRes <- function(funnelRes,lengthOut) {

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


#' plot funnel plot
#'
#' @param funnelRes funnel plot object
#' @param identify all
#' @param label none
#' @param lengthOut 500
#'
plot.funnelRes <- function(funnelRes,identify="all",label="none",lengthOut=500,...) {

  # calculate control limits for plotting
  rngLimits  <- plotLimits(funnelRes,lengthOut=lengthOut)

  # long form
  rngLimits <- data.table::melt(data.table::setDT(rngLimits),
                                  id.vars = c("n"), variable.name = "limit")
  rngLimits$limit_id <- paste0("ctrl_",sub("[a-z]*_","",rngLimits$limit))

  if(identify[1] == "outliers") {
    outlierVars <- names(funnelRes$results)[grep(pattern = "inside",x = names(funnelRes$results))]
    outlierRows <- rowSums(!funnelRes$results[outlierVars,drop=FALSE])
    funnelRes$results <- funnelRes$results[outlierRows,]
  } else if(identify[1] != "all") {
    stopifnot(identify %in% funnelRes$results$id)
    funnelRes$results <- funnelRes$results[funnelRes$results$id %in% identify,]
  }
  # the plot
  pp <- ggplot2::ggplot(data=funnelRes$results,ggplot2::aes_string(x="n",y="prop_adj")) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = funnelRes$target) +
    ggplot2::geom_line(data=rngLimits,
                       ggplot2::aes_string(x="n",y="value",group="limit",col="limit_id"))
  if(label[1] == "outliers") {
    outlierVars <- names(funnelRes$results)[grep(pattern = "inside",x = names(funnelRes$results))]
    outlierRows <- (rowSums(!funnelRes$results[,outlierVars,drop=FALSE]) >= 1)
    labelData <- funnelRes$results[outlierRows,]
    pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y="prop_adj",label="id"),size=4)
  } else if(label[1] != "none") {
    stopifnot(label %in% funnelRes$results$id)
    labelData <- funnelRes$results[funnelRes$results$id %in% label,]
    pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y="prop_adj",label="id"),size=4)
  }
  pp
}

