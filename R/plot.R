#' helper function for plotting the funnel plots
#'
#' @export
plotLimits <- function(funnelRes,...) UseMethod("plotLimits")

#'
#'
#'
#' @export
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
#' @export
plot.funnelRes <- function(funnelRes,identify="all",lengthOut=500,...) {

  # calculate control limits for plotting
  rngLimits  <- plotLimits(funnelRes,lengthOut=lengthOut)

  # long form
  rngLimits <- data.table::melt(data.table::setDT(rngLimits),
                                  id.vars = c("n"), variable.name = "limit")
  rngLimits$limit_id <- paste0("ctrl_",sub("[a-z]*_","",rngLimits$limit))

  if(identify[1] != "all") {
    stopifnot(identify %in% funnelRes$results$id)
    funnelRes$results <- funnelRes$results[funnelRes$results$id %in% identify,]
  }

  # the plot
  ggplot2::ggplot(data=funnelRes$results,ggplot2::aes(x=n,y=prop_adj)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = funnelRes$target) +
    ggplot2::geom_line(data=rngLimits,
                       ggplot2::aes(x=n,y=value,group=limit,col=limit_id))
}

