#########################################################################
# Funnel plot functions
# - Plot a funnel plot
#########################################################################

#' helper function for plotting the funnel plots
#'
#' Determines the position of the control limits at N values.
#'
#' @param funnelRes funnel plot object
#' @param lengthOut resolution of control limits. Number of points sampled to construct limit lines. default = 500.
#'
plotLimits <- function(funnelRes,lengthOut) {
  ## checks on input args
  assertthat::assert_that(all(class(funnelRes) %in% c("list","funnelRes")))
  assertthat::assert_that(is.integer(lengthOut))

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
plot.funnelRes <- function(x,identify="all",label="none",lengthOut=500L,...) {
  ## checks on input args
  assertthat::assert_that(all(class(x) %in% c("list","funnelRes")))
  assertthat::assert_that(any(c(identify %in% x$results$id, identify %in% c("all","outliers"))))
  assertthat::assert_that(any(c(label %in% x$results$id, label %in% c("none","outliers"))))
  assertthat::assert_that(is.integer(lengthOut))
  assertthat::assert_that(lengthOut > 0)
  assertthat::assert_that(length(unique(x$results$n)) > 1)

  # calculate control limits for plotting
  rngLimits  <- plotLimits(x,lengthOut=lengthOut)

  # long form
  rngLimits <- data.table::melt(data.table::setDT(rngLimits),
                                  id.vars = c("n"), variable.name = "limit")
  #rngLimits$expected <- rngLimits$n*x$target
  rngLimits$limit_id <- paste0("ctrl_",sub("[a-z]*_","",rngLimits$limit))

  if(identify[1] == "outliers") {
    outlierVars <- names(x$results)[grep(pattern = "inside",x = names(x$results))]
    outlierRows <- rowSums(!x$results[outlierVars,drop=FALSE])
    x$results <- x$results[outlierRows,]
  } else if(identify[1] != "all") {
    stopifnot(any(identify %in% x$results$id))
    x$results <- x$results[x$results$id %in% identify,]
  }
  # what to plot
  if (x$control$standardised == FALSE) {
    outcome <- "prop_adj"
    yintercept <- x$target
    precision <- "n"
  } else if (x$control$standardised == TRUE) {
    outcome <- "observed_expected"
    yintercept <- 1
    precision <- "n"
    #precision <- "expected"
  }
  # the plot
  pp <- ggplot2::ggplot(data=x$results,ggplot2::aes_string(x=precision,y=outcome)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = yintercept) +
    ggplot2::geom_line(data=rngLimits,
                       ggplot2::aes_string(x=precision,y="value",
                         group="limit",col="limit_id"))
  if(label[1] == "outliers") {
    outlierVars <- names(x$results)[grep(pattern = "inside",x = names(x$results))]
    outlierRows <- (rowSums(!x$results[,outlierVars,drop=FALSE]) >= 1)
    labelData <- x$results[outlierRows,]
    pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x="n",y=outcome,label="id"),size=4)
  } else if(label[1] != "none") {
    labelData <- x$results[x$results$id %in% label,]
    pp <- pp + ggplot2::geom_text(data=labelData,ggplot2::aes_string(x=precision,
      y=outcome,label="id"),size=4)
  }
  pp
}


#' A ggtheme
#'
#' A minimalist ggtheme for drawing funnel plots
#'
#' @param base_size base font size
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect elements
#'
#' @export
theme_funnel <- function(base_size = 11, base_family = "",
                     base_line_size = base_size / 22,
                     base_rect_size = base_size / 22) {
  # Starts with theme_grey and then modify some parts
  ggplot2::theme_grey(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) +
    ggplot2::theme(
      # white background and dark border
      panel.background = ggplot2::element_rect(fill = "gray96", colour = NA),
      legend.key       = ggplot2::element_rect(fill = "white", colour = NA),
      complete = TRUE
    )
}
