#' helper function for plotting the funnel plots
#'
#'
ctrl_limits <- function(pfunnel,...) UseMethod("ctrl_limits")
ctrl_limits.pfunnel <- function(pfunnel,length_out,...) {

  grouped <- pfunnel$grouped

  # limits
  mn <- min(pfunnel$grouped$n)
  mx <- max(pfunnel$grouped$n)
  rng <- floor(seq(from = mn,to = mx,length.out=length_out))
  limits <- pfunnel$control_limits$limits

  # inflation factor
  if(pfunnel$over_dispersion$control == TRUE) {
    inflation_factor <- pfunnel$over_dispersion$stats["inflation_factor"]
  } else {
    inflation_factor <- 1
  }

  # limits
  upper <- vector(mode = "list",length = length(limits))
  lower <- vector(mode = "list",length = length(limits))

  # calculate control limits
  for(ii in 1:length(limits)) {
    if(pfunnel$control_limits["method"] == "normal") {
      # distribution quantiles
      zz_lower <- qnorm(limits[ii]/2)
      zz_upper <- qnorm(1-(limits[ii]/2))
      # control limits
      if (pfunnel$over_dispersion$random_effects == TRUE) {
        tau2 <- pfunnel$over_dispersion$stats["tau2"]
        vv <- (pfunnel$target*(1-pfunnel$target))/rng
        upper[[ii]] <- pfunnel$target + zz_upper*sqrt(vv + tau2)
        lower[[ii]] <- pfunnel$target + zz_lower*sqrt(vv + tau2)
      } else {
        vv <- (pfunnel$target*(1-pfunnel$target))/rng
        upper[[ii]] <- pfunnel$target + zz_upper*inflation_factor*sqrt(vv)
        lower[[ii]] <- pfunnel$target + zz_lower*inflation_factor*sqrt(vv)
      }
    } else if(pfunnel$control_limits["method"] == "exact") {
      # continuity adjusted control limits
      upper[[ii]] <- cont_adjust_bin(1-(limits[ii]/2),rng,pfunnel$target)
      lower[[ii]] <- cont_adjust_bin((limits[ii]/2),rng,pfunnel$target)
      # over-dispersion adjustment
      upper[[ii]] <- pfunnel$target + inflation_factor*(upper[[ii]] - pfunnel$target)
      lower[[ii]] <- pfunnel$target + inflation_factor*(lower[[ii]] - pfunnel$target)
    }
  }
  tmp <- as.character(100 * (1-(limits)))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper_",tmp)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower_",tmp)

  ctrl_limits <- data.frame(n = rng, upper, lower)
  ctrl_limits
}


#' plot funnel plot
#'
#'
plot.pfunnel <- function(pfunnel,identify="all",length_out=500,...) {

  # calculate control limits for plotting
  ctrl_limits <- ctrl_limits(pfunnel,length_out=length_out)

  # long form
  ctrl_limits <- data.table::melt(data.table::setDT(ctrl_limits),
                                  id.vars = c("n"), variable.name = "limit")
  ctrl_limits$limit_id <- paste0("ctrl_",sub("[a-z]*_","",ctrl_limits$limit))

  if(identify[1] != "all") {
    stopifnot(identify %in% pfunnel$grouped$group)
    pfunnel$grouped <- pfunnel$grouped[pfunnel$grouped$group %in% identify,]
  }

  # the actual plot
  ggplot2::ggplot(data=pfunnel$grouped,ggplot2::aes(x=n,y=prop_obs)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = pfunnel$target) +
    ggplot2::geom_line(data=ctrl_limits,
                       ggplot2::aes(x=n,y=value,group=limit,col=limit_id))
}

#' plot the funnel shape
#'
#'
ggfunnel <- function(pfunnel,...) UseMethod("ggfunnel")
ggfunnel.pfunnel <- function(pfunnel,length_out=500,...) {

  # calculate control limits for plotting
  ctrl_limits <- ctrl_limits(pfunnel,length_out=length_out)

  # long form
  ctrl_limits <- data.table::melt(data.table::setDT(ctrl_limits),
                                  id.vars = c("n"), variable.name = "limit")
  ctrl_limits$limit_id <- paste0(sub("[a-z]*_","",ctrl_limits$limit),"% limits")

  # the actual plot
  ggplot2::ggplot(data=ctrl_limits,ggplot2::aes(x=n,y=value,group=limit,col=limit_id)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(title = NULL))
}

