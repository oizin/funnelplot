#########################################################################
# Funnel plot functions
#
#########################################################################


#' Calculate proportion of
#'
#'
grouped_proportion <- function(observed,group,expected=NULL,target=NULL) {
  # checks
  stopifnot(is.factor(group))
  stopifnot(is.numeric(observed))
  stopifnot(is.numeric(expected)|is.null(expected))

  # calculate target
  if(is.null(target)) {
    target <- sum(observed)/length(observed)
  }

  # outcome calculations by clinic
  nn <- tapply(observed,group,length)
  obs_grpd <- tapply(observed,group,sum)
  prop_obs <- (obs_grpd/nn)
  vv <- (target*(1-target))/nn

  # prepare output
  if(!is.null(expected)) {
    adj_grpd <- tapply(expected,group,sum)
    prop_adj <- (obs_grpd/adj_grpd)*target  # standardised proportion
    grouped <- data.frame(group = levels(group), n = nn,
                          observed = obs_grpd, expected = adj_grpd,
                          prop_obs = prop_obs, prop_adj = prop_adj,
                          var = vv,row.names = NULL)
  } else {
    grouped <- data.frame(group = levels(group), n = nn,
                          observed = obs_grpd,
                          prop_obs = prop_obs,
                          var = vv,row.names = NULL)
  }

  # add class
  out <- list(grouped = grouped, target = target)
  class(out) <- c(class(out),"pfunnel")

  # return
  out
}


#'
#'
#'
#'
z_scores <- function(pfunnel,...) UseMethod("z_scores")
z_scores.pfunnel <- function(pfunnel,type="naive") {
  # checks
  stopifnot(is.list(pfunnel))
  stopifnot(names(pfunnel) %in% c("grouped","target"))
  stopifnot(type %in% c("naive","adjusted"))

  grouped <- pfunnel$grouped
  stopifnot(names(grouped) %in% c("group","n","observed","prop_obs","var"))
  target <- pfunnel$target

  # calculations
  if(type == "naive") {
    zz <- (grouped$prop_obs - target)/sqrt(grouped$var)
  } else {
    zz <- (grouped$prop_adj - target)/sqrt(grouped$var)
  }

  # return
  zz
}


#' Check for over-dispersion
#'
#'
#'
over_dispersion <- function(pfunnel,w=NULL) UseMethod("over_dispersion")
over_dispersion.pfunnel <- function(pfunnel,w=NULL) {

  grouped <- pfunnel$grouped
  target <- pfunnel$target

  # z-scores
  if(is.null(pfunnel$grouped$prop_adj)) {
    zz <- (grouped$prop_obs - target)/sqrt(grouped$var)
  } else {
    zz <- (grouped$prop_adj - target)/sqrt(grouped$var)
  }
  nn <- length(zz)

  # winsorisation
  if(!is.null(w)) {
    stopifnot(w < 1 & w >= 0)
    zz <- sort(zz)
    upper <- ceiling(nn*(1-w))
    lower <- floor(nn*(w))
    zz[zz > zz[upper]] <- zz[upper]
    zz[zz < zz[lower]] <- zz[lower]
  }

  # chi-square naive
  chi <- (1/nn)*sum(zz^2)
  prob <- 2*pchisq(chi,1,lower.tail = FALSE)

  # return
  c(stat = chi, pval = prob, inflation_factor = sqrt(chi))
}

#' Create the funnel plot
#'
#'
funnel <- function(outcome,group,limits=0.05,normal_approx=FALSE,casemix_adj=NULL,
                   w=NULL,target=NULL,ctrl_over_disp=TRUE) {
  # checks
  stopifnot(is.factor(group))
  stopifnot(is.numeric(outcome))
  stopifnot(is.data.frame(casemix_adj)|is.null(casemix_adj)|is.matrix(casemix_adj))

  ## casemix adjustment
  if(!is.null(casemix_adj)) {
    if(is.data.frame(casemix_adj)) {
      stopifnot(!is.na(casemix_adj))
      nms <- names(casemix_adj)
      form <- as.formula(paste0("~",paste0(nms,collapse="+")))
      mm <- model.matrix(form,data=casemix_adj)
    } else {
      mm <- casemix_adj
    }
    adj_mod <- glm.fit(x = mm, y = outcome,family=binomial(link="logit"))
    expected <- fitted(adj_mod)
    pfunnel <- grouped_proportion(outcome,group,expected,target)
  } else {
    pfunnel <- grouped_proportion(outcome,group,target=target)
  }

  # calculate proportions
  target <- pfunnel$target
  grouped <- pfunnel$grouped

  # inflation factor
  od_res <- over_dispersion(pfunnel,w=w)
  if(ctrl_over_disp == TRUE){
    inflation_factor <- od_res["inflation_factor"]
  } else {
    inflation_factor <- 1
  }

  if(normal_approx == TRUE) {
    # distribution quantiles
    zz_lower <- qnorm(limits/2)
    zz_upper <- qnorm(1-(limits/2))
    # control limits
    upper <- target + zz_upper*inflation_factor*sqrt(grouped$var)
    lower <- target + zz_lower*inflation_factor*sqrt(grouped$var)
  } else {
    # continuity adjusted control limits
    upper <- cont_adjust_bin(1-(limits/2),pfunnel$grouped$n,target)
    lower <- cont_adjust_bin((limits/2),grouped$n,target)
    # over-dispersion adjustment
    upper <- pfunnel$target + inflation_factor*(upper - pfunnel$target)
    lower <- pfunnel$target + inflation_factor*(lower - pfunnel$target)
  }

  # determine clinics inside/outside limits
  incontrol <- (pfunnel$grouped$prop_obs <= upper &
                  pfunnel$grouped$prop_obs >= lower)

  # output
  grouped <- data.frame(grouped,upper=upper,lower=lower,incontrol=incontrol,
                        row.names=NULL)
  control_limits = c(normal_approx=normal_approx, limits=limits)

  # add class
  out <- list(grouped = grouped, target = target,
              over_dispersion = list(control = ctrl_over_disp, stats = od_res),
              control_limits = control_limits,
              winsorisation = w)
  class(out) <- c(class(out),"pfunnel")

  # return
  out
}

#' helper function for plotting the funnel plots
#'
#'
ctrl_limits <- function(pfunnel) UseMethod("ctrl_limits")
ctrl_limits.pfunnel <- function(pfunnel) {

  grouped <- pfunnel$grouped

  # limits
  mn <- min(pfunnel$grouped$n)
  mx <- max(pfunnel$grouped$n)
  rng <- floor(seq(from = mn,to = mx,length.out = 500))
  limits <- pfunnel$control_limits["limits"]

  # inflation factor
  if(pfunnel$over_dispersion$control == TRUE) {
    inflation_factor <- pfunnel$over_dispersion$stats["inflation_factor"]
  } else {
    inflation_factor <- 1
  }

  # calculate control limits
  if(pfunnel$control_limits["normal_approx"] == TRUE) {
    # normal approximation control limits
    # distribution quantiles
    zz_lower <- qnorm(limits/2)
    zz_upper <- qnorm(1-(limits/2))
    # control limits
    vv <- (pfunnel$target*(1-pfunnel$target))/rng
    upper <- pfunnel$target + zz_upper*inflation_factor*sqrt(vv)
    lower <- pfunnel$target + zz_lower*inflation_factor*sqrt(vv)
  } else {
    # continuity adjusted binomial control limits
    upper <- cont_adjust_bin(1-(limits/2),rng,pfunnel$target)
    lower <- cont_adjust_bin((limits/2),rng,pfunnel$target)
    # inflation adjustment
    upper <- pfunnel$target + inflation_factor*(upper - pfunnel$target)
    lower <- pfunnel$target + inflation_factor*(lower - pfunnel$target)
  }
  ctrl_limits <- data.frame(n = rng, upper = upper, lower = lower)
  ctrl_limits
}


#' plot funnel plot
#'
#'
ggplot <- ggplot2::ggplot
ggplot.pfunnel <- function(pfunnel,...) {

  # calcalate control limits for plotting
  ctrl_limits <- ctrl_limits(pfunnel)

  # the actual plot
  ggplot2::ggplot(data=pfunnel$grouped,ggplot2::aes(x=n,y=prop_obs)) +
    ggplot2::geom_point(ggplot2::aes(col=incontrol)) +
    ggplot2::geom_hline(yintercept = pfunnel$target) +
    ggplot2::geom_line(data=ctrl_limits,ggplot2::aes(x=n,y=lower)) +
    ggplot2::geom_line(data=ctrl_limits,ggplot2::aes(x=n,y=upper))
}

#'
#'
#'
#'
cont_adjust_bin <- function(limit, n, target) {
  rp <- qbinom(p = limit,size = n,prob = target)
  t1 <- pbinom(q = rp,size = n,prob = target)
  t2 <- dbinom(x = rp,size = n,prob = target)
  aa <- (t1-limit)/t2
  ctrl_limit <- (rp-aa)/n
  ctrl_limit
}



