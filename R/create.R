#########################################################################
# Funnel plot functions
#
#########################################################################



#' Dispersion parameters
#'
#' @param pfunnel
#' @param w
#'
#'
dispersion <- function(funnelData,trim=NULL) UseMethod("dispersion")
dispersion.funnelData <- function(funnelData,trim=NULL) {

  #
  data <- funnelData$data
  target <- funnelData$target

  # z scores
  zz <- (data$prop_adj - target)/sqrt(data$var)
  nn <- length(zz)

  # chi-square test for over-dispersion
  inflationFactorSq <- chi <- (1/nn)*sum(zz^2)
  prob <- 2*pchisq(chi,1,lower.tail = FALSE)

  # winsorisation
  if(!is.null(trim)) {
    stopifnot(trim < 1 & trim >= 0)
    zz <- sort(zz)
    upper <- ceiling(nn*(1-trim))
    lower <- floor(nn*(trim))
    zz[zz > zz[upper]] <- zz[upper]
    zz[zz < zz[lower]] <- zz[lower]
    inflationFactorSq <- (1/nn)*sum(zz^2)
  }

  # method of moments RE approach
  pp <- data$prop_obs
  vv <- pp*(1 - pp)*(1/nn)
  ww <- 1/vv
  effectVar <- (nn*inflationFactorSq - (nn - 1)) / (sum(ww) - (sum(ww^2)/sum(ww)))

  # return
  c(stat = chi, pval = prob, inflationFactor = sqrt(inflationFactorSq),
    effectVar = effectVar)
}




#' Create the funnel plot
#'
#' @param formula a two-sided formula y ~ x1 + x2 | group, if no covariates use y ~ 1 | group
#'
#' @param target a list of values that define how the funnel plot should be constructed.
#' See pointTarget and distTarget.
#' @param data a data frame containing the variables named in formula
#'
#' CHANGE normal_approx TO METHOD = "exact", "normal" AND SO ON!
#'
funnel <- function(formula, target=pointTarget(), data) {
  ## checks
  stopifnot(all.vars(formula) %in% c(names(data)))

  ## edit formula
  sepFormula <- getFunnelFormula(formula)
  groupVar <- sepFormula$group
  outcomeVar <- sepFormula$outcome
  newFormula <- sepFormula$newForm

  ## casemix adjustment
  adjMod <- glm(newFormula,family=binomial(link="logit"),data=data)
  expected <- fitted(adjMod)
  observed <- data[[outcomeVar]]
  group <- data[[groupVar]]

  ## calculation of points for funnel
  funnelData <- groupOutcomes(observed,expected,group)

  # calculate
  # 1. overdispersion relative to assumed distribution
  # 2. parameters for random effects distribution
  dispRes <- dispersion(funnelData, trim = target$trim)

  if (target$method == "point") {
    # are we controlling for over-dispersion?
    if(target$crtlOverDisp == TRUE){
      inflationFactor <- dispRes["inflationFactor"]
    } else {
      inflationFactor <- 1
    }
    ctrlLimits <- pointLimits(funnelData,target$limits,inflationFactor,target$normalApprox)
    lower <- ctrlLimits$lower
    upper <- ctrlLimits$upper

  } else if (target$method == "distribution") {
    inflationFactor <- odRes["effectVar"]
    ctrlLimits <- randomEffectLimits(funnelData,effectVar)
    lower <- ctrlLimits$lower
    upper <- ctrlLimits$upper
  }

  # determine clinics inside/outside limits
  incontrol <- vector(mode = "list",length=length(target$limits))
  for(ii in 1:length(target$limits)) {
    incontrol[[ii]] <- (funnelData$data$prop_obs <= upper[[ii]] &
                          funnelData$data$prop_obs >= lower[[ii]])
  }
  incontrol <- data.frame(incontrol)
  tmp <- as.character(100* (1-(target$limits)))
  names(incontrol) <- paste0("inside",tmp)

  ## output
  results <- data.frame(funnelData$data,
                        upper,
                        lower,
                        incontrol,
                        row.names=NULL)
  control_limits = list(method=target$method, limits=target$limits)
  # add class
  out <- list(results = results, target = target,
              dispersion = odRes,
              control_limits = control_limits
              #adj_model = adj_mod,
              #casemix_adj = casemix_adj,
              #outcome = outcome
  )
  class(out) <- c(class(out), "funnelRes")

  # return
  out
}

pointLimits <- function(funnelData,limits,inflationFactor,normalApprox) {
  # vectors to store results
  upper <- vector(mode = "list",length = length(limits))
  lower <- vector(mode = "list",length = length(limits))

  for(ii in 1:length(limits)) {
    if(normalApprox == TRUE) {
      # distribution quantiles
      zz_lower <- qnorm(limits[ii]/2)
      zz_upper <- qnorm(1-(limits[ii]/2))
      # control limits
      upper[[ii]] <- funnelData$target + zz_upper*inflationFactor*sqrt(funnelData$grouped$var)
      lower[[ii]] <- funnelData$target + zz_lower*inflationFactor*sqrt(funnelData$grouped$var)
    } else if (normalApprox == FALSE) {
      # continuity adjusted control limits
      upper[[ii]] <- cont_adjust_bin(1-(limits[ii]/2),funnelData$grouped$n,funnelData$target)
      lower[[ii]] <- cont_adjust_bin((limits[ii]/2),funnelData$grouped$n,funnelData$target)
      # over-dispersion adjustment
      upper[[ii]] <- funnelData$target + inflationFactor*(upper[[ii]] - funnelData$target)
      lower[[ii]] <- funnelData$target + inflationFactor*(lower[[ii]] - funnelData$target)
    }
  }
  limitChr <- as.character(100 * (1 - limits))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper", limitChr)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower", limitChr)

  # return
  list(lower = lower, upper = upper)
}

#'
#'
#'
#'
randomEffectLimits <- function(funnelData,limits,effectVar) {

  # vectors to store results
  upper <- vector(mode = "list",length = length(limits))
  lower <- vector(mode = "list",length = length(limits))

  zz_lower <- qnorm(limits[ii]/2)
  zz_upper <- qnorm(1-(limits[ii]/2))

  upper[[ii]] <- funnelData$target + zz_upper*sqrt(funnelData$grouped$var + effectVar)
  lower[[ii]] <- funnelData$target + zz_lower*sqrt(funnelData$grouped$var + effectVar)

  tmp <- as.character(100* (1-(limits)))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper_",tmp)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower_",tmp)


  list(lower = lower, upper = upper)


}



#' Calculate proportion of
#'
#'
groupOutcomes <- function(observed,expected,group,target=NULL) {
  # checks
  stopifnot(is.factor(id))
  stopifnot(is.numeric(observed))
  stopifnot(is.numeric(expected)|is.null(expected))

  # calculate target
  if(is.null(target)) {
    target <- sum(observed)/length(observed)
  }

  # outcome calculations by clinic
  nn <- tapply(observed, id, length)
  obs_grpd <- tapply(observed, id, sum)
  prop_obs <- (obs_grpd/nn)
  std_err <- sqrt((target*(1-target))/nn)

  # prepare output
  adj_grpd <- tapply(expected, id, sum)
  prop_adj <- (obs_grpd/adj_grpd)*target  # standardised proportion
  data <- data.frame(id = levels(id), n = nn,
                        observed = obs_grpd, expected = adj_grpd,
                        prop_obs = prop_obs, prop_adj = prop_adj,
                        std_err = std_err, row.names = NULL)


  # add class
  out <- list(data = data, target = target)
  class(out) <- c(class(out), "funnelData")

  # return
  out
}


#' print function
#'
#'
print.pfunnel <- function(pfunnel) {
  print(pfunnel$grouped)
}

#' summary function
summary.pfunnel <- function(pfunnel) {

  # number of institutions
  n <- nrow(pfunnel$grouped)
  cat("n:",n,"\n")

  # number of clinics outside various limits
  ctrl_names <- grep("inside",names(pfunnel$grouped),value=TRUE)
  print_names <- sub("inside_","",ctrl_names)
  for(ii in 1:length(ctrl_names)) {
    print_val <- paste0("outside ",print_names[ii],"% control limits:")
    out <- sum(pfunnel$grouped[ctrl_names[ii]] == FALSE)
    out <- sum(pfunnel$grouped[ctrl_names[ii]] == FALSE)

    cat(print_val,out,"\n")
  }
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



