#########################################################################
# Funnel plot functions
#
#########################################################################


#' Dispersion parameters
#'
#' @param pfunnel
#' @param w
#'
#' @export
dispersion <- function(funnelData,trim=NULL) UseMethod("dispersion")

#' @export
dispersion.funnelData <- function(funnelData,trim=NULL) {

  #
  data <- funnelData$data
  target <- funnelData$target

  # z scores
  zz <- (data$prop_adj - target)/data$std_err0
  N <- length(zz)

  # chi-square test for over-dispersion
  inflationFactorSq <- chi <- (1/N)*sum(zz^2)
  prob <- 2*pchisq(chi,1,lower.tail = FALSE)

  # winsorisation
  if(!is.null(trim)) {
    stopifnot(trim < 1 & trim >= 0)
    zz <- sort(zz)
    upper <- ceiling(N*(1-trim))
    lower <- floor(N*(trim))
    zz[zz > zz[upper]] <- zz[upper]
    zz[zz < zz[lower]] <- zz[lower]
    inflationFactorSq <- (1/N)*sum(zz^2)
  }

  # method of moments RE approach
  n <- data$n
  pp <- data$prop_obs
  vv <- pp*(1 - pp)*(1/n)
  ww <- 1/vv
  effectVar <- (N*inflationFactorSq - (N - 1)) / (sum(ww) - (sum(ww^2)/sum(ww)))

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
#' @export
funnel <- function(formula, control=pointTarget(), data) {
  ## checks
  stopifnot(all.vars(formula) %in% c(names(data)))

  ## edit formula
  sepFormula <- getFunnelFormula(formula)
  idVar <- sepFormula$id
  outcomeVar <- sepFormula$outcome
  newFormula <- sepFormula$newForm

  ## casemix adjustment
  adjMod <- glm(newFormula,family=binomial(link="logit"),data=data)
  expected <- fitted(adjMod)
  observed <- data[[outcomeVar]]
  id <- as.factor(data[[idVar]])

  ## calculation of points for funnel
  funnelData <- groupOutcomes(observed,expected,id)

  # calculate
  # 1. overdispersion relative to assumed distribution
  # 2. parameters for random effects distribution
  dispRes <- dispersion(funnelData, trim = control$trim)

  if (control$method == "point") {
    # are we controlling for over-dispersion?
    if(control$crtlOverDisp == TRUE){
      inflationFactor <- dispRes["inflationFactor"]
    } else {
      inflationFactor <- 1
    }
    hypTestRes <- normHypTest(obs = funnelData$data$prop_adj,
                              exp = funnelData$target,
                              stdErr = funnelData$data$std_err0)
    ctrlLimits <- pointLimits(target=funnelData$target,
                              n=funnelData$data$n,
                              N=nrow(funnelData$data),
                              inflationFactor,
                              control,
                              pval = hypTestRes$pval)
    lower <- ctrlLimits$lower
    upper <- ctrlLimits$upper

  } else if (control$method == "distribution") {
    effectVar <- dispRes["effectVar"]
    ctrlLimits <- randomEffectLimits(target=funnelData$target,
                                     n=funnelData$data$n,
                                     effectVar,
                                     control)
    lower <- ctrlLimits$lower
    upper <- ctrlLimits$upper
  }

  # determine clinics inside/outside limits
  incontrol <- vector(mode = "list",length=length(control$limits))
  for(ii in 1:length(control$limits)) {
    incontrol[[ii]] <- (funnelData$data$prop_adj <= upper[[ii]] &
                          funnelData$data$prop_adj >= lower[[ii]])
  }
  incontrol <- data.frame(incontrol)
  tmp <- as.character(100* (1-(control$limits)))
  names(incontrol) <- paste0("inside",tmp)

  ## output
  results <- data.frame(funnelData$data,
                        upper,
                        lower,
                        incontrol,
                        hypTestRes,
                        row.names=NULL)
  control_limits = list(method=control$method, limits=control$limits)
  # add class
  out <- list(results = results,
              formula = formula,
              control = control,
              dispersion = dispRes,
              ctrlLimits = ctrlLimits,
              target = funnelData$target,
              data = data
              #adj_model = adj_mod,
              #casemix_adj = casemix_adj,
              #outcome = outcome
  )
  class(out) <- c(class(out), "funnelRes")

  # return
  out
}

#' @export
bonferroni <- function(limits, N, ...) {
  new_limits <- limits/N
  new_limits
}

#'
#'
#'
#' @export
pointLimits <- function(target,n,N,inflationFactor,control,pval) {
  # vectors to store results
  upper <- vector(mode = "list",length = length(control$limits))
  lower <- vector(mode = "list",length = length(control$limits))

  # multiplicity adjustment
  limitChr <- as.character(100 * (1 - control$limits))
  if(control$multAdj != "none") {
    control$limits <- do.call(control$multAdj,
                              args = list(limits=control$limits,N=N,pval=pval))
  }

  for(ii in 1:length(control$limits)) {
    if(control$normalApprox == TRUE) {
      # standard error for approximations
      stdErr0 <- sqrt(target*(1-target)*(1/n))
      # distribution quantiles
      zzlower <- qnorm(control$limits[ii]/2)
      zzupper <- qnorm(1-(control$limits[ii]/2))
      # control limits
      lower[[ii]] <- target + zzlower*inflationFactor*stdErr0
      upper[[ii]] <- target + zzupper*inflationFactor*stdErr0
    } else if (control$normalApprox == FALSE) {
      # continuity adjusted control limits
      upper[[ii]] <- cont_adjust_bin(1-(control$limits[ii]/2),n,target)
      lower[[ii]] <- cont_adjust_bin((control$limits[ii]/2),n,target)
      # over-dispersion adjustment
      upper[[ii]] <- target + inflationFactor*(upper[[ii]] - target)
      lower[[ii]] <- target + inflationFactor*(lower[[ii]] - target)
    }
  }
  upper <- data.frame(upper)
  names(upper) <- paste0("upper", limitChr)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower", limitChr)
  # return
  out <- list(lower = lower, upper = upper)
  out
}

#' Calculate ...
#'
#' @param target ...
#' @param n
#' @param effectVar
#' @param control
#'
#' @export
randomEffectLimits <- function(target,n,effectVar,control) {

  # vectors to store results
  upper <- vector(mode = "list",length = length(control$limits))
  lower <- vector(mode = "list",length = length(control$limits))

  zzlower <- qnorm(control$limits[ii]/2)
  zzupper <- qnorm(1-(control$limits[ii]/2))

  var0 <- target*(1-target)*(1/n)
  upper[[ii]] <- target + zzupper*sqrt(var0 + effectVar)
  lower[[ii]] <- target + zzlower*sqrt(var0 + effectVar)

  tmp <- as.character(100* (1-(control$limits)))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper_",tmp)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower_",tmp)
  # return
  out <- list(lower = lower, upper = upper)
  out
}



#' Calculate proportion of
#'
#' @export
groupOutcomes <- function(observed,expected,id,target=NULL) {
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
  std_err0 <- sqrt((target*(1-target))/nn)

  # prepare output
  adj_grpd <- tapply(expected, id, sum)
  prop_adj <- (obs_grpd/adj_grpd)*target  # standardised proportion
  data <- data.frame(id = levels(id), n = nn,
                        observed = obs_grpd, expected = adj_grpd,
                        prop_obs = prop_obs, prop_adj = prop_adj,
                        std_err0 = std_err0, row.names = NULL)

  # add class
  out <- list(data = data, target = target)
  class(out) <- c(class(out), "funnelData")

  # return
  out
}



#'
#'
#'
#' @export
cont_adjust_bin <- function(limit, n, target) {
  rp <- qbinom(p = limit,size = n,prob = target)
  t1 <- pbinom(q = rp,size = n,prob = target)
  t2 <- dbinom(x = rp,size = n,prob = target)
  aa <- (t1-limit)/t2
  ctrl_limit <- (rp-aa)/n
  ctrl_limit
}

#'
#'
#' @export
normHypTest <- function(obs, exp, stdErr) {
  zz <- (obs - exp)/stdErr
  pval <- pnorm(abs(zz),lower.tail = FALSE)*2
  data.frame(zz, pval)
}

#'
#'
#'
#'

fdr <- function(limits,N,pval) {
  out <- numeric(length = length(limits))
  for(alpha in 1:length(limits)) {
    orderedPos <- order(pval)
    criticalPEach <- (orderedPos*limits[alpha])/N
    criticalPAll <- max(pval[pval < criticalVal])
    out[i] <- criticalPAll
  }
  out
}

