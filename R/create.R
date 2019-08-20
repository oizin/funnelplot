#########################################################################
# Funnel plot functions
#
#########################################################################

#' Dispersion parameters
#'
#' Calculate the inflation parameters for situations where overdispersion is present.
#'
#' @param funnelData funnel plot data
#' @param trim winsorisation
#'
dispersion <- function(funnelData,trim=NULL) UseMethod("dispersion")

#' Dispersion parameters
#'
#' Calculate the inflation parameters for situations where overdispersion is present.
#'
#' @param funnelData funnel plot data
#' @param trim winsorisation
#'
dispersion.funnelData <- function(funnelData,trim=NULL) {

  #
  data <- funnelData$data
  target <- funnelData$target

  # z scores
  zz <- (data$prop_adj - target)/data$std_err0
  N <- length(zz)

  # chi-square test for over-dispersion
  inflationFactorSq <- chi <- (1/N)*sum(zz^2)
  prob <- 2*stats::pchisq(chi,1,lower.tail = FALSE)

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

#' Check the input to funnel
#'
#' Checks whether the input to the funnel function is correct.
#'
#' @param x formula for funnel
#' @param var_names names of variables in dataset provided to funnel
#'
check_formula <- function(x,var_names) {
  tmp <- as.list(x)
  tmp <- as.character(tmp[[3]])
  vars_in_form <- trimws(unlist(strsplit(tmp[2],"+",fixed = TRUE)))
  assertthat::assert_that(tmp[1] == "|")
  assertthat::assert_that(vars_in_form == "1" | all(vars_in_form %in% var_names))
}

#' Table of results for risk adjusted funnel plots.
#'
#' Perform risk adjusted cluster (institution) comparison for the purpose of outlier discovery.
#'
#' @param formula a two-sided formula. Either of the form y ~ x1 + x2 | group,
#' or if no covariates use y ~ 1 | group
#' @param data a data frame containing the variables named in formula
#' @param control a description of how the funnel plot target and control limits should be constructed. See the control functions
#' pointTarget() and distTarget().
#'
#' @return funnel returns an object of class "funnelRes". The functions \code{outliers} and \{plot} are used to obtain a summary of the
#' outlying clusters (institutions) and produce a funnel plot.\cr
#' An object of class "funnelRes" is a list containing at least the following components:\cr
#'
#'  \itemize{
#'  \item{\code{results}:}{ a data frame with each row summarising the performance of a cluster (institution) on the metric of interest}
#'  \item{\code{dispersion}:}{ the dispersion statistic and p-value, and over-inflation factor (on the standard error scale)}
#'  \item{\code{target}:}{ the average value of the performance metric in the population}
#' }
#'
#'
#' @section Author(s):
#' The package is based on Spiegelhalter (2005) and Jones, Ohlssen & Spiegelhalter (2008). All errors in implementation and are the responsibility of the package authors (Oisin Fitzgerald).
#'
#' @section References:
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
#'
#' Spiegelhalter, D. J. (2005). Funnel plots for comparing institutional performance. Statistics in medicine, 24(8), 1185-1202.
#'
#' Jones, H. E., Ohlssen, D. I., & Spiegelhalter, D. J. (2008). Use of the false discovery rate when comparing multiple health care providers. Journal of clinical epidemiology, 61(3), 232-240.
#'
#' @export
funnel <- function(formula, control=pointTarget(), data) {
  ## checks on input args
  # prelim
  var_names <- all.vars(formula)
  y_vals <- unique(data[[var_names[1]]])
  # checks
  assertthat::assert_that(is.data.frame(data))
  assertthat::assert_that(class(formula) == "formula")
  assertthat::assert_that(all(var_names %in% c(names(data))))
  assertthat::assert_that(all(y_vals %in% c(0,1)))
  assertthat::assert_that(check_formula(formula,var_names))

  ## edit formula
  sepFormula <- getFunnelFormula(formula)
  idVar <- sepFormula$id
  outcomeVar <- sepFormula$outcome
  newFormula <- sepFormula$newForm

  ## casemix adjustment
  adjMod <- stats::glm(newFormula,family=stats::binomial(link="logit"),data=data)
  expected <- stats::fitted(adjMod)
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
                              stdErr = inflationFactor*funnelData$data$std_err0)
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
  tmp <- as.character(100*(1-(control$limits)))
  names(incontrol) <- paste0("inside",tmp)

  ## output
  results <- data.frame(funnelData$data,
                        upper,
                        lower,
                        incontrol,
                        hypTestRes,
                        row.names=NULL)

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

#' Bonferroni multiplicity adjustment.
#'
#' Adjusts control limits to correct for multiple comparisons to a target.
#'
#' @param limits alpha levels
#' @param N number of tests
#'
bonferroni <- function(limits, N) {
  new_limits <- limits/N
  new_limits
}

#' Calculate control limits
#'
#' Calculate the control limits for use in plotting and comparison of individual clusters (institutions) to a point target.
#'
#' @param target institution target value (a proportion)
#' @param n precision values at which to calculate limit
#' @param N number of clusters (institutions)
#' @param inflationFactor Unexplained variation
#' @param control description of task
#' @param pval p-values for FDR calculations
#'
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
    if (control$standardised == TRUE) {
      if (control$normalApprox == TRUE) {
        # expected births
        e0 <- n*target
        stdErr0 <- sqrt(1/e0)
        # distribution quantiles
        zzlower <- stats::qnorm(control$limits[ii]/2)
        zzupper <- stats::qnorm(1-(control$limits[ii]/2))
        # control limits
        lower[[ii]] <- 1 + zzlower*inflationFactor*stdErr0
        upper[[ii]] <- 1 + zzupper*inflationFactor*stdErr0
      } else if (control$normalApprox == FALSE) {
        stop("only normal approx currently implemented for standardised rates")
      }
    } else if (control$standardised == FALSE) {
      if (control$normalApprox == TRUE) {
        # standard error for approximations
        stdErr0 <- sqrt(target*(1-target)*(1/n))
        # distribution quantiles
        zzlower <- stats::qnorm(control$limits[ii]/2)
        zzupper <- stats::qnorm(1-(control$limits[ii]/2))
        # control limits
        lower[[ii]] <- target + zzlower*inflationFactor*stdErr0
        upper[[ii]] <- target + zzupper*inflationFactor*stdErr0
      } else if (control$normalApprox == FALSE) {
        # continuity adjusted control limits
        upper[[ii]] <- contAdjustBin(1-(control$limits[ii]/2),n,target)
        lower[[ii]] <- contAdjustBin((control$limits[ii]/2),n,target)
        # over-dispersion adjustment
        upper[[ii]] <- target + inflationFactor*(upper[[ii]] - target)
        lower[[ii]] <- target + inflationFactor*(lower[[ii]] - target)
      }
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

#' Calculate control limits
#'
#' Calculate the control limits for use in plotting and comparison of individual clusters (institutions) to a distribution target.
#'
#' @param target cluster (institution) mean distribution mean target value
#' @param effectVar cluster (institution) mean distribution variance
#' @inheritParams pointLimits
#'
randomEffectLimits <- function(target,n,effectVar,control) {

  # vectors to store results
  upper <- vector(mode = "list",length = length(control$limits))
  lower <- vector(mode = "list",length = length(control$limits))

  for(ii in 1:length(control$limits)) {
    zzlower <- stats::qnorm(control$limits[ii]/2)
    zzupper <- stats::qnorm(1-(control$limits[ii]/2))

    var0 <- target*(1-target)*(1/n)
    upper[[ii]] <- target + zzupper*sqrt(var0 + effectVar)
    lower[[ii]] <- target + zzlower*sqrt(var0 + effectVar)
  }

  tmp <- as.character(100*(1-(control$limits)))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper_",tmp)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower_",tmp)
  # return
  out <- list(lower = lower, upper = upper)
  out
}



#' Calculate summary values per cluster
#'
#' Calculates summary measures of per cluster (institution) performance
#'
#' @param observed observed outcome (a vector)
#' @param expected expected outcome (a vector)
#' @param id cluster (institution) identifier
#' @inheritParams pointLimits
#'
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
  exp_grpd <- tapply(expected, id, sum)
  prop_obs <- (obs_grpd/nn)
  std_err0 <- sqrt((target*(1-target))/nn)

  # prepare output
  adj_grpd <- tapply(expected, id, sum)
  prop_adj <- (obs_grpd/adj_grpd)*target  # standardised proportion
  data <- data.frame(id = levels(id), n = nn,
                        observed = obs_grpd, expected = adj_grpd,
                        observed_expected = obs_grpd/exp_grpd,
                        prop_obs = prop_obs, prop_adj = prop_adj,
                        std_err0 = std_err0, row.names = NULL)

  # add class
  out <- list(data = data, target = target)
  class(out) <- c(class(out), "funnelData")

  # return
  out
}



#' Continuity adjusted binomial limits
#'
#' Adjust binomial distribution control limits for lack of continuity. See #4 of A.1.1 of Spiegelhalter (2005).
#'
#' @param limit a numeric vector containing the (1-limit)100\% values for the control limits.
#' @param n number of clusters (institutions) being compared.
#' @inheritParams pointLimits
#'
contAdjustBin <- function(limit, n, target) {
  rp <- stats::qbinom(p = limit,size = n,prob = target)
  t1 <- stats::pbinom(q = rp,size = n,prob = target)
  t2 <- stats::dbinom(x = rp,size = n,prob = target)
  aa <- (t1-limit)/t2
  ctrl_limit <- (rp-aa)/n
  ctrl_limit
}

#' Normal assuming hypothesis test on grouped data
#'
#' Compare institutions to a target using a normality assumption.
#'
#' @param obs observed
#' @param exp null hypothesis
#' @param stdErr standard error
#'
normHypTest <- function(obs, exp, stdErr) {
  zz <- (obs - exp)/stdErr
  pval <- stats::pnorm(abs(zz),lower.tail = FALSE)*2
  data.frame(zz, pval)
}

#' False discovery rate multiplicity adjustment
#'
#' Adjust control limits to correct for multiple comparisons to a target.
#'
#' @param limits a numeric vector containing the (1-limit)100\% values for the control limits.
#' @inheritParams pointLimits
#'
fdr <- function(limits,N,pval) {
  out <- numeric(length = length(limits))
  for(alpha in 1:length(limits)) {
    orderedPos <- order(pval)
    criticalPEach <- (orderedPos*limits[alpha])/N
    criticalPAll <- max(pval[pval < criticalPEach])
    out[alpha] <- criticalPAll
  }
  out
}

