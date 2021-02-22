#########################################################################
# Funnel plot functions
# - Create the data for a funnel plot
#########################################################################


#' Dispersion parameters
#'
#' Calculate the inflation parameters for situations where overdispersion is present.
#'
#' @param funnelData funnel plot data
#' @param trim winsorisation
#'
dispersion <- function(funnelData,trim=NULL) {

  # move up one level
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


#' Casemix adjustment
#'
#' @param formula model formula in the form outcome ~ covariates
#' @param data dataset
#' @param params adjParams
#'
calcExpected <- function(formula,data,params) {

  if (params$method == "use_all") {
    adjMod <- do.call(params$model, list(data = data,formula = formula))
    expected <- predict(adjMod,newdata=data,type="response")
    evaluation <- NA
  } else if (params$method == "out_of_fold") {
    # make folds
    N <- nrow(data)
    per_samp <- floor(N/params$nfolds)
    samp_i <- 1:N
    avail <- rep(TRUE,length(samp_i))
    make_folds <- function(x) {
      out <- sample(x[avail],size = per_samp,replace=FALSE)
      avail[out] <<- FALSE
      out
    }
    cv_leave_out_i <- replicate(n = params$nfolds, make_folds(samp_i))

    # make folds
    cv_folds <- vector(mode = "list", length = params$nfolds)
    for(leave_out in 1:ncol(cv_leave_out_i)) {
      train_i <- as.numeric(cv_leave_out_i[,-leave_out])
      test_i <- as.numeric(cv_leave_out_i[,leave_out])
      cv_folds[[leave_out]] <- list(train = train_i, test = test_i)
    }

    # evaluate fold times (cols = folds)
    brier <- numeric(length(cv_folds))
    accuracy <- numeric(length(cv_folds))
    auc_roc <- numeric(length(cv_folds))

    expected <- c()
    test_index <- c()
    outcomeVar <- all.vars(formula)[1]

    for(fold in 1:length(cv_folds)) {
      # fit model
      adj_mod_i <- do.call(params$model, list(data = data[cv_folds[[fold]]$train,],formula = formula))

      # predictions
      test_index <- c(test_index,cv_folds[[fold]]$test)
      expected_i <- predict(adj_mod_i,newdata=data[cv_folds[[fold]]$test,],type="response",formula=formula)
      expected <- c(expected,expected_i)

      # evaluation metrics
      binary_pred_i <- classify(expected_i,cutoff = 0.5)
      true_out_i <- data[cv_folds[[fold]]$test,outcomeVar]
      brier[fold] <- mean((expected_i - true_out_i)^2)
      accuracy[fold] <- accuracy(binary_pred_i,true_out_i)
      auc_roc[fold] <- auc_roc(expected_i,true_out_i)
    }
    expected <- expected[order(test_index)]
    outcome <- data[[outcomeVar]]
    no_info_rate <- max(mean(outcome),1-mean(outcome))
    evaluation <- data.frame(run = c(1:params$nfolds,"overall"),
               brier = c(brier,mean(brier)),
               accuracy = c(accuracy,mean(accuracy)),
               auc_roc = c(auc_roc,mean(auc_roc)),
               no_info_rate = c(rep(NA,params$nfolds),no_info_rate))
  }
  list(expected=expected,evaluation=evaluation)
}


#' Check the input to funnel
#'
#' Checks whether the input to the funnel function is correct.
#'
#' @param formula formula for funnel
#' @param var_names names of variables in dataset provided to funnel
#'
check_formula <- function(formula,var_names) {
  tmp <- as.list(formula)
  tmp <- as.character(tmp[[3]])
  vars_in_form <- trimws(unlist(strsplit(tmp[2],"+",fixed = TRUE)))
  if (all(vars_in_form != "1")) {
    assertthat::assert_that(length(all.vars(formula)) >= 3,
      msg = "formula must be of the form `y ~ covariates | cluster` or `y ~ 1 | cluster`")
    assertthat::assert_that(tmp[1] == "|" & all(all.vars(formula) %in% var_names),
      msg = "formula must be of the form y ~ covariates | cluster or y ~ 1 | cluster")
  } else {
    assertthat::assert_that(tmp[1] == "|" & all(vars_in_form == "1"),
      msg = "formula must be of the form y ~ covariates | cluster or y ~ 1 | cluster")
  }
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
#' @param adj a description of the algorithm that will be used in calculating the expected outcome for each observation for the purpose of casemix adjustment.
#' See adjModel() for more details.
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
#' The package is based on Spiegelhalter (2005) and Jones, Ohlssen & Spiegelhalter (2008). All errors in implementation and are the responsibility of the package author (Oisin Fitzgerald).
#'
#' @section References
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
#'
#' Spiegelhalter, D. J. (2005). Funnel plots for comparing institutional performance. Statistics in medicine, 24(8), 1185-1202.
#'
#' Jones, H. E., Ohlssen, D. I., & Spiegelhalter, D. J. (2008). Use of the false discovery rate when comparing multiple health care providers. Journal of clinical epidemiology, 61(3), 232-240.
#'
#' @export
funnelModel <- function(formula, control=pointTarget(), adj=adjParams(), data) {
  ## checks on input args
  # prelim
  var_names <- all.vars(formula)
  y_vals <- unique(data[[var_names[1]]])
  # checks
  assertthat::assert_that(is.data.frame(data))
  assertthat::assert_that(class(formula) == "formula")
  assertthat::assert_that(all(y_vals %in% c(0,1)))
  assertthat::assert_that(check_formula(formula,names(data)))

  ## edit formula
  sepFormula <- getFunnelFormula(formula)
  idVar <- sepFormula$id
  outcomeVar <- sepFormula$outcome
  newFormula <- sepFormula$newForm

  ## casemix adjustment
  model_res <- calcExpected(newFormula,data,adj)
  expected <- model_res$expected
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
    hypTestRes <- rep(NA,nrow(funnelData$data))
  }

  # determine clinics inside/outside limits
  incontrol <- vector(mode = "list",length=length(control$limits))
  for(ii in 1:length(control$limits)) {
    if (control$standardised == FALSE) {
      incontrol[[ii]] <- (funnelData$data$prop_adj <= upper[[ii]] &
          funnelData$data$prop_adj >= lower[[ii]])
    } else if (control$standardised == TRUE) {
      incontrol[[ii]] <- (funnelData$data$observed_expected <= upper[[ii]] &
          funnelData$data$observed_expected >= lower[[ii]])
    }
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
              data = data,
              adj_model_perf = model_res$evaluation
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
#' @param ... additional args
#'
bonferroni <- function(limits, N,...) {
  new_limits <- limits/N
  new_limits
}

#' Calculate control limits
#'
#' Calculate the control limits for use in plotting and comparison of individual clusters (institutions) to a point target.
#'
#' @param target institution target value (a proportion)
#' @param n precision values at which to calculate limit (vector)
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
        # expected outcome
        e0 <- n*target
        stdErr0 <- sqrt(1/e0 - 1/n)  # 1/n adjustment
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
    if (control$standardised == TRUE) {
      zzlower <- stats::qnorm(control$limits[ii]/2)
      zzupper <- stats::qnorm(1-(control$limits[ii]/2))
      e0 <- n*target
      var0 <- (1/e0 - 1/n)
      upper[[ii]] <- 1 + zzupper*sqrt(var0 + effectVar/target^2)
      lower[[ii]] <- 1 + zzlower*sqrt(var0 + effectVar/target^2)
    } else if (control$standardised == FALSE) {
      zzlower <- stats::qnorm(control$limits[ii]/2)
      zzupper <- stats::qnorm(1-(control$limits[ii]/2))
      var0 <- target*(1-target)*(1/n)
      upper[[ii]] <- target + zzupper*sqrt(var0 + effectVar)
      lower[[ii]] <- target + zzlower*sqrt(var0 + effectVar)
    }
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
    criticalTmp <- pval[pval < criticalPEach]
    if (length(criticalTmp) == 0) {
      criticalPAll <- limits[alpha]
    } else {
      criticalPAll <- max(criticalTmp)
    }
    out[alpha] <- criticalPAll
  }
  out
}

