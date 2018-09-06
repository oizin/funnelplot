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
#' @param pfunnel
#' @param w
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

  # chi-square test for over-dispersion
  inflation_factor2 <- chi <- (1/nn)*sum(zz^2)
  prob <- 2*pchisq(chi,1,lower.tail = FALSE)

  # winsorisation
  if(!is.null(w)) {
    stopifnot(w < 1 & w >= 0)
    zz <- sort(zz)
    upper <- ceiling(nn*(1-w))
    lower <- floor(nn*(w))
    zz[zz > zz[upper]] <- zz[upper]
    zz[zz < zz[lower]] <- zz[lower]
    inflation_factor2 <- (1/nn)*sum(zz^2)
  }

  # method of moments RE approach
  pp <- pfunnel$grouped$prop_obs
  vv <- pp*(1 - pp)*(1/nn)
  ww <- 1/vv
  tau2 <- (nn*inflation_factor2 - (nn - 1)) / (sum(ww) - (sum(ww^2)/sum(ww)))

  # return
  c(stat = chi, pval = prob, inflation_factor = sqrt(inflation_factor2),
    tau2 = tau2)
}


#'
#'
#'
#'
predict_ <- function(mod, newdata) {
  linear <- newdata %*% mod$coefficients
  prob <- exp(linear) / (1 + exp(linear))
  as.numeric(prob)
}






#' Casemix adjustment
#'
#'
eval_indirect_adj <- function(pfunnel,method="cv",folds=10) {
  outcome <- pfunnel$outcome
  casemix_adj <- pfunnel$casemix_adj

  # move this from here....
  if(is.data.frame(casemix_adj)) {
    stopifnot(!is.na(casemix_adj))
    nms <- names(casemix_adj)
    form <- as.formula(paste0("~",paste0(nms,collapse="+")))
    mm <- model.matrix(form,data=casemix_adj)
  } else {
    mm <- casemix_adj
  }

  # note dependency to above
  nn <- nrow(mm)

  if (method == "cv") {
    # make folds
    per_samp <- floor(nn/folds)
    samp_i <- 1:nn
    avail <- rep(TRUE,length(samp_i))
    make_folds <- function(x) {
      out <- sample(x[avail],size = per_samp,replace=FALSE)
      avail[out] <<- FALSE
      out
    }
    cv_leave_out_i <- replicate(n = folds, make_folds(samp_i))

    # make folds
    cv_folds <- vector(mode = "list", length = folds)
    for(leave_out in 1:ncol(cv_leave_out_i)) {
      train_i <- as.numeric(cv_leave_out_i[,-leave_out])
      test_i <- as.numeric(cv_leave_out_i[,leave_out])
      cv_folds[[leave_out]] <- list(train = train_i, test = test_i)
    }

    # evaluate fold times (cols = folds)
    brier <- numeric(length(cv_folds))
    accuracy <- numeric(length(cv_folds))
    pseudoR2 <- numeric(length(cv_folds))
    auc_roc <- numeric(length(cv_folds))
    auc_pr <- numeric(length(cv_folds))

    for(fold in 1:length(cv_folds)) {
      # fit null model
      null_mod_i <- glm.fit(x = mm[cv_folds[[fold]]$train,1,drop=FALSE],
                            y = outcome[cv_folds[[fold]]$train],
                            family=binomial(link="logit"))
      class(null_mod_i) <- c(class(null_mod_i),"glm")
      # fit model
      adj_mod_i <- glm.fit(x = mm[cv_folds[[fold]]$train,],
                           y = outcome[cv_folds[[fold]]$train],
                           family=binomial(link="logit"))
      class(adj_mod_i) <- c(class(adj_mod_i),"glm")
      # predictions
      predict_out <- predict_(adj_mod_i, mm[cv_folds[[fold]]$test,])
      true_out <- outcome[cv_folds[[fold]]$test]
      binary_pred <- classify(predict_out,cutoff = 0.5)
      # evaluation metrics
      brier[fold] <- mean((predict_out - true_out)^2)
      accuracy[fold] <- accuracy(binary_pred,true_out)
      pseudoR2[fold] <- 1-(logLik(adj_mod_i)/logLik(null_mod_i))
      auc_roc[fold] <- auc_roc(predict_out,true_out)
      auc_pr[fold] <- auc_pr(predict_out,true_out)
    }
  }

  # no information rate
  no_info_rate <- max(mean(true_out),1-mean(true_out))

  # no information brier
  no_info_brier <- mean((mean(true_out) - true_out)^2)

  # outcomes
  list(no_info_brier = no_info_brier, brier = mean(brier),
       no_info_rate = no_info_rate,
       accuracy = mean(accuracy), pseudoR2 = mean(pseudoR2),
       auc_roc = mean(auc_roc), auc_pr = mean(auc_pr))
}


#'
#'
#'
#'
#'
pointTarget <- function(limits = 0.05, normalApprox = TRUE, crtlOverDisp = FALSE,
                        w = 0, multAdj = "none") {
  list(limits = limits, normalApprox = normalApprox, crtlOverDisp = crtlOverDisp,
       multAdj = multAdj, w = w)
}

#'
#'
#'
#'
distTarget <- function(limits = 0.05, w = 0) {
  list(limits = limits, w = w)
}


#' Create the funnel plot
#'
#' @param mult_adj multiple testing adjustment
#'
#' CHANGE normal_approx TO METHOD = "exact", "normal" AND SO ON!
#'
funnel <- function(y,group,x=NULL,target=pointTarget()) {
  # checks
  stopifnot(is.factor(group))
  stopifnot(length(unique(group)) != length(group))
  stopifnot(is.numeric(outcome))
  stopifnot(is.data.frame(casemix_adj)|is.null(casemix_adj)|is.matrix(casemix_adj))
  stopifnot(all(limits < 1))
  stopifnot(all(limits > 0))
  stopifnot(method %in% c("exact","normal"))
  #stopifnot(ctrl_over_disp %in% c("none","multiplicative","additive"))

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
    adj_mod <- NULL
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

  # limits
  upper <- vector(mode = "list",length = length(limits))
  lower <- vector(mode = "list",length = length(limits))

  for(ii in 1:length(limits)) {
    if(method == "normal") {
      # distribution quantiles
      zz_lower <- qnorm(limits[ii]/2)
      zz_upper <- qnorm(1-(limits[ii]/2))
      # control limits
      if (random_effects == TRUE) {
        tau2 <- od_res["tau2"]
        upper[[ii]] <- target + zz_upper*sqrt(grouped$var + tau2)
        lower[[ii]] <- target + zz_lower*sqrt(grouped$var + tau2)
      } else {
        upper[[ii]] <- target + zz_upper*inflation_factor*sqrt(grouped$var)
        lower[[ii]] <- target + zz_lower*inflation_factor*sqrt(grouped$var)
      }
    } else if(method == "exact") {
      # continuity adjusted control limits
      upper[[ii]] <- cont_adjust_bin(1-(limits[ii]/2),pfunnel$grouped$n,target)
      lower[[ii]] <- cont_adjust_bin((limits[ii]/2),grouped$n,target)
      # over-dispersion adjustment
      upper[[ii]] <- pfunnel$target + inflation_factor*(upper[[ii]] - pfunnel$target)
      lower[[ii]] <- pfunnel$target + inflation_factor*(lower[[ii]] - pfunnel$target)
    }
  }
  tmp <- as.character(100* (1-(limits)))
  upper <- data.frame(upper)
  names(upper) <- paste0("upper_",tmp)
  lower <- data.frame(lower)
  names(lower) <- paste0("lower_",tmp)

  # determine clinics inside/outside limits
  incontrol <- vector(mode = "list",length=length(limits))
  for(ii in 1:length(limits)){
    incontrol[[ii]] <- (pfunnel$grouped$prop_obs <= upper[[ii]] &
                    pfunnel$grouped$prop_obs >= lower[[ii]])
  }
  incontrol <- data.frame(incontrol)
  names(incontrol) <- paste0("inside_",tmp)

  # output
  grouped <- data.frame(grouped,
                        upper,
                        lower,
                        incontrol,
                        row.names=NULL)
  control_limits = list(method=method, limits=limits)

  # add class
  out <- list(grouped = grouped, target = target,
              over_dispersion = list(control = ctrl_over_disp,
                                     stats = od_res,
                                     winsorisation = w,
                                     random_effects = random_effects),
              control_limits = control_limits,
              adj_model = adj_mod,
              casemix_adj = casemix_adj,
              outcome = outcome)
  class(out) <- c(class(out),"pfunnel")

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



