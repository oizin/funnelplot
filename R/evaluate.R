#' Casemix adjustment
#'
#' @export
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
#' @export
predict_ <- function(mod, newdata) {
  linear <- newdata %*% mod$coefficients
  prob <- exp(linear) / (1 + exp(linear))
  as.numeric(prob)
}
