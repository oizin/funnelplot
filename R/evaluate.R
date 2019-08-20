#' Casemix adjustment
#'
#' @param funnelRes funnel plot object
#' @param method cross-validation
#' @param folds 10
#'
#' @export
evalCasemixAdj <- function(funnelRes,method="cv",folds=10) {

  # checks
  assertthat::assert_that(method == "cv",msg = "only cross validation currently supported")
  assertthat::assert_that(is.numeric(folds))
  assertthat::assert_that("funnelRes" %in% class(funnelRes))

  ## edit formula
  sepFormula <- getFunnelFormula(funnelRes$formula)
  idVar <- sepFormula$id
  outcomeVar <- sepFormula$outcome
  newFormula <- sepFormula$newForm
  nullFormula <- sepFormula$nullForm

  # note dependency to above
  N <- nrow(funnelRes$data)

  if (method == "cv") {
    # make folds
    per_samp <- floor(N/folds)
    samp_i <- 1:N
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
      null_mod_i <- stats::glm(nullFormula,
                        data = funnelRes$data[cv_folds[[fold]]$train,],
                        family=stats::binomial(link="logit"))
      # fit model
      adj_mod_i <- stats::glm(newFormula,
                       data = funnelRes$data[cv_folds[[fold]]$train,],
                       family=stats::binomial(link="logit"))
      # predictions
      predict_out <- stats::predict(adj_mod_i, funnelRes$data[cv_folds[[fold]]$test,],type="response")
      true_out <- funnelRes$data[cv_folds[[fold]]$test,outcomeVar]
      binary_pred <- classify(predict_out,cutoff = 0.5)
      # evaluation metrics
      brier[fold] <- mean((predict_out - true_out)^2)
      accuracy[fold] <- accuracy(binary_pred,true_out)
      pseudoR2[fold] <- 1-(stats::logLik(adj_mod_i)/stats::logLik(null_mod_i))
      auc_roc[fold] <- auc_roc(predict_out,true_out)
      auc_pr[fold] <- auc_pr(predict_out,true_out)
    }
  }

  # no information rate
  outcome <- funnelRes$data[[outcomeVar]]
  no_info_rate <- max(mean(outcome),1-mean(outcome))

  # no information brier
  no_info_brier <- mean((mean(outcome) - outcome)^2)

  # outcomes
  list(no_info_brier = no_info_brier, brier = mean(brier),
       no_info_rate = no_info_rate,
       accuracy = mean(accuracy), pseudoR2 = mean(pseudoR2),
       auc_roc = mean(auc_roc), auc_pr = mean(auc_pr))
}


