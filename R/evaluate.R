#########################################################################
# Funnel plot functions
# - How good is the case mix adjustment?
#########################################################################


#' Evaluate ability of casemix adjustment variables to predict outcome
#'
#' @param funnelRes a funnel plot object
#' @param method "cv" for cross-validation
#' @param folds an integer (>= 2) indicating the number of cross validation folds
#'
#' @section Note:
#' Current implementation is for test purposes and may change substantially.
#'
#' @export
evalCasemixAdj <- function(funnelRes,method="cv",folds=5L) {

  # checks
  assertthat::assert_that(method == "cv",msg = "only cross validation currently supported")
  assertthat::assert_that(is.numeric(folds))
  assertthat::assert_that("funnelRes" %in% class(funnelRes))
  assertthat::assert_that(nrow(funnelRes$data)/folds > 2)

  ## edit formula
  sepFormula <- getFunnelFormula(funnelRes$formula)
  idVar <- sepFormula$id
  outcomeVar <- sepFormula$outcome
  newFormula <- sepFormula$newForm

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
    auc_roc <- numeric(length(cv_folds))

    for(fold in 1:length(cv_folds)) {
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
      auc_roc[fold] <- auc_roc(predict_out,true_out)
    }
  }

  # no information rate
  outcome <- funnelRes$data[[outcomeVar]]
  no_info_rate <- max(mean(outcome),1-mean(outcome))

  # return
  data.frame(run = c(1:folds,"overall"),
             brier = c(brier,mean(brier)),
             accuracy = c(accuracy,mean(accuracy)),
             auc_roc = c(auc_roc,mean(auc_roc)),
             no_info_rate = c(rep(NA,folds),no_info_rate))
}


