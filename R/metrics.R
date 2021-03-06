#########################################################################
# Funnel plot functions
# - How reliable is a predicted outcome?
#########################################################################

#' classify a probability
#'
#' @param prob probability of event
#' @param cutoff threshold above which events are classified as 1, else 0.
#'
classify <- function(prob, cutoff = 0.5) {
  pred <- rep(0, length(prob))
  pred[prob >= cutoff] <- 1
  pred
}

#' labels data like a confusion matrix
#'
#' @param pred prediction (binary)
#' @param truth observed value
#'
confusion <- function(pred,truth) {
  labelled <- rep("TN",length(pred))
  labelled[(pred == truth) & (truth == 1)] <- "TP"
  labelled[(pred != truth) & (truth == 0)] <- "FP"
  labelled[(pred != truth) & (truth == 1)] <- "FN"
  labelled
}

#' accuracy
#'
#' @inheritParams confusion
#'
accuracy <- function(pred,truth) {
  sum(pred == truth)/length(pred)
}

#' error rate
#'
#' @inheritParams confusion
#'
err <- function(pred,truth) {
  1 - accuracy(pred,truth)
}

#' recall, sensitivity, true positive rate
#'
#' @inheritParams confusion
#'
tpr <- function(pred,truth) {
  tp <- sum((pred == truth) & (truth == 1))
  fn <- sum((pred != truth) & (truth == 1))
  tp / (tp + fn)
}

#' specificity
#'
#' @inheritParams confusion
#'
sp <- function(pred,truth) {
  tn <- sum((pred == truth) & (truth == 0))
  fp <- sum((pred != truth) & (truth == 0))
  tn / (tn + fp)
}

#' fpr
#'
#' @inheritParams confusion
#'
fpr <- function(pred,truth) {
  fp <- sum((pred != truth) & (truth == 0))
  tn <- sum((pred == truth) & (truth == 0))
  fp / (tn + fp)
}

#' precision, positive predictive value
#'
#' @inheritParams confusion
#'
ppv <- function(pred,truth) {
  tp <- sum((pred == truth) & (truth == 1))
  fp <- sum((pred != truth) & (truth == 0))
  tp / (tp + fp)
}

#' Receiver operating characteristic curve
#'
#' @param prob prediction (probability)
#' @param length_out 100
#' @inheritParams confusion
#'
roc <- function(prob,truth,length_out=100) {
  cutoffs <- seq(0,1,length.out=length_out)
  res <- data.frame(cutoffs, tpr = NA, fpr = NA)
  for (i in 1:length_out) {
    pred <- classify(prob,cutoff = cutoffs[i])
    res$tpr[i] <- tpr(pred,truth)
    res$fpr[i] <- fpr(pred,truth)
  }
  res
}

#' Area under the receiver operating characteristic curve
#'
#' @inheritParams roc
#'
auc_roc <- function(prob,truth,length_out=100) {
  roc_df <- roc(prob,truth,length_out=length_out)
  out <- numeric(length_out-1)
  for (i in 1:length(out)) {
    y <- rev(roc_df$tpr)
    x <- rev(roc_df$fpr)
    out[i] <- ((y[i] + y[i+1])/2)*(x[i+1]-x[i])
  }
  sum(out,na.rm=TRUE)
}


#' precision recall curves
#'
#' @inheritParams roc
#'
pr <- function(prob,truth,length_out=100) {
  cutoffs <- seq(0,1,length.out=length_out)
  res <- data.frame(cutoffs, ppv = NA, tpr = NA)
  for (i in 1:length_out) {
    pred <- classify(prob,cutoff = cutoffs[i])
    res$ppv[i] <- ppv(pred,truth)
    res$tpr[i] <- tpr(pred,truth)
  }
  res
}
