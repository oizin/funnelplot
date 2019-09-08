#' List institutional outliers
#'
#' Return a ?type vector of all outliers
#'
#' @param funnelRes funnel plot object
#' @param limits Optional argument indicating which subset of the control limits should be used for selecting outlying institutions.
#'
#' @export
outliers <- function(funnelRes,limits=NULL) {
  ctrl_names <- grep("inside",names(funnelRes$results),value=TRUE)
  print_names <- sub("inside_","",ctrl_names)
  nctrl <- length(funnelRes$control$limits)
  outlierIndex <- (rowSums(funnelRes$results[,ctrl_names,drop=FALSE]) < nctrl)
  funnelRes$results[outlierIndex,"id"]
}

#' print function
#'
#' @param x funnel plot object.
#' @param ... other arguments to print.
#'
#' @export
print.funnelRes <- function(x,...) {
  tmp <- summary(x)

  cat("Number of outliers:\n")
  cat(tmp$n,"\n\n")

  cat("Summary of outliers:\n")
  print(tmp$outliers,row.names = FALSE)
  cat("\n")

  cat("Model:\n")
  cat(tmp$model,"\n")

}

#' Cluster performance summary
#'
#' Summarise the results of comparing the observed performance of clusters to their expected performance
#'
#' @param object funnel plot object.
#' @param ... other arguments to print.
#'
#' @export
summary.funnelRes <- function(object,...) {
  insideVars <- grepl("inside",names(object$results))
  insideVars <- names(object$results)[insideVars]
  lmtVars <- grepl("upper|lower",names(object$results))
  lmtVars <- names(object$results)[lmtVars]
  out <- vector(mode = "list", length = 3)
  nctrl <- length(object$control$limits)

  # number of outliers
  tmp <- unlist(lapply(object$results[insideVars],sum))
  out[[1]] <- rep(nrow(object$results),nctrl) - tmp

  # identities of outliers
  outlierIndex <- (rowSums(object$results[,insideVars,drop=FALSE]) < nctrl)
  out[[2]] <- object$results[outlierIndex,
    c("id","n","observed","expected","prop_obs","prop_adj",lmtVars,insideVars)]
  row.names(out[[2]]) <- NULL

  # model
  f_vars <- all.vars(object$formula)
  p <- length(f_vars)
  if(p > 2) {
    out[[3]] <- f_vars[-c(1,p)]
  } else {
    out[[3]] <- NA
  }

  names(out) <- c("n","outliers","model")

  out
}

