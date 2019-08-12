#' List institutional outliers
#'
#' @param funnelRes funnel plot object
#'
#' @export
outliers <- function(funnelRes) {
  ctrl_names <- grep("inside",names(funnelRes$results),value=TRUE)
  print_names <- sub("inside_","",ctrl_names)
  outlierIndex <- (rowSums(funnelRes$results[,ctrl_names,drop=FALSE]) == 0)
  funnelRes$results[outlierIndex,]
}

#' print function
#'
#' @param x funnel plot object.
#' @param ... other arguments to print.
#'
#' @export
print.funnelRes <- function(x,...) {
  print(format(x$results, digits=3))
}
