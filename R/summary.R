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
#' @param funnelRes funnel plot object
#'
#' @export
print.funnelRes <- function(funnelRes) {
  print(format(funnelRes$results, digits=3))
}
