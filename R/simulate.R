#########################################################################
# Funnel plot functions
# - Simulate example data to improve understanding
#########################################################################

#' Sample a certain number of observations
#'
#' @inheritParams generate_example
#'
sample_nobs <- function(N,n_min,n_max,sample_func) {
  if (is.null(sample_func)) {
    n_i <- round(stats::runif(n = N, min = n_min, max = n_max), 0)

  } else {
    n_i <- sample_func(n = N)

    if(any(n_i < n_min)) {
      replace <- (n_i < n_min)
      n <- sum(replace)
      n_i[replace] <- sample_nobs(n,sample_func)
    }
  }
  n_i
}

#' Sample means
#'
#' @inheritParams generate_example
#'
sample_means <- function(N,mean,sd) {
  mean_i <- stats::rnorm(N, mean = mean, sd = sd)
  if (any(mean_i <= 0)) {
    replace <- (mean_i <= 0)
    mean_i[replace] <- sample_means(N = sum(replace),mean,sd)
  }
  mean_i
}

#' Simulate data for a funnel plot
#'
#' @param N number of clusters (institutions)
#' @param n_min limit on the number of observations within the smallest cluster (institution)
#' @param n_max limit on the  of observations within the largest cluster (institution)
#' @param mean mean of cluster (institution) mean distribution
#' @param sd standard deviation of cluster (institution) mean distribution
#' @param sample_func function to sample jth cluster (institution) size in a form other than Unif(n_min,n_max)
#'
#' @section Note:
#' Current implementation is for test pruposes and may change substantially.
#'
#' @export
generate_example <- function(N, n_min, n_max, mean, sd = 0, sample_func = NULL) {

  if (stats::pnorm(0, mean, sd) > 0.0) {
    warning("unreliable results where negative support on random effect distribution exceeds 0")
  }
  # sample hyperparameters (sample sizes, cluster means, #success, #failure)
  n_i <- sample_nobs(N,n_min,n_max,sample_func=sample_func)
  mean_i <- sample_means(N, mean, sd)
  m_i <- stats::rbinom(n = N, size = n_i, prob = mean_i)
  d_i <- n_i - m_i

  # observation level data
  dij <- unlist(lapply(d_i,function(x) rep(0,x)))
  mij <- unlist(lapply(m_i,function(x) rep(1,x)))
  out <- data.frame(id = c(rep(1:N,times=d_i),rep(1:N,times=m_i)),
             y = c(dij,mij))
  out
}


#' Example sampling function for simulate
#'
#' @param n number of samples to take
#'
#' @export
eg_sample_func <- function(n) stats::rnbinom(n, size = 1.1404360,mu = 775.0000000)
