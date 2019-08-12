
#' Sample a certain number of observations
#'
#' @inheritParams simulate
#'
#' @export
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
#' @inheritParams simulate
#'
#' @export
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
#' @export
simulate <- function(N, n_min, n_max, mean, sd, sample_func = NULL) {

  if (stats::pnorm(0, mean, sd) > 0.0) {
    warning("unreliable results where negative support on random effect distribution exceeds 0")
  }

  n_i <- sample_nobs(N,n_min,n_max,sample_func=sample_func)
  n_min <- min(n_i)
  n_max <- max(n_i)

  mean_i <- sample_means(N, mean, sd)

  m_i <- stats::rbinom(n = N, size = n_i, prob = mean_i)
  y_i <- m_i/n_i
  se_i <- sqrt(y_i*(1 - y_i)*(1/n_i))

  out <- list(data = data.frame(n = n_i, y = y_i, se = se_i),
              params = c(n_min = n_min, n_max = n_max, mean = mean, sd = sd),
              mean = mean_i)
  class(out) <- c(class(out), "fsim")
  out
}


#' Plot the simulated data
#'
#' @param x the simulation results
#' @param effect whether to plot the "true" or "estimated" cluster (institution) outcome value
#' @param method whether to assume a "point" or "distribution" target
#' @param ... other arguments to plot.fsim
#'
#' @export
plot.fsim <- function(x,effect="true",method="dist",...) {

  # variance calculations
  rng <- floor(seq(x$params["n_min"],  x$params["n_max"], length.out = 100))
  target <- x$params["mean"]
  se <- sqrt(target*(1-target)*(1/rng))

  # estimated quantities
  if (effect == "estimated") {
    # random effect calculations
    N <- nrow(x$data)
    s0 <- sqrt(target*(1-target)*(1/x$data$n))
    z <- (x$data$y - target)/s0
    theta <- (1/N)*sum(z*z)
    w <- 1/(x$data$se^2)
    w <- w[w != Inf]
    tau2 <- ((N*theta) - (N - 1)) / (sum(w) - (sum(w*w)/sum(w)))
    # for the plot
    v_comb <- tau2 + se*se
  } else {
    v_comb <- se*se + x$params["sd"]*x$params["sd"]
    v_noncomb <- se*se
  }


  if(method == "dist") {
    # calculate limits
    upper_95 <- target + 1.96*sqrt(v_comb)
    lower_95 <- target - 1.96*sqrt(v_comb)
    upper_99 <- target + 2.58*sqrt(v_comb)
    lower_99 <- target - 2.58*sqrt(v_comb)

    pts <- data.frame(n = x$data$n, y = x$data$y)
    lmts <- data.frame(n = rng, upper_95, lower_95, upper_99, lower_99)
    lmts <- tidyr::gather(lmts,"limit","y",-1)
    lmts$limit_id <- c(rep("95%", length(rng)*2), rep("99%", length(rng)*2))

    out <- ggplot2::ggplot(data=pts,ggplot2::aes_(x="n",y="y")) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = x$params["mean"]) +
      ggplot2::geom_line(data=lmts,
                         ggplot2::aes_(x="n",y="y",group="limit",col="limit_id"),linetype=2)

  } else if (method == "point") {
    upper_95 <- target + 1.96*sqrt(v_noncomb)
    lower_95 <- target - 1.96*sqrt(v_noncomb)
    upper_99 <- target + 2.58*sqrt(v_noncomb)
    lower_99 <- target - 2.58*sqrt(v_noncomb)

    pts <- data.frame(n = x$data$n, y = x$data$y)
    lmts <- data.frame(n = rng, upper_95, lower_95, upper_99, lower_99)
    lmts <- tidyr::gather(lmts,"limit","y",-1)
    lmts$limit_id <- c(rep("95%", length(rng)*2), rep("99%", length(rng)*2))

    out <- ggplot2::ggplot(data=pts,ggplot2::aes_(x="n",y="y")) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = x$params["mean"]) +
      ggplot2::geom_line(data=lmts,
                         ggplot2::aes_(x="n",y="y",group="limit",col="limit_id"),linetype=2)
  }
  out
}


#' Example sampling function for simulate
#'
#' @param n number of samples to take
#'
#' @export
eg_sample_func <- function(n) stats::rnbinom(n, size = 1.1404360,mu = 775.0000000)
