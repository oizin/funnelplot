
#'
#'
#' @export
sample_nobs <- function(N,sample_func) {
  if (is.null(sample_func)) {
    n_i <- round(runif(n = N, min = n_min, max = n_max), 0)

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

#' @export
sample_means <- function(N,mean,sd) {
  mean_i <- rnorm(N, mean = mean, sd = sd)
  if (any(mean_i <= 0)) {
    replace <- (mean_i <= 0)
    mean_i[replace] <- sample_means(N = sum(replace),mean,sd)
  }
  mean_i
}

#' simulate random effects
#'
#'
#' @export
simulate <- function(N, n_min, n_max, mean, sd, sample_func = NULL) {

  if (pnorm(0, mean, sd) > 0.0) {
    warning("unreliable results where negative support on random effect distribution exceeds 0")
  }

  n_i <- sample_nobs(N,sample_func=sample_func)
  n_min <- min(n_i)
  n_max <- max(n_i)

  mean_i <- sample_means(N, mean, sd)

  m_i <- rbinom(n = N, size = n_i, prob = mean_i)
  y_i <- m_i/n_i
  se_i <- sqrt(y_i*(1 - y_i)*(1/n_i))
  out <- list(data = data.frame(n = n_i, y = y_i, se = se_i),
              params = c(n_min = n_min, n_max = n_max, mean = mean, sd = sd),
              mean = mean_i)
  class(out) <- c(class(out), "fsim")
  out
}


#' @export
plot.fsim <- function(sim,effect="true") {

  # variance calculations
  rng <- floor(seq(sim$params["n_min"],  sim$params["n_max"], length.out = 100))
  target <- sim$params["mean"]
  se <- sqrt(target*(1-target)*(1/rng))

  # estimated quantities
  if (effect == "estimated") {
    # random effect calculations
    N <- nrow(sim$data)
    s0 <- sqrt(target*(1-target)*(1/sim$data$n))
    z <- (sim$data$y - target)/s0
    theta <- (1/N)*sum(z*z)
    w <- 1/(sim$data$se^2)
    w <- w[w != Inf]
    tau2 <- ((N*theta) - (N - 1)) / (sum(w) - (sum(w*w)/sum(w)))
    # for the plot
    v_comb <- tau2 + se*se
  } else {
    v_comb <- se*se + sim$params["sd"]*sim$params["sd"]
    v_noncomb <- se*se
  }

  # calculate limits
  upper_95 <- target + 1.96*sqrt(v_comb)
  lower_95 <- target - 1.96*sqrt(v_comb)
  upper_99 <- target + 2.58*sqrt(v_comb)
  lower_99 <- target - 2.58*sqrt(v_comb)

  # the plot
  par(mfrow = c(1,2))
  # accurate method
  minmax <- c(min(sim$data$y) - 0.1, max(sim$data$y) + 0.1)
  plot(sim$data[,c("n","y")], pch = 20, ylim = minmax, main = "target: distribution")
  lines(rng, upper_95, lty = 2, col = "red")
  lines(rng, lower_95, lty = 2, col = "red")
  lines(rng, upper_99, lty = 2, col = "blue")
  lines(rng, lower_99, lty = 2, col = "blue")
  abline(h = target, col = "gray")
  legend("topright",legend = c("95% limits", "99% limits"), col = c("red", "blue"),
         lty = 2, bty = "n")
  # inaccurate method
  minmax <- c(min(sim$data$y) - 0.1, max(sim$data$y) + 0.1)
  plot(sim$data[,c("n","y")], pch = 20, ylim = minmax, main = "target: point")
  lines(rng, target + 1.96*sqrt(v_noncomb), lty = 2, col = "red")
  lines(rng, target - 1.96*sqrt(v_noncomb), lty = 2, col = "red")
  lines(rng, target + 2.58*sqrt(v_noncomb), lty = 2, col = "blue")
  lines(rng, target - 2.58*sqrt(v_noncomb), lty = 2, col = "blue")
  abline(h = target, col = "gray")
  legend("topright",legend = c("95% limits", "99% limits"), col = c("red", "blue"),
         lty = 2, bty = "n")
  par(mfrow = c(1,1))
}


# plot(simulate(100, 50, 1000, 0.22,  sqrt(0.007711165 )),effect="estimated")

# f = function(n) rnbinom(n, size = 1.1404360,mu = 775.0000000)
# plot(tt <- simulate(100, 50, 1000, 0.22,  .15, sample_func = f),effect="estimated")

