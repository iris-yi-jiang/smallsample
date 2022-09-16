#' Uniform Points
#'
#' @description This function generates points targeting a uniform
#'    distribution.
#'
#' @param n Number of points to generate.
#' @param l Lower bound.
#' @param u Upper bound.
#' @return N-vector containing the sampled points.
#' @export
#'
#' @examples
#' gen_unif(n = 10, l = 2, u = 5)
#'
gen_unif <- function(n, l = 0, u = 1) {
  if (n <= 0) {
    stop("number of points should be greater than 0")
  }
  if (u < l) {
    stop("upper bound has to be greater than lower bound")
  }
  x <- c(NA)
  length(x) <- n
  b <- seq(from = l, to = u, length = (n + 1))
  for (i in 1:n) {
    x[i] <- stats::runif(n = 1, min = b[i], max = b[i + 1])
  }
  return(sample(x))
}

#' Normal Points
#'
#' @description This function generates points targeting a normal
#'    distribution.

#' @param n Number of points to generate.
#' @param mu Mean.
#' @param sig Standard deviation.
#' @return N-vector containing the sampled points.
#' @export
#'
#' @examples
#' gen_norm(n = 20, mu = 2, sig = 3)
#'
gen_norm <- function(n, mu = 0, sig = 1) {
  if (n <= 0) {
    stop("number of points should be greater than 0")
  }
  if (sig <= 0) {
    stop("standard deviation cannot be negative")
  }
  u <- gen_unif(n, l = 0, u = 1)
  x <- stats::qnorm(u, mean = mu, sd = sig)
  return(x)
}

#' Chi-Squared Points
#'
#' @description This function generates points targeting a Chi-Squared
#'    distribution.
#'
#' @param n Number of points to generate.
#' @param df Degrees of freedom.
#' @return N-vector containing the sampled points.
#' @export
gen_chisq <- function(n, df) {
  if (n <= 0) {
    stop("number of points should be greater than 0")
  }
  if (df <= 0) {
    stop("degrees of freedom has to be non-negative but can be non-integer")
  }
  u <- gen_unif(n, l = 0, u = 1)
  x <- stats::qchisq(u, df = df)
  return(x)
}

#' Exponential Points
#'
#' #' @description This function generates points targeting an Exponential
#'    distribution.
#'
#' @param n Number of points to generate.
#' @param rate Rate parameter.
#' @return N-vector containing the sampled points.
#' @export
gen_exp <- function(n, rate = 1) {
  if (n <= 0) {
    stop("number of points should be greater than 0")
  }
  if (rate <= 0) {
    stop("rate parameter should be positive")
  }
  u <- gen_unif(n, l = 0, u = 1)
  x <- stats::qexp(u, rate = rate)
  return(x)
}

gen_gh <- function(n, a = 0, b = 1, g = 0, h = 0) {
  u <- gen_unif(n, l = 0, u = 1)
  z <- stats::qnorm(u, mean = 0, sd = 1)
  if (g != 0) {
    gz <- (exp(g * z) - 1) / (g * z)
  } else {
    gz <- 1
  }
  hz <- exp((h * z^2) / 2)
  x <- a + b * gz * hz * z
  return(x)
}

gen_gk <- function(n, a = 0, b = 1, g = 0, k = 0) {
  u <- gen_unif(n, l = 0, u = 1)
  z <- stats::qnorm(u, mean = 0, sd = 1)
  gz <- 1 + 0.8 * (1 - exp(-g * z)) / (1 + exp(-g * z))
  kz <- (1 + z^2)^k
  x <- a + b * gz * kz * z
  return(x)
}
