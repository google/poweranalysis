# Copyright 2022 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Log-likelihood Function
#'
#' Given a variable \code{y} in the parent environment (not passed into this
#' function), return an R function that computes the log-likelihood of \code{y}
#' for one of the following parametric distributions, specified in the argument
#' \code{dist}:
#' (1) One of three symmetric distributions: the normal (\code{"norm"}),
#'     logistic (\code{"logis"}), and Cauchy (\code{"cauchy"}) distributions, or
#' (2) One of two right-skewed distributions: the lognormal (\code{"lnorm"}) and
#'     gamma (\code{"gamma"}) distributions.
#'
#' @param dist One of the distribution keyword strings specified above.
#' @return A function taking a vector \code{par} of parameters and a vector
#'   \code{y} of observed data to compute the log-likelihood of the observed
#'   data in \code{y} for the parametric distribution specified in \code{dist}.
#'
#' @export

LogLikFun <- function(dist = c("norm", "logis", "cauchy", "lnorm", "gamma")) {
  dist <- match.arg(dist)
  DensityFun <- match.fun(paste0("d", dist))
  out <- function(par, y) {
    return(sum(DensityFun(y, par[1], par[2], log = TRUE)))
  }
  return(out)
}

#' Start Values of Log-likelihood Parameters
#'
#' Given a numeric vector \code{y}, estimate the start values of the parameters
#' for maximizing the log-likelihood of \code{y} for one of the following
#' parametric distributions specified for the \code{LogLikFun} function:
#' (1) Normal (\code{"norm"}) and logistic (\code{"logis"}): The sample mean and
#'     standard deviation of \code{y}.
#' (2) Cauchy (\code{"cauchy"}): The sample median and interquartile range (IQR)
#'     of \code{y}.
#' (3) Lognormal (\code{"lnorm"}): The sample mean and standard deviation of
#'     \code{log(y)}.
#' (4) Gamma (\code{"gamma"}): The method-of-moments (MOM) estimators of the
#'     shape and rate parameters of \code{y}, based on the sample mean and
#'     standard deviation of \code{y}.
#'
#' For the normal and lognormal distributions, the start value of the first
#' parameter, the sample mean, is exactly the maximum likelihood estimate (MLE),
#' while the start value of the second parameter, the sample standard deviation,
#' is exactly the MLE if \code{sd_mle = TRUE} and asymptotically the MLE if
#' \code{sd_mle = FALSE}. Hence, not much work is required in estimating the
#' MLE for these two distributions, since the MLEs has nice analytic closed
#' forms, which can be used as start values.
#'
#' @param y A numeric vector containing the observed data for estimating the
#'   MLE of the parameters.
#' @param dist One of the distribution keyword strings specified above.
#' @param sd_mle Whether to use the MLE of the standard deviation and variance
#'   for the sample standard deviation in the start values. If \code{TRUE},
#'   the typical sample standard deviation and variance estimators are
#'   multiplied by, respectively, the factors \code{sqrt((n - 1) / n} and
#'   \code{(n - 1) / n}. If \code{FALSE}, the typical standard deviation and
#'   variance estimators are used, where the latter is an unbiased estimator of
#'   of the population variance.
#' @return A vector containing the start values for estimating the MLE of the
#'   parameters under the assumed parametric model.
#' @export

StartPars <- function(y, dist = c(
                        "norm", "logis", "cauchy", "lnorm", "gamma"
                      ), sd_mle = TRUE) {
  dist <- match.arg(dist)
  if (sd_mle) {
    n <- length(y)
    sd_factor <- sqrt((n - 1) / n)
  } else {
    sd_factor <- 1
  }
  out <- switch(dist,
    norm = c(mean(y), sd(y) * sd_factor),
    logis = c(mean(y), sd(y) * sd_factor),
    cauchy = c(median(y), IQR(y)),
    lnorm = c(mean(log(y)), sd(log(y)) * sd_factor),
    gamma = c(
      mean(y)^2 / (var(y) * sd_factor^2), var(y) * sd_factor^2 / mean(y)
    )
  )
  return(out)
}

#' Maximum Likehood Estimates of Model Parameters
#'
#' Given a numeric vector \code{y}, compute the maximum likelihood estimates
#' (MLEs) for one of the parametric distributions specified for the
#' \code{LogLikFun} and \code{StartPars} functions. See documentation for these
#' two functions for details.
#'
#' @param y A numeric vector containing the observed data for estimating the
#'   MLE of the parameters.
#' @param dist One of the distribution keyword strings specified in the
#'   \code{LogLikFun} function.
#' @param sd_mle Whether to use the MLE for the sample standard deviation and
#'   variance estimators in the start values. See documentation in the
#'   \code{StartPars} function for details.
#' @param ... Additional arguments to pass into the \code{base::optim} function,
#'   e.g., the \code{method} argument for the optimization method. See
#'   documentationn for the \code{base::optim} for details.
#' @return A vector containing the MLEs of the parameters under the assumed
#'   parametric model.
#' @export

Mle <- function(y, dist = c("norm", "lnorm", "logis", "cauchy", "gamma"),
                sd_mle = TRUE, ...) {
  dist <- match.arg(dist)
  start_pars <- StartPars(y, dist = dist, sd_mle = sd_mle)
  log_lik_fun <- LogLikFun(dist = dist)
  # Note that the y argument is the additional argument of the log_lik_fun
  # function, which takes par as its first argument.
  out <- optim(
    par = start_pars, fn = log_lik_fun, control = list(fnscale = -1), y = y, ...
  )$par
  return(out)
}

#' Sparsity Function Evaluated at Quantile
#'
#' Given a numeric vector \code{y} of observations and \code{tau} for the
#' quantile probability, estimate the sparsity function, evaluated at the
#' \code{tau}-th quantile.
#'   The sparsity function, a term coined by Tukey (1965), refers the inverse of
#' the density function. As the density function is the slope of the
#' distribution function, the sparsity function, as its inverse, is the slope
#' of the quantile function, the inverse of the distribution function. That is,
#' given a density function \code{f} of distribution function \code{F}, the
#' sparsity function is simply \code{1 / f}, the slope of the quantile function
#' \code{Q}, the inverse of \code{F}.
#'   An important fact is that if the sample observations are iid, the
#' asymptotic variance of a sample quantile, a difference between two sample
#' quantiles (of the same probability), or any regression coefficient of a
#' quantile regression is proportional to the sparsity function evaluated at the
#' quantile of interest. In particular, it is proportional to the variance term
#' \code{tau * (1 - tau) * (1 / f(Q(tau)))}, where \code{1 / f(Q(tau))} is the
#' sparsity evaluated at the quantile of interest.
#'   The result for the asymptotic variance of a sample quantile for iid
#' has been known for some time and can be found in, e.g., Cramer (1946). The
#' result for any regression coefficient of a quantile regression is first shown
#' in Koenker and Bassett (1978). Hence, estimating the asymptotic variance of a
#' sample quantile, a sample quantile difference, or a coefficient estimator in
#' a quantile regression depends in large part on the estimation of the sparsity
#' \code{1 / f(Q(tau))}.
#'   This function offers three approaches for estimating the sparsity:
#' (1) A nonparametric approach using the Siddiqui-Bloch-Gastwirth (SBG) method
#'     of Siddiqui (1960) and Block and Gaswirth (1968), implemented for the
#'     \code{se = "iid"} option of the \code{quantreg::summary.rq} function.
#' (2) An alternative nonparametric approach based on the inferring the
#'     sparsity from the width of a nonparametric quantile confidence interval
#'     (CI), based on the CI's asymptotic properties. This idea was first
#'     proposed by Hao et al. (2019), and to our knowledge, there isn't any
#'     mention of this in existing literature. We will refer to this as the CI
#'     method.
#' (3) A parametric approach, using one of the parametric distributions defined
#'     in the \code{LogLikFun} function. In this approach
#'   We recommend using the nonparametric approach, when the underlying
#' distribution of the data is unknown, which is typically the case. The
#' parametric approach can be used if we know the underlying distribution or
#' as validation for the nonparametric approach. For the nonparametric
#' approach, we recommend using the SBG method over the CI method, especially
#' for small and moderate sample sizes, since the former has better properties
#' based on simulation exercises and is the default method used for the
#' \code{se = "iid"} option of the \code{quantreg::summary.rq} function for
#' sparsity estimation. However, when the sample size is large, e.g., >= 100K
#' observations, the SBG method may be computationally intensive, due to the
#' need to fit a quantile regression on large data. Since the difference
#' between the SBG and CI methods for large sample sizes may be small, one may
#' prefer the CI method over the SBG method due to computational speed.
#'
#' @param y Either numeric vector containing the observed response values or
#'   quantile regression residuals for estimating the sparsity, or
#'   a \code{quantreg::rq} object from which the quantile regression residuals
#'   will be extracted. Note that for the SBG method, the numeric vector
#'   \code{y} will need to be residualized first, i.e., converted into quantile
#'   regression residuals. This can be done using before passing the values into
#'   the argument \code{y} or using the \code{residualize = TRUE} option if the
#'   \code{y} values passed into the function are not quantile regression
#'   residuals.
#' @param tau The quantile probability for the estimated quantile.
#' @param method Either \code{"npar_sbg"} for the nonparmaetric SBG method, or
#'   \code{"npar_ci"} for the nonparametric CI method, or \code{"par"} for the
#'   parametric method.
#' @param ci_method The method for obtaining the nonparametric quantile CI,
#'   equal to the value of the \code{method} argument of the
#'   \code{QuantileCI} function. See the documentation for the
#'   \code{QuantileCI} function for details.
#' @param dist The parametric distribution for the parametric method, equal to
#'   the value of the \code{dist} argument for the \code{LogLikFun} function.
#'   See the documentation for the \code{LogLikFun} function for details.
#' @param residualize Whether to residualize \code{y} using a simple
#'   location-model quantile regression.
#' @param p The number of parameters estimated in the linear conditional
#'   quantile of a quantile regression, if \code{y} is a vector of quantile
#'   regression residuals.
#' @param hs For the SBG method (\code{method = "npar_sbg"}), whether to use the
#'   method of Hall and Sheather (1988) to define the bandwidth to estimate
#'   sparsity. If \code{hs = FALSE}, the Bofinger (1975) method is used.
#' @param sig_level The significance level used to define (1) the z quantiles
#'   in the bandwidth selected using either the Hall and Sheather (1988) method
#'   or the Bofinger (1975) method, for the SBG method, and (2) the confidence
#'   level of the interval for the CI method.
#' @param residualize Whether to residualize \code{y} using a simple
#'   location-model quantile regression, if \code{y} is a numeric vector.
#'   This is needed only for the SBG method, but we provide the option for the
#'   other methods for benchmarking purposes.
#' @return A scalar estimate of the sparsity evaluated at the quantile of
#'   interest.
#' @references
#' \itemize{
#' \item Bloch, Daniel A., and Gastwirth, Joseph L. (1968). "On a Simple
#'   Estimate of the Reciprocal of the Density Function." \emph{Annals of
#'   Mathematical Statistics,} 39: 1083-1085.
#' \item Bofinger, Eve (1975). "Estimation of a Density Function Using Order
#'   Statistics." \emph{Australian and New Zealand Journal of Statistiscs}, 17:
#'   1-7.
#' \item Cramer, Harald (1946). \emph{Mathematical Methods of Statistics}.
#'   Princeton University Press.
#' \item Hall, Peter and Simon J. Sheather (1988). "On the Distribution of a
#'   Studentized Quantile." \emph{Journal of the Royal Statistical Society
#'   Series B}, 50: 381-391.
#' \item Hao, Siteng, Tianhong He, and Meeyoung Park (2019). "Statistical
#'   Methods of Percentile Estimation and Inference in the Page Load Latency
#'   Study." Google Internal Document.
#' \item Koenker, Roger and Gilbert Bassett, Jr. (1978). "Regression Quantiles."
#'   \emph{Econometrica}, 46(1): 33-50.
#' \item Siddiqui, M. M. (1960). "Distribution of Quantiles in Samples from a
#'   Bivariate Population." \emph{Journal of Research of the National Bureau of
#'   Standards,} 64B(3): 145-150.
#' \item Tukey, John Wilder (1965). "Which Part of the Sample Contains the
#'   Information?" \emph{Proceedings of the National Academy of Sciences of the
#'   United States of America}, 53(1): 127-134.
#' }
#' @export

Sparsity <- function(y, tau = 0.5, method = c("npar_sbg", "npar_ci", "par"),
                     ci_method = c("interpolate", "exact", "normal_approx"),
                     dist = c("norm", "logis", "cauchy", "lnorm", "gamma"),
                     p = 1, hs = TRUE, sig_level = 0.05, pars = NULL,
                     residualize = FALSE, tol = .Machine$double.eps^0.5) {
  method <- match.arg(method)
  # Extract residuals of a quantile regression under two scenarios:
  # (1) If y is a quantreg::rq object, set p to the number of coefficients
  #     estimated and set y to equal the residuals of the quantile regression.
  # (2) If y is a numeric vector, and residualize == TRUE, perform a simple
  #     location-model quantile regression and set y equal to the residuals and
  #     p to 1, for the intercept coefficient.
  if (inherits(y, "rq")) {
    p <- length(coef(y))
    y <- resid(y)
  } else if (residualize) {
    location_model <- quantreg::rq(y ~ 1, tau = tau)
    y <- y - quantile(y, probs = tau)
    p <- 1
  }
  n <- length(y)
  if (method == "npar_sbg") {
    # This section is taken almost verbatim from the code to estimate sparsity
    # for the \code{se = "iid"} option of the \code{quantreg::summary.rq}
    # function, i.e., mainly from this section of the function code:
    # http://google3/third_party/R/packages/quantreg/quantreg/R/quantreg.R;l=1093-1098;rcl=217341861
    pz <- sum(abs(y) < tol)
    bw <- quantreg::bandwidth.rq(p = tau, n = n, hs = hs, alpha = sig_level)
    h <- max(p + 1, ceiling(n * bw))
    ir <- (pz + 1):(h + pz + 1)
    ord_y <- sort(y[order(abs(y))][ir])
    xt <- ir / (n - p)
    out <- as.numeric(quantreg::rq(ord_y ~ xt)$coef[2])
  } else if (method == "npar_ci") {
    ci_method <- match.arg(ci_method)
    out <- QuantileCI(
      y,
      tau = tau, alternative = "two_sided", method = ci_method,
      conf_level = 1 - sig_level
    ) %>%
      TwoSidedCI2SE(conf_level = 1 - sig_level) %>%
      SE2SD(n = n) %>%
      SD2Sparsity(tau = tau)
  } else if (method == "par") {
    dist <- match.arg(dist)
    if (any(y < 0) && dist %in% c("lnorm", "gamma")) {
      stop_msg <- glue("Cannot fit {dist} distribution on negative values.")
      stop(stop_msg)
    }
    if (is.null(pars)) {
      if (dist %in% c("norm", "lnorm")) {
        pars <- StartPars(y, dist = dist, sd_mle = FALSE)
      } else {
        pars <- Mle(y, dist = dist)
      }
    }
    text_exp <- glue("
      1 / (
        tau %>%
          q{dist}(pars[1], pars[2]) %>%
          d{dist}(pars[1], pars[2])
      )
    ")
    out <- eval(parse(text = text_exp))
  }
  return(out)
}

#' Standard Deviation and Standard Error for a Sample Quantile
#'
#' Given a numeric vector \code{y} of observations and \code{tau} for the
#' quantile probability, estimate the following asymptotic measures:
#' (1) The asymptotic population standard deviation (SD) for the sample
#'     \code{tau}-th quantile, given by
#'     \code{tau * (1 - tau) * (1 / f(Q(tau)))}, where \code{1 / f(Q(tau))} is
#'     the sparsity evaluated at the quantile of interest.
#' (2) The asymptotic standard error (SE) for the sample quantile, defined as
#'     the asymptotic population SD divided by \code{sqrt(n)}, i.e., the
#'     square root of the sample size.
#'   Note that the asymptotic population SD of the sample quantile is denoted by
#' \code{omega} in Equation 3.2.1 in Koenker (2005). See documentation for the
#' \code{Sparsity} function for details.
#'
#' @param y A numeric vector containing the observed response values or
#'   quantile regression residuals for estimating the asymptotic SD or SE.
#' @param ... Additional arguments to pass into the \code{Sparsity} function for
#'   the \code{QuantileSD} function, or into the \code{QuantileSD} function for
#'   the \code{QuantileSE} function.
#' @return A scalar estimate of the asymptotic population SD of the sample
#'   quantile for the \code{QuantileSD} function and the asymptotic SE for the
#'   \code{QuantileSE} function.
#' @references
#' \itemize{
#' \item Koenker, Roger (2005). \emph{Quantile Regression}. Cambridge University
#'   Press.
#' }
#' @export

QuantileSD <- function(y, tau = 0.5, ...) {
  sparsity <- Sparsity(y = y, tau = tau, ...)
  out <- sparsity * sqrt(tau * (1 - tau))
  return(out)
}
QuantileSE <- function(y, ...) {
  out <- QuantileSD(y = y, ...) / sqrt(length(y))
  return(out)
}

#' Confidence Interval for a Sample Quantile
#'
#' Estimate the confidence interval (CI) of a sample quantile using order
#' statistics, with the following three variants:
#' (1) Asymptotic normal approximation using uninterpolated order statistics.
#' (2) Exact method using uninterpolated order statistics.
#' (3) Exact method using interpolated order statistics, first proposed by
#'     Hettmansperger and Sheather (1986) for the median and later extended to
#'     arbitrary quantiles by Nyblom (1992).
#'   This is currently implemented as a wrapper function around the following
#' three functions in the \code{EnvStats} that computes the sample quantile CIs
#' for, respectivel, the above three methods:
#' (1) \code{EnvStats::ci.qnpar.normal.approx} for the normal approximation
#'     method.
#' (2) \code{EnvStats::ci.qnpar.exact} for the exact method without
#'     interpolation.
#' (3) \code{EnvStats::ci.qnpar.interpolate} for the exact method with
#'     interpolation.
#'
#' @param y A numeric vector containing observed data to estimate the sample
#'   quantile CI for.
#' @param tau The quantile probability for the sample quantile.
#' @param method One of the three methods for obtaining the nonparametric
#'   quantile CI, i.e., \code{"interpolate"} for the exact method with
#'   interpolation, \code{"exact"} for the exact method without interpolation,
#'   and \code{"normal_approx"} for using an asymptotic normal approximation for
#'   the binomial CI.
#' @param alternative Whether the CI corresponds to one- or two-sided hypothesis
#'   test. Defaults to \code{"two_sided"} for a two-sided test.
#' @param conf_level The confidence level for the CI. This is approximate for
#'   the exact method.
#' @return A two-entry numeric vector containing lower and upper bounds for the
#'   sample quantile CI estimate.
#' @references
#' \itemize{
#' \item Hettmansperger, Thomas P., and Sheather, Simon J. (1986). "Confidence
#'   Intervals Based on Interpolated Order Statistics." \emph{Statistics &
#'   Probability Letters,} 4:75-79.
#' \item Nyblom, Jukka (1992). "Note on Interpolated Order Statistics."
#'   \emph{Statistics & Probability Letters,} 14: 129-131.
#' }
#' @export

QuantileCI <- function(y, tau = 0.5, method = c(
                         "interpolate", "exact", "normal_approx"
                       ), alternative = c("two_sided", "one_sided"),
                       conf_level = 0.95, min_coverage = TRUE) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  # The following code is modified from a section in the
  # \code{EnvStats::ci.qnpar} function:
  # http://google3/third_party/R/packages/EnvStats/EnvStats/R/ci.qnpar.R;l=53-60;rcl=399163923
  arg_list <- list(
    x = sort(y), n = length(y), p = tau, lb = -Inf, ub = Inf,
    ci.type = gsub("_", "-", alternative, fixed = TRUE),
    approx.conf.level = conf_level
  )
  if (method == "exact") {
    arg_list <- c(arg_list, list(min.coverage = min_coverage, tol = 0))
  }
  out <- do.call(
    paste("ci.qnpar", gsub("_", ".", method, fixed = TRUE), sep = "."),
    arg_list
  )$ci.limits %>%
    as.numeric()
  return(out)
}

#' Converting a Two-sided Confidence Interval to a Sparsity Estimate
#'
#' Three intermediate functions to convert a two-sided symmetric normal
#' confidence interval (CI) to a sparsity estimate:
#' (1) Converting from the two-sided CI to the standard error (SE)
#' (2) Converting from an SE to a standard deviation (SD)
#' (3) Converting from an asymptotic quantile SD to the sparsity evaluated at
#'     the quantile of interest.
#'
#' @param ci A two-entry numeric vector containing the limits of two-sided CI.
#' @param conf_level The confidence level of the two-sided CI.
#' @param se The value of the SE.
#' @param n The number of observations.
#' @param sd The value of the SD.
#' @param tau The probability for the quantile of interest.
#' @return The SE for the \code{TwoSidedCI2SE} function, the SD for the
#'    \code{SE2SD} function, and the sparsity evaluated at the quantile of
#'    interest for \code{SD2Sparsity} function.
#' @export

TwoSidedCI2SE <- function(ci = qnorm(0.025) * c(-1, 1),
                          conf_level = 0.95) {
  ci_half_width <- (ci[2] - ci[1]) / 2
  out <- ci_half_width / qnorm((1 + conf_level) / 2)
  return(out)
}
SE2SD <- function(se = 1, n = 100) {
  out <- se * sqrt(n)
  return(out)
}
SD2Sparsity <- function(sd = 1, tau = 0.50) {
  out <- sd / sqrt(tau * (1 - tau))
  return(out)
}

#' Intra-cluster Correlation
#'
#' Compute the intra-cluster correlation (ICC) for a simple random intercept
#' model without covariates.
#'
#' @param data A data frame containing the response variable in \code{y} and
#'   and the cluster id variable in \code{cluster_id}.
#' @param y The string name of the response variable.
#' @param cluster_id The string name of the cluster id variable.
#' @return The ICC of the response variable \code{y} given the clustering
#'   specified by \code{cluster_id}.
#' @export

ICC <- function(data, y, cluster_id = NULL) {
  if (is.null(cluster_id)) {
    out <- 0
  } else {
    model_formula <- glue("{y} ~ (1 | {cluster_id})") %>%
      as.formula()
    varcor <- summary(lme4::lmer(model_formula, data = data))$varcor
    out <- as.numeric(varcor[[1]] / attr(varcor, "sc")^2)
  }
  return(out)
}

#' Variance-stabilizing Transformations
#'
#' Apply the appropriate variance-stabilizing transformation (VST) or its
#' inverse for specific parametric distributions, implemented for the estimation
#' of the following parameters:
#' (1) Poisson mean: The square-root transformation.
#' (2) Binomial proportion: The arcsine transformation, also known as the
#'     arcsine-square-root transformation or the angular transformation.
#'
#' @param y A numeric vector containing data to apply the VST to.
#' @param dist The distribution of the underlying data.
#' @param inverse Whether to apply the inverse of the VST instead of the VST.
#' @param const The proportionality constant for the VST.
#' @return The VST or the inverse VST of the data in \code{x}.
#' @export

VarStabTransf <- function(y, dist = c("poisson", "binomial"), const = NULL,
                          inverse = FALSE) {
  if (is.null(const)) {
    const <- 2
  }
  if (inverse) {
    y <- y / const
    out <- switch(dist,
      poisson = y^2,
      binomial = sin(y)^2
    )
  } else {
    out <- const * switch(dist,
      poisson = sqrt(y),
      binomial = asin(sqrt(y))
    )
  }
  return(out)
}
PoissonVST <- function(y, ...) {
  return(VarStabTransf(y = y, dist = "poisson", ...))
}
BinomialVST <- function(y, ...) {
  return(VarStabTransf(y = y, dist = "binomial", ...))
}
ArcsineTransf <- function(y, ...) {
  return(VarStabTransf(y = y, dist = "binomial", ...))
}
InvArcsineTransf <- function(y, ...) {
  return(VarStabTransf(y = y, dist = "binomial", inverse = TRUE, ...))
}
