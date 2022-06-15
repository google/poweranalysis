#' Power and Precision Analysis
#'
#' Conduct a power analysis or precision analysis to determine the minimum
#' required sample size (MRSS) given a minimum detectable effect (MDE), or
#' vice versa, i.e., determine the MDE given an MRSS. Note that in the context
#' of a precision analysis, the MDE should be interpreted as the maximum error
#' margin (MME) instead of a detectable effect.
#'   While the name of the function \code{PowerPrecisionAnalysis} seems
#' long-winded, we typically use one of the many subvariants of the function,
#' as follows:
#' (1) \code{PowerAnalysis} and \code{PrecisionAnalysis} for, respectively,
#'     power analysis and precision analysis.
#' (2) \code{MinRequiredSampleSize}, \code{MinSampleSize}, or \code{MRSS} for
#'     performing a power/precision analysis to determine the MRSS.
#' (3) \code{MinDetectableEffect} or \code{MDE} for performing a power/precision
#'     analysis to determine the MDE.
#'   Note the following convention for the MRSS values:
#' (1) For the one-sample case, the MRSS values represent the sample size, i.e.,
#'     the total number of observations from the single sample.
#' (2) For the two-sample case, the MRSS values are in terms of the total number
#'     of observations (from the two samples) divided by two. For a balanced
#'     allocation, this is simply the number of observations per arm, for the
#'     two arms. For an unbalanced allocation, we can think of this as the
#'     average number of observations per arm.
#'
#' @param stat0 A numeric vector of the sample statistics for various
#'   conversion metrics, either sample means when \code{stat_type == "mean"} or
#'   sample quantiles when \code{stat_type == "quantile"}. For a one-sample
#'   analysis, the vector entries are the assumed values of the underlying
#'   parameter, either mean or quantile, under the null hypothesis, instead of
#'   the sample statistic. For a two-sample analysis, these are values of the
#'   sample statistic as estimates of the underlying parameter of the Control
#'   group population.
#' @param sd0 A numeric vector of the estimates of population standard
#'   deviations (SDs) for various conversion metrics, used only when
#'   \code{resp_type == "continuous"}, where we assume that the population SD is
#'   independent from the mean when \code{stat_type == "mean"} and the
#'   population SD is independent from the population quantile when
#'   \code{stat_type == "quantile"}. This could hold for quantiles, because
#'   we assume that changes to the quantile are due location shifts to the
#'   distribution, so that density or sparsity evaluated at the quantile of
#'   interest stays the same. When \code{resp_type == "binary"}, we compute
#'   the population SD in the function instead, because it is dependent on the
#'   value of \code{effect}.
#' @param cor0 A numeric vector of multiple correlation R values for the
#'   regression of the corresponding conversion on covariates. These could be
#'   the correlations between the post-period and pre-period metrics for the
#'   conversion for the regression of the post-period conversion on the
#'   pre-period conversion. But they could be the square-root of the R-squared
#'   values for the regression of the conversion on any set of covariates.
#'   See documentation for \code{AncovaDeff} for details.
#' @param cor0_treat A numeric vector of multiple correlation R values for the
#'   regression of the treatment variable on covariates, applicable only in a
#'   two-sample analysis, typically when there are imbalances between the
#'   study groups. See documentation for \code{AncovaDeff} for details.
#' @param n A numeric vector of MRSS values to determine the corresponding MDEs
#'   in a power analysis or precision analysis.
#' @param effect A numeric vector of MDE values to determine the corresponding
#'   MRSS values in a power analysis, or a numeric vector of maximum
#'   margin-of-error (MME) values in a precision analysis.
#' @param stat The power analysis or precision analysis statistic to determine,
#'   either \code{"mrss"} for the MRSS or \code{"mde"} for the MDE.
#' @param analysis The type of analysis conducted, either \code{"power"} for
#'   power analysis or \code{"precision"} for precision analysis.
#' @param resp_type The type of response variable for the conversion metrics,
#'   either \code{"binary"} for binary metrics or \code{"continuous"} for
#'   continuous metrics.
#' @param effect_type The type of effect for the MDE, either \code{"absolute"}
#'   for an absolute effect, \code{"relative"} for a relative effect, and
#'   \code{"effect_size"} for an effect size.
#' @param method The method used for the power analysis and/or precision
#'   analysis, either \code{"unconditional"} for the unconditional method,
#'   \code{"conditional"} for the conditional method, and \code{"asin"} for
#'   arcsine transformation method. Note that the unconditional method is
#'   available for all use cases. The conditional method is valid only for
#'   performing power analysis on conversions with binary responses. We
#'   allow for determining the MRSS only with the conditional method. While we
#'   we can determine the MDE as well using the conditional method, this feature
#'   is currently not implemented yet but may be implemented in the future.
#'   The arcsine method is valid is valid only for conversions with binary
#'   responses. We currently don't offer this method for precision analysis
#'   or for determining the MDE. This may change in the future. For more
#'   details, see http://go/cloud-sample-size.
#' @param alternative Whether the underlying hypothesis test for the power
#'   precision analysis is one- or two-sided. Defaults to \code{"two_sided"} for
#'   a two-sided test.
#' @param num_samples Either \code{"two_sample"} for a two-sample analysis or
#'   \code{"one_sample"} for a one-sample analysis. Note that the two-sample
#'   analysis is more common than the one-sample analysis, but we include the
#'   latter for completeness' sake.
#' @param sig_level The significance level for the power or precision analysis.
#'   Defaults to \code{0.05}.
#' @param power The minimum power required for the power analysis. Defaults to
#'   \code{0.80}.
#' @param phi The proportion of observations in one of the arms, either
#'   Treatment or Control, for a two-sample analysis. Defaults to \code{0.50},
#'   i.e., a balanced allocation for both arms. This determines the design
#'   effect for arm imbalance, given by \code{1 / (4 * phi * (1 - phi))}. For
#'   details, see documentation for the \code{ImbalanceDeff} function.
#' @param cluster_stats A named list of statistics for evaluating the cluster
#'   design effect, with \code{mean} denoting the mean of cluster sizes,
#'   \code{sd} denoting the standard deviation (SD) of cluster sizes, and
#'   \code{icc}, the intra-cluster correlation (ICC), i.e., the within-cluster
#'   correlation. See documentation for \code{ClusterDeff} for details.
#' @param cluster_deff_type The type of cluster design effect, either
#'   \code{"min"} for the lower bound or \code{"max"} for the upper bound.
#'   See documentation for the \code{ClusterDeff} function for details.
#' @return Either a numeric vector of MRSS values or MDE/MME values.
#' @references
#' \itemize{
#' \item Chow, Shien-Chung, Jun Shao, and Hansheng Wang (2003). \emph{Sample
#'   Size Calculations in Clinical Research.} Marcel Dekker: New York.
#' \item "Cloud Sample Size" (http://go/cloud-sample-size)
#' \item "Minimum Detectable Effect" (http://go/minimum-detectable-effect)
#' }
#' @export

PowerPrecisionAnalysis <- function(stat0, sd0 = NULL, cor0 = 0, cor0_treat = 0,
                                   n = 1e4, effect = 0, estimand = c(
                                     "mrss", "mde"
                                   ), analysis = c("power", "precision"),
                                   resp_type = c("binary", "continuous"),
                                   stat_type = c("mean", "quantile"),
                                   effect_type = c(
                                     "relative", "absolute", "effect_size"
                                   ), method = c(
                                     "unconditional", "conditional", "asin"
                                   ), alternative = c("two_sided", "one_sided"),
                                   num_samples = c("two_sample", "one_sample"),
                                   sig_level = 0.05, power = 0.80, phi = 0.50,
                                   cluster_stats = list(
                                     mean = 1, sd = 0, icc = 0
                                   ), cluster_deff_type = c("min", "max")) {
  # Read and check arguments
  analysis <- match.arg(analysis)
  resp_type <- match.arg(resp_type)
  stat_type <- match.arg(stat_type)
  effect_type <- match.arg(effect_type)
  alternative <- match.arg(alternative)
  num_samples <- match.arg(num_samples)
  cluster_deff_type <- match.arg(cluster_deff_type)
  if (resp_type == "continuous" ||
    (analysis == "precision" && method == "conditional")) {
    method <- "unconditional"
  }
  # Compute multipliers
  # (1) The sample-number multiplier based on the number of samples or groups.
  sample_num_multiplier <- ifelse(num_samples == "two_sample", 2, 1)
  # (2) The z-quantile multiplier using upper z quantiles
  #     based on significance level and power
  z_quantile_multiplier <- ZQuantileMultiplier(
    analysis = analysis, alternative = alternative, sig_level = sig_level,
    power = power
  )
  # Compute population standard deviation
  if (effect_type != "effect_size") {
    if (resp_type == "binary") {
      # Use mean0, since the statistic is always a mean for binary metrics
      mean0 <- stat0
      # Compute mean under alternative hypothesis and pooled mean
      mean1 <- switch(effect_type,
        absolute = mean0 + effect,
        relative = mean0 * (1 + effect),
      )
      mean_pooled <- 0.5 * (mean0 + mean1)
      # Compute population standard deviation
      sd_unconditional <- sqrt(
        0.5 * (mean0 * (1 - mean0) + mean1 * (1 - mean1))
      )
      if (method == "unconditional") {
        sd <- sd_unconditional
      } else if (method == "conditional") {
        z_upper_sig <- ZUpperSig(
          alternative = alternative, sig_level = sig_level
        )
        z_beta <- qnorm(power)
        w <- z_upper_sig / (z_upper_sig + z_beta)
        sd_pooled <- sqrt(mean_pooled * (1 - mean_pooled))
        sd <- w * sd_pooled + (1 - w) * sd_unconditional
      } else if (method == "asin") {
        sd <- 1
      }
    } else if (resp_type == "continuous") {
      sd <- sd0
    }
  }
  # Compute design effects for arm imbalance, ANCOVA, and clustering. Note that:
  # (1) The imbalance design effect is 1 for a one-sample analysis, since it
  #     only applies to a two-sample analysis.
  # (2) The ANCOVA and cluster design effects will be computed only if the test
  #     is a mean, not a quantile, and if the respective correlations are
  #     nonzero, i.e., the multiple correlation for ANCOVA design effect and the
  #     ICC for the cluster design effect.
  imbalance_deff <- ifelse(
    num_samples == "two_sample", ImbalanceDeff(phi = phi), 1
  )
  ancova_deff <- ifelse(
    stat_type == "mean" && (cor0 != 0 || cor0_treat != 0),
    AncovaDeff(r = cor0, r_treat = cor0_treat), 1
  )
  cluster_deff <- ifelse(
    stat_type == "mean" && cluster_stats$icc != 0, ClusterDeff(
      mean = cluster_stats$mean, sd = cluster_stats$sd,
      icc = cluster_stats$icc, type = cluster_deff_type
    ), 1
  )
  deff <- imbalance_deff * ancova_deff * cluster_deff
  # Compute minimum required sample size (MRSS) or the minimum detectable effect
  # size (MDES) with design effect adjustments
  if (estimand == "mrss") {
    if (resp_type == "binary" && method == "asin") {
      mdes <- ArcsineTransf(mean1) - ArcsineTransf(mean0)
    } else {
      mdes <- switch(effect_type,
        absolute = effect / sd,
        relative = effect * stat0 / sd,
        effect_size = effect
      )
    }
    out <- sample_num_multiplier * (z_quantile_multiplier / mdes)^2 * deff
  } else if (estimand == "mde") {
    # Compute the MDES with design effect
    mdes <- sqrt(sample_num_multiplier / n) * z_quantile_multiplier * sqrt(deff)
    # Compute the MDE for various effect sizes
    if (effect_type == "effect_size") {
      out <- mdes
    } else {
      if (resp_type == "binary") {
        # Note that for the conditional method, estimation = "numerical"
        # regardless of the value of the estimation argument.
        out <- MdePropTest(
          n = n, p0 = mean0, sig_level = sig_level, power = power, deff = deff,
          alternative = alternative, num_samples = num_samples, method = method,
          effect_type = effect_type, analysis = analysis,
          estimation = "analytic"
        )
      } else if (resp_type == "continuous") {
        sd_factor <- switch(effect_type,
          absolute = sd,
          relative = sd / stat0
        )
        out <- mdes * sd_factor
      }
    }
  }
  return(out)
}
PowerAnalysis <- function(...) {
  return(PowerPrecisionAnalysis(analysis = "power", ...))
}
PrecisionAnalysis <- function(...) {
  return(PowerPrecisionAnalysis(analysis = "precision", ...))
}
MinRequiredSampleSize <- function(...) {
  return(PowerPrecisionAnalysis(estimand = "mrss", ...))
}
MinSampleSize <- function(...) {
  return(PowerPrecisionAnalysis(estimand = "mrss", ...))
}
MRSS <- function(...) {
  return(PowerPrecisionAnalysis(estimand = "mrss", ...))
}
MinDetectableEffect <- function(...) {
  return(PowerPrecisionAnalysis(estimand = "mde", ...))
}
MDE <- function(...) {
  return(PowerPrecisionAnalysis(estimand = "mde", ...))
}

#' Power Analysis for the Binomial Proportion Test
#'
#' An enhanced version of the \code{power.prop.test} function in base R, which
#' according to its documentation functions to: "Compute the power of the
#' two-sample test for proportions, or determine parameters to obtain a target
#' power." The \code{power.prop.test} function uses the conditional method,
#' but \code{PowerPropTest} function enhances it to allow for the unconditional
#' method, as well as a one-sample test. The code for the most part is taken
#' from \code{power.prop.test} but with some modifications for syntax and
#' implementing the unconditional method. Similarly, the following documentation
#' has been taken almost verbatim from that of \code{power.prop.test}.
#'
#' @param n The number of observations (per group).
#' @param p1 The probability in one group for a two-sample analysis, or the
#'   probability under the null hypothesis for a one-sample analysis.
#' @param p2 The probability in other group for a two-sample analysis, or the
#'   probability under the alternative hypothesis for a one-sample analysis.
#' @param sig_level The significance level (Type I error probability).
#' @param power The power of the test (1 minus Type II error probability).
#' @param deff The design effect applied to the sample size. Defaults to 1.
#' @param alternative Either one- or two-sided test. Can be abbreviated.
#' @param num_samples Either one- or two-sample analysis.
#' @param strict Use strict interpretation in two-sided case.
#' @param method The method for inference, either conditional or unconditional.
#' @param tol The numerical tolerance used in root finding, the default
#'   providing (at least) four significant digits.
#' @return Object of class "power.htest", a list of the arguments (including the
#'   computed one) augmented with method and note elements.
#' @export

PowerPropTest <- function(n = NULL, p1 = NULL, p2 = NULL, sig_level = 0.05,
                          power = NULL, deff = 1, alternative = c(
                            "two_sided", "one_sided"
                          ), num_samples = c("two_sample", "one_sample"),
                          strict = FALSE, method = c(
                            "conditional", "unconditional"
                          ), tol = .Machine$double.eps^0.25) {
  if (sum(sapply(list(n, p1, p2, power, sig_level), is.null)) != 1) {
    stop(
      "exactly one of 'n', 'p1', 'p2', 'power', and 'sig_level' must be NULL"
    )
  }
  if (!is.null(sig_level) && !is.numeric(sig_level) ||
    any(0 > sig_level | sig_level > 1)) {
    stop("'sig_level' must be numeric in [0, 1]")
  }
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  tside <- switch(alternative,
    one_sided = 1,
    two_sided = 2
  )
  # Define the expression for power
  p_body <- "
    qu <- qnorm(sig_level / tside, lower.tail = FALSE)
    d <- abs(p1 - p2)
    q1 <- 1 - p1
    q2 <- 1 - p2
    v1 <- p1 * q1
    v2 <- p2 * q2
  "
  if (method == "conditional") {
    if (num_samples == "two_sample") {
      p_body <- paste(p_body, "
        pbar <- (p1 + p2) / 2
        qbar <- 1 - pbar
        vbar <- pbar * qbar
      ")
      quantile_value <- paste(
        "(sqrt(n / deff) * d - qu * sqrt(2 * vbar)) /",
        "sqrt(v1 + v2)"
      )
    } else if (num_samples == "one_sample") {
      quantile_value <- "(sqrt(n / deff) * d - qu * sqrt(v1)) / sqrt(v2)"
    }
  } else if (method == "unconditional") {
    if (num_samples == "two_sample") {
      quantile_value <- "sqrt(n / deff) * d / sqrt(v1 + v2) - qu"
    } else if (num_samples == "one_sample") {
      quantile_value <- "sqrt(n / deff) * d / sqrt(v2) - qu"
    }
  }
  prob_value <- paste0("pnorm(", quantile_value, ")")
  if (strict && tside == 2) {
    prob_value <- paste0(
      prob_value, " + ",
      paste0("pnorm(", quantile_value, ", lower.tail = FALSE)")
    )
  }
  p_body <- paste(
    p_body, "
    ", prob_value
  )
  # Estimate the power analysis estimand
  if (is.null(power)) {
    power <- eval(parse(test = p_body))
  } else if (is.null(n)) {
    n <- uniroot(
      function(n) {
        return(eval(parse(text = p_body)) - power)
      },
      c(1, 1e+07),
      tol = tol, extendInt = "upX"
    )$root
  } else if (is.null(p1)) {
    p1 <- uniroot(
      function(p1) {
        return(eval(parse(text = p_body)) - power)
      },
      c(0, p2),
      tol = tol, extendInt = "yes"
    )$root
    if (p1 < 0) {
      warning("No p1 in [0, p2] can be found to achieve the desired power")
    }
  } else if (is.null(p2)) {
    p2 <- uniroot(
      function(p2) {
        return(eval(parse(text = p_body)) - power)
      },
      c(p1, 1),
      tol = tol, extendInt = "yes"
    )$root
    if (p2 > 1) {
      warning("No p2 in [p1, 1] can be found to achieve the desired power")
    }
  } else if (is.null(sig_level)) {
    sig_level <- uniroot(
      function(sig_level) {
        return(eval(parse(text = p_body)) - power)
      },
      c(1e-10, 1 - 1e-10),
      tol = tol, extendInt = "upX"
    )$root
    if (sig_level < 0 || sig_level > 1) {
      warning(
        "No significance level [0, 1] can be found to achieve the desired power"
      )
    }
  } else {
    stop("internal error", domain = NA)
  }
  out <- structure(
    list(
      n = n, p1 = p1, p2 = p2, sig.level = sig_level, power = power,
      alternative = alternative, note = "n is number in *each* group",
      method = "Two-sample comparison of proportions power calculation"
    ),
    class = "power.htest"
  )
  return(out)
}

#' Minimum Detectable Effect for the Binomial Proportion Test
#'
#' Computes the minimum detectable effect (MDE) for a binomial proportion test
#' with the following options:
#' (1) A relative effect or an absolute effect.
#' (2) The unconditional, conditional, and arcsine transformation methods.
#' (3) Estimation using analytic formulas or numerical root-finding procedures.
#' (4) Power analysis or precision analysis.
#' (5) Two- or one-sample analysis.
#' (6) Two- or one-sided test.
#'
#' @param n The number of observations (per group).
#' @param p0 The binomial proportion in the baseline group for a two-sample
#'   analysis, or a binomial proportion under the null hypothesis for a
#'   one-sample analysis.
#' @param sig_level The significance level (Type I error probability).
#' @param power The power of the test (1 minus Type II error probability).
#' @param deff The design effect applied to the sample size. Defaults to 1.
#' @param alternative Either one- or two-sided test. Can be abbreviated.
#' @param num_samples Either one- or two-sample analysis.
#' @param effect_type Whether a relative effect or an absolute effect.
#' @param method The method for inference, either the unconditional,
#'   conditional, or arcsine transformation method.
#' @param estimation The estimation procedure, either using an analytic formula
#'   or a numerical root-finding procedure. The conditional method currently
#'   implements the numerical procedure only, without the option for an
#'   analytic formula.
#' @param analysis The analysis type, whether a power analysis or a precision
#'   analysis. Note that the precision analysis is currently available only for
#'   the unconditional and arcsine methods using the analytic estimation
#'   procedure.
#' @return The MDE, either the absolute MDE or the minimum detectable relative
#'   effect (MDRE).
#' @export

MdePropTest <- function(n = NULL, p0 = NULL, sig_level = 0.05,
                        power = 0.80, deff = 1,
                        alternative = c("two_sided", "one_sided"),
                        num_samples = c("two_sample", "one_sample"),
                        effect_type = c("relative", "absolute"),
                        method = c("unconditional", "conditional", "asin"),
                        estimation = c("analytic", "numerical"),
                        analysis = c("power", "precision")) {
  alternative <- match.arg(alternative)
  num_samples <- match.arg(num_samples)
  effect_type <- match.arg(effect_type)
  method <- match.arg(method)
  estimation <- match.arg(estimation)
  analysis <- match.arg(analysis)
  z_quantile_multiplier <- ZQuantileMultiplier(
    analysis = analysis, alternative = alternative, sig_level = sig_level,
    power = power
  )
  if (method == "unconditional") {
    if (estimation == "analytic") {
      a0 <- n / (deff * (z_quantile_multiplier)^2)
      A <- (a0 + 1) * p0
      B <- 2 * p0 - 1
      C <- (p0 - 1) * ifelse(num_samples == "two_sample", 2, 1)
      mdre <- (-B + sqrt(B^2 - 4 * A * C)) / (2 * A) # positive root
    } else if (estimation == "numerical") {
      mdre <- PowerPropTest(
        n = n, p1 = p0, p2 = NULL, sig_level = sig_level,
        power = power, alternative = alternative, num_samples = num_samples,
        method = method, strict = FALSE, deff = deff
      )$p2 / p0 - 1
    }
  } else if (method == "conditional") {
    mdre <- PowerPropTest(
      n = n, p1 = p0, p2 = NULL, sig_level = sig_level,
      power = power, alternative = alternative, num_samples = num_samples,
      method = method, strict = FALSE, deff = deff
    )$p2 / p0 - 1
  } else if (method == "asin") {
    if (estimation == "analytic") {
      sample_num_multiplier <- ifelse(num_samples == "two_sample", 2, 1)
      mdes <- sqrt(sample_num_multiplier / n) * z_quantile_multiplier *
        sqrt(deff)
    } else if (estimation == "numerical") {
      pwr_fun <- switch(num_samples,
        two_sample = pwr::pwr.2p.test,
        one_sample = pwr::pwr.p.test
      )
      alternative_dot <- switch(alternative,
        two_sided = "two.sided",
        one_sided = "greater"
      )
      mdes <- pwr_fun(
        h = NULL, n = n / deff, power = power, sig.level = sig_level,
        alternative = alternative_dot
      )$h
    }
    mdre <- InvArcsineTransf(ArcsineTransf(p0) + mdes) / p0 - 1
  }
  out <- switch(effect_type,
    absolute = mdre * p0,
    relative = mdre
  )
  return(out)
}

#' Z-Quantile Multiplier
#'
#' Helper function to compute the multiplier of z-quantiles, i.e., quantiles of
#' the standard normal distribution, for computing the minimum detectable
#' effect (MDE) in a power analysis or the maximum margin or error (MME) in a
#' precision analysis, assuming asymptotic normality of the test statistic,
#' with options for two- and one-sided alternatives. This multiplier is squared
#' for computing the minimum required sample size (MRSS) in power or precision
#' analysis.
#'
#' @param analysis The type of analysis conducted, either \code{"power"} for
#'   power analysis or \code{"precision"} for precision analysis.
#' @param alternative Whether the underlying hypothesis test for the power
#'   precision analysis is one- or two-sided. Defaults to \code{"two_sided"} for
#'   a two-sided test.
#' @param sig_level The significance level for the power or precision analysis.
#'   Defaults to \code{0.05}.
#' @param power The minimum power required for power analysis. Defaults to
#'   \code{0.80}. Not applicable for precision analysis.
#' @return The z-quantile multiplier for power or precision analysis.
#' @references
#' \itemize{
#' \item Chow, Shien-Chung, Jun Shao, and Hansheng Wang (2003). \emph{Sample
#'   Size Calculations in Clinical Research.} Marcel Dekker: New York.
#' }
#' @export

ZQuantileMultiplier <- function(analysis = c("power", "precision"),
                                alternative = c("two_sided", "one_sided"),
                                sig_level = 0.05, power = 0.80) {
  analysis <- match.arg(analysis)
  z_upper_sig <- ZUpperSig(alternative = alternative, sig_level = sig_level)
  z_beta <- qnorm(power)
  out <- switch(analysis,
    power = z_upper_sig + z_beta,
    precision = z_upper_sig
  )
  return(out)
}

#' Z Upper Quantile for Significance
#'
#' Helper function to compute the standard normal upper-quantile for a
#' specified significance level \code{sig_level}:
#' (1) \code{qnorm(1 - sig_level)} for a one-sided hypothesis test, and
#' (2) \code{qnorm(1 - sig_level / 2)} for a two-sided hypothesis test.
#'
#' @param alternative Whether the underlying hypothesis test for the power
#'   precision analysis is one- or two-sided. Defaults to \code{"two_sided"} for
#'   a two-sided test.
#' @param sig_level The significance level for the power or precision analysis.
#'   Defaults to \code{0.05}.
#' @return The z upper-quantile for the significance level.
#' @export

ZUpperSig <- function(alternative = c("two_sided", "one_sided"),
                      sig_level = 0.05) {
  alternative <- match.arg(alternative)
  out <- qnorm(
    sig_level / ifelse(alternative == "two_sided", 2, 1),
    lower.tail = FALSE
  )
  return(out)
}

#' Imbalance Design Effect
#'
#' Compute the design effect (relative variance) from an imbalance in number of
#' observations between two study arms or groups, typically, Treatment and
#' Control, relative to a balanced design.
#'   Given \code{phi} as the proportion in one of the two arms, the design
#' effect is given by \code{1 / (4 * phi * (1 - phi))}, which is applied as a
#' post-hoc multiplier in two ways:
#' (1) The design effect is multiplied to the minimum required sample size
#'     (MRSS), which itself is a relative variance measure.
#' (2) The square root of the design effect, a measure of relative standard
#'     standard deviation, is multiplied to the minimum detectable effect (MDE)
#'     and the maximum margin of error (MME), which themselves are relative
#'     standard deviation measures.
#'
#' @param phi The proportion of observations in one of the arms, either
#'   Treatment or Control, for a two-sample analysis. Defaults to \code{0.50},
#'   i.e., a balanced allocation for both arms.
#' @return The imbalance design effect for the study with two arms.
#' @export

ImbalanceDeff <- function(phi = 0.50) {
  stopifnot(phi >= 0, phi <= 1)
  out <- 1 / (4 * phi * (1 - phi))
  return(out)
}

#' Analysis-of-Covariance (ANCOVA) Design Effect
#'
#' For a two-sample design, compute the design effect from using an analysis
#' of covariance (ANCOVA) relative to an analysis of variance (ANOVA) design,
#' using the approximation of Borm et al. (2007). The ANCOVA design is analyzed
#' using a regression model on treatment and additional covariates, while the
#' ANOVA design is analyzed using regression model on treatment assignmen
#' alone.
#'   The design effect can be used in a one-sample design as well, even though
#' the name ANCOVA would typically not be used in this setting. We essentially
#' compare the error variance of a regression model on covariates with the
#' error variance of a regression model with no covariates, i.e., an
#' intercept-only model, to obtain the design effect.
#'   The design effect is the multiplier for the minimum required sample size
#' (MRSS), and the square root of the design effect is the multiplier for
#' the minimum detectable effect (MDE) and the maximum margin of error (MME).
#' See the documentation for \code{ImbalanceDeff} for details.
#'   Note that for a two-sample analysis, there could be imbalances between
#' the Treatment and Control groups, e.g., due to the study being an
#' observational without randomized group assignment, or due to a low sample
#' size. In this setting, the design effect needs to be multiplied by the
#' variance inflation factor (VIF) for the treatment effect, which is a function
#' of the R-squared of the regression of the treatment variable on all
#' covariates.
#'
#' @param r The multiple correlation from the regression of the response
#'   on covariates, not including the treatment variable, in either the two- or
#'   one-sample design. If there is only one covariate, this reduces to the
#'   correlation between the response and the covariate. For a two-sample
#'   analysis, this should be estimated from the Control or baseline group data
#'   only, so that there is no bias is estimating this parameter from the
#'   presence of a treatment effect.
#' @param r_treat The multiple correlation from the regression of the
#'   treatment variable on all the covariates, applicable only in a two-sample
#'   analysis. This is typically close to zero and negligible for a randomized
#'   experiment. But for an observational study, this could be large enough
#'   so that we need to account for it.
#' @return The ANCOVA design effect for the study design.
#' @references
#' \itemize{
#' Borm, George F., Jaap Fransen, and Wim A. J. G. Lemmens (2007). "A Simple
#'   Sample Size Formula for Analysis of Covariance in Randomized Clinical
#'   Trials." \emph{Journal of Clinical Epidemiology,} 60: 1234-1238.
#' }
#' @export

AncovaDeff <- function(r = 0, r_treat = 0) {
  stopifnot(r >= 0, r <= 1, r_treat >= 0, r_treat <= 1)
  out <- (1 - r^2) / (1 - r_treat^2)
  return(out)
}

#' Cluster Design Effect
#'
#' Estimate the design effect in clustered data from having within-cluster
#' correlation relative to a zero-correlation design, with the following two
#' options:
#' (1) Minimum design effect, given by \code{1 + (mean - 1) * icc}, where
#'     \code{mean} is the mean of the cluster sizes and \code{icc}, the
#'     intra-cluster correlation (ICC), i.e., the within-cluster correlation.
#' (2) Maximum design effect, given by
#'     \code{1 + ((cv^2 * (1 - 1 / k) + 1) * mean - 1) * icc}, where
#'     \code{cv} is the coefficient of variation (CV) for the cluster sizes and
#'     \code{k}, the number of clusters. If \code{k} is \code{NULL}, we assume
#'     that \code{k} is \code{Inf}, so that the design effect becomes
#'     \code{1 + (cv^2 * mean - 1) * icc}.
#' Eldridge et al. (2006) provided the above lower and upper bounds of the
#' cluster design effect, because the actual cluster sizes are typically not
#' known before study begins.
#'
#' @param cluster_size A numeric vector of cluster sizes for the study.
#' @param mean The mean of the cluster sizes. Used if
#'   \code{cluster_size} is \code{NULL}.
#' @param sd The standard deviation (SD) of the cluster sizes, mainly for
#'   computing the CV, i.e., the SD divided by the mean. Used if
#'   \code{cluster_size} is \code{NULL}.
#' @param icc The intra-cluster correlation (ICC)
#' @param k The number of clusters.
#' @param type Either the minimum design effect (\code{"min"}) or the maximum
#'   (\code{"max"}). TODO(shyuemingloh@): Implement other design effect types
#'   when the actual cluster sizes are provided.
#' @return The cluster design effect for the study.
#' @references
#' \itemize{
#' Eldridge, Sandra M., Deborah Ashby, and Sally Kerry (2006). "Sample Size for
#'   Cluster Randomized Trials: Effect of Coefficient of Variation of Cluster
#'   Size and Analysis Method." \emph{International Journal of Epidemiology,}
#'   35(5): 1292-1300.
#' }
#' @export

ClusterDeff <- function(cluster_size = NULL, mean = 1, sd = 0, k = NULL,
                        icc = 0, type = c("min", "max")) {
  type <- match.arg(type)
  if (!is.null(cluster_size)) {
    mean <- mean(cluster_size)
    sd <- sd(cluster_size)
  }
  mean_multiplier <- switch(type,
    min = 1,
    max = {
      cv <- sd / mean
      cv_multiplier <- ifelse(is.null(k), 1, 1 - 1 / k)
      cv^2 * cv_multiplier + 1
    }
  )
  out <- 1 + (mean * mean_multiplier - 1) * icc
  return(out)
}

#' Sample Size and Duration Conversion
#'
#' Convert from sample size to duration or vice versa, i.e., duration to
#' sample size, based on an assumed linear relationship between the two,
#' with parameters determine from interpolating and extrapolating a line
#' through the origin and the sample size at a nonzero time point.
#'
#' @param n A numeric vector of sample sizes. If not \code{NULL}, then the
#'   sample size values will be converted to duration values.
#' @param duration A numeric vector of duration in weeks. If \code{n} is
#'   \code{NULL} and \code{duration} is not \code{NULL}, the duration values
#'   will be converted to sample size values.
#' @param n0 The sample size at a particular nonzero time point, i.e., at
#'   \code{t0}. In the language of Cloud experiments, this is the exposure
#'   traffic size within a \code{t0}-week period. Defaults to 10,000 units.
#' @param t0 The nonzero time point to observe the sample size \code{n0}.
#'   Defaults to 4 weeks.
#' @param arms The number of arms for the experiment. Defaults to 2 arms.
#' @param traffic_prop The proportion of exposure traffic available for
#'   experimentation. Defaults to \code{MendelStudyTrafficProportion()}. See
#'   documentation for \code{TrafficProportion} for details.
#' @param round_up Whether to round up the returned value.
#' @return A numeric vector of sample size or duration values, depending on
#'   on the values of \code{n} and \code{duration}.
#' @export

SampleSizeDurationConversion <- function(n = NULL, duration = NULL, n0 = 1e4,
                                         t0 = 4, arms = 2,
                                         traffic_prop =
                                           StudyTrafficProportion(),
                                         round_up = FALSE) {
  k0 <- traffic_prop * n0 / (t0 * arms)
  if (!is.null(n)) {
    out <- n / k0
  } else if (!is.null(duration)) { # nolint
    out <- k0 * duration
  } else {
    out <- NULL
  }
  if (round_up) {
    out <- ceiling(out)
  }
  return(out)
}
SampleSizeToDuration <- function(n, ...) {
  stopifnot(!is.null(n))
  out <- SampleSizeDurationConversion(n = n, duration = NULL, ...)
  return(out)
}
DurationToSampleSize <- function(duration, ...) {
  stopifnot(!is.null(duration))
  out <- SampleSizeDurationConversion(duration = duration, n = NULL, ...)
  return(out)
}

#' Traffic Proportion for Study and Launch
#'
#' The proportion of traffic available to a study—i.e., an A/B experiment—or a
#' launch. In the context of many online experiments, including those at Google,
#' a portion of the traffic is held back from experimentation, typically 15%, so
#' that the proportion available for experimentation in a study is typically
#' 85%. For a launch, the proportion is 100%.
#'
#' @param type Whether a study (\code{"study"}) or launch (\code{"launch"}).
#' @return The proportion of traffic available for experimentation.
#' @export

TrafficProportion <- function(type = c("study", "launch")) {
  type <- match.arg(type)
  out <- switch(type,
    study = 0.85,
    launch = 1
  )
  return(out)
}
StudyTrafficProportion <- function() {
  return(TrafficProportion(type = "study"))
}
LaunchTrafficProportion <- function() {
  return(TrafficProportion(type = "launch"))
}

#' Power and Precision Analysis Data
#'
#' Preprocess data for power or precision analysis, with the following key
#' variables:
#' (1) Name of conversion metric.
#' (2) Response type, whether binary or continuous.
#' (3) Statistic type, whether mean or a quantile.
#' (4) Number of observations for the power or precision analysis.
#' (5) Value of the statistic for the analysis, either mean or quantile.
#' (6) Population standard deviation, without covariate adjustment.
#' (7) Population standard deviation, with covariate adjustment.
#' (8) Multiple correlation for the regression of the metric on the covariates,
#'     available only for an analysis using the mean statistic.
#' (9) Optionally, statistics for estimating the design effect for clustered
#'     data, including the mean and standard deviation of the cluster sizes
#'     and the estimated intra-cluster correlation (ICC) for the metric given
#'     the clusters.
#'   The arguments \code{y_binary}, \code{y_continuous_mean}, and
#' \code{y_continuous_quantile} take named list of numeric vector containing
#' the names of, respectively, binary metric, continuous metrics using the
#' mean statistic, and continuous metrics using the quantile statistic.
#' The list names are names of exposure variables. The exposure name specified
#' in the argument \code{all_exposure} is a keyword for using all observations.
#'   Note that we have three different names for this function, i.e.,
#' \code{PowerPrecisionAnalysisData}, \code{PowerAnalysisData}, and
#' \code{PrecisionAnalysisData}. All three are identical, since the input data
#' prepared for power analysis is the same as that for precision analysis.
#'
#' @param data A data frame containing disaggregated data to estimate inputs
#'   for power or precision analysis.
#' @param y_binary A named list of numeric vectors containing the names of the
#'   binary metrics.
#' @param y_continuous_mean A named list of numeric vectors containing the names
#'   of continuous metrics that use a mean statistic.
#' @param y_continuous_quantile A named list of numeric vector containing the
#'   names of continuous metrics that use a quantile statistic.
#' @param ancova_deff Whether the statistics to estimate the analysis of
#'   covariance (ANCOVA) design effect, namely, the regression-adjusted standard
#'   deviation \code{sd_adj} and the regression multiple correlation \code{r},
#'   should be computed.
#' @param x A numeric vector containing the names of all covariates used
#'   for the regression adjustment of all metrics.
#' @param include_pre Whether to include the pre-period version of the metrics
#'   as covariates in addition to those specified in \code{x}.
#' @param post_suffix The suffix used for the metrics to indicate that they
#'   measured in the post-period. These would be replaced by \code{pre_suffix}
#'   to obtain the names of the pre-period metrics if
#'   \code{include_pre = TRUE}.
#' @param pre_suffix The suffix used for the names of the pre-period versions
#'   of the metrics.
#' @param quantile_probs The quantile probabilities for the quantiles to
#'   conduct power or precision analysis for.
#' @param quantile_type The type of quantile estimated for the \code{stat_value}
#'   column of the data frame, basically the value passed into the \code{type}
#'   argument of the \code{stats::quantile} function.
#' @param cluster_deff Whether the statistics to estimate the cluster design
#'   effect, namely the cluster size mean \code{cs_mean}, the cluster size
#'   standard deviation \code{cs_sd}, and the intracluster correlation
#'   \code{icc}, should be computed.
#' @param cluster_id String name of the cluster ID variable.
#' @param all_exposure String keyword to refer to using all observations as
#'   the exposure condition.
#' @param ... Additional arguments passed into the \code{QuantileSD} function
#'   to estimate the quantile population SD.
#' @return A data frame containing the key columns mentioned in the description
#'   above.
#' @export

PowerPrecisionAnalysisData <- function(data, y_binary = NULL,
                                       y_continuous_mean = NULL,
                                       y_continuous_quantile = NULL,
                                       ancova_deff = FALSE, x = NULL,
                                       include_pre = FALSE,
                                       post_suffix = "_post",
                                       pre_suffix = "_pre",
                                       cluster_deff = FALSE,
                                       cluster_id = NULL,
                                       quantile_probs = c(
                                         0.5, 0.90, 0.95, 0.99
                                       ), quantile_type = 7,
                                       all_exposure = "all",
                                       ...) {
  # Check for outcomes
  if (is.null(y_binary) && is.null(y_continuous_mean) &&
    is.null(y_continuous_quantile)) {
    stop(
      "One of y_binary or y_continuous_mean or y_continuous_quantile ",
      "must be non-empty."
    )
  }
  # Define vector of all outcomes with accompanying information on exposure,
  # response type, statistic type, and pre-period metrics.
  if (is.null(y_binary)) {
    exposure_binary <- NULL
  } else {
    exposure_binary <- rep(names(y_binary), times = lengths(y_binary))
    y_binary <- as.character(unlist(y_binary))
  }
  if (is.null(y_continuous_mean)) {
    exposure_continuous_mean <- NULL
  } else {
    exposure_continuous_mean <- rep(
      names(y_continuous_mean),
      times = lengths(y_continuous_mean)
    )
    y_continuous_mean <- as.character(unlist(y_continuous_mean))
  }
  if (is.null(y_continuous_quantile)) {
    exposure_continuous_quantile <- NULL
  } else {
    exposure_continuous_quantile <- rep(
      names(y_continuous_quantile),
      times = lengths(y_continuous_quantile)
    )
    y_continuous_quantile <- as.character(unlist(y_continuous_quantile))
  }
  y <- c(y_binary, y_continuous_mean, y_continuous_quantile)
  exposure <- c(
    exposure_binary, exposure_continuous_mean, exposure_continuous_quantile
  )
  resp_type <- c(
    rep("binary", length(y_binary)),
    rep("continuous", length(c(y_continuous_mean, y_continuous_quantile)))
  )
  stat_type <- c(
    rep("mean", length(c(y_binary, y_continuous_mean))),
    rep("quantile", length(y_continuous_quantile))
  )
  if (include_pre) {
    y_pre <- gsub(post_suffix, pre_suffix, y, fixed = TRUE)
  } else {
    y_pre <- NULL
  }
  # Show at least the median if there are continuous metrics with quantile
  # statistic but no quantile probabilities specified
  if (is.null(quantile_probs) && !is.null(y_continuous_quantile)) {
    quantile_probs <- 0.5
  }
  # Estimate cluster mean and sd for all exposure conditions if cluster
  # information is provided.
  if (!is.null(cluster_id)) {
    cluster_size <- list()
    unique_exposures <- unique(exposure)
    for (i in seq_along(unique_exposures)) {
      if (unique_exposures[i] == all_exposure) {
        exposure_data <- data
      } else {
        exposure_data <- data[data[[unique_exposures[i]]] == 1, ]
      }
      cluster_data <- exposure_data %>%
        group_by_at(cluster_id) %>%
        summarize(cluster_size = n())
      cluster_size[[unique_exposures[i]]] <- list(
        mean = mean(cluster_data$cluster_size),
        sd = sd(cluster_data$cluster_size)
      )
    }
  }
  # Iterate through all outcomes to compute relevant stats
  out <- list()
  for (i in seq_along(y)) {
    if (exposure[i] == all_exposure) {
      data_i <- data
    } else {
      data_i <- data[data[[exposure[i]]] == 1, ]
    }
    if (stat_type[i] == "mean") {
      stat <- "mean"
      value <- mean(data_i[[y[i]]])
      tau <- NA
    } else {
      stat <- paste0("q", 100 * quantile_probs)
      value <- quantile(
        data_i[[y[i]]],
        probs = quantile_probs, type = quantile_type
      ) %>%
        as.numeric()
      tau <- quantile_probs
    }
    stats <- data.frame(
      exposure = exposure[i],
      metric = y[i],
      resp_type = resp_type[i],
      stat = stat,
      tau = tau,
      n = nrow(data_i),
      stat_value = value,
      sd = 0
    )
    if (ancova_deff) {
      stats <- cbind(stats, sd_adj = NA, r = NA)
    }
    if (cluster_deff) {
      stats <- cbind(stats, cs_mean = NA, cs_sd = NA, icc = NA)
    }
    covars <- union(y_pre[i], x)
    for (j in seq(nrow(stats))) {
      if (stats[j, "stat"] == "mean") {
        stats[j, "sd"] <- sd(data_i[[y[i]]])
        if (!is.null(covars)) {
          m <- lm(CreateFormula(y[i], covars), data = data_i)
          stats[j, "sd_adj"] <- summary(m)$sigma
          stats[j, "r"] <- sqrt(summary(m)$r.squared)
        }
        if (!is.null(cluster_id)) {
          stats$cs_mean <- cluster_size[[exposure[i]]]$mean
          stats$cs_sd <- cluster_size[[exposure[i]]]$sd
          stats$icc <- ICC(data_i, y = y[i], cluster_id = cluster_id)
        }
      } else {
        stats[j, "sd"] <- QuantileSD(
          data_i[[y[i]]],
          tau = stats[j, "tau"], residualize = TRUE, ...
        )
        if (!is.null(covars)) {
          m <- quantreg::rq(
            CreateFormula(y[i], covars),
            tau = stats[j, "tau"], data = data_i
          )
          stats[j, "sd_adj"] <- QuantileSD(m, tau = stats[j, "tau"], ...)
          stats[j, "r"] <- NA
        }
      }
    }
    out[[i]] <- stats
  }
  # Create output data frame with variables defined for percentizing
  percent_vars <- c("stat_value", "sd")
  if (ancova_deff) {
    percent_vars <- c(percent_vars, "sd_adj", "r")
  }
  if (cluster_deff) {
    percent_vars <- c(percent_vars, "icc")
  }
  out <- do.call(rbind, out) %>%
    dplyr::select(-tau) %>%
    Percentize(x = percent_vars, digits = 2)
  return(out)
}
PowerAnalysisData <- function(...) {
  return(PowerPrecisionAnalysisData(...))
}
PrecisionAnalysisData <- function(...) {
  return(PowerPrecisionAnalysisData(...))
}

#' Power and Precision Analysis Table
#'
#' Create table of power/precision analysis results, either a table of
#' minimum durations, or a table of minimum detectable effects (MDEs) for a
#' power analysis or a table of maximum margins of error (MMEs) for a precision
#' analysis.
#'
#' @param data A data frame containing inputs for the power analysis, i.e.,
#'   columns for the metric names and the metric means, standard deviations, and
#'   ANCOVA/pre-post correlations.
#' @param exposure A string containing the name of the column of exposure
#'   variables.
#' @param metric A string containing the name of the column of names for
#'   various conversion metrics.
#' @param resp_type A string containing the name of the column of defining the
#'   response type, whether binary or continuous. See documentation for
#'   \code{PowerPrecisionAnalysis} for details.
#' @param stat_type A string containing the name of the column of defining the
#'   statistic type, whether a mean or a quantile. See documentation for
#'   \code{PowerPrecisionAnalysis} for details.
#' @param stat0 A string containing the name of the column of statistics for the
#'   conversion metrics.
#' @param sd0 A string containing the name of the column of standard deviations
#'   for the conversion metrics.
#' @param cor0 A string containing the name of the column of multiple
#'   correlation R values for the regression of the corresponding conversion on
#'   covariates. These could be the correlations between the post-period and
#'   pre-period metrics for the conversion for the regression of the post-period
#'   conversion on the pre-period conversion. But they could be the square-root
#'   of the R-squared values for the regression of the conversion on any set of
#'   covariates.
#' @param n0 A string containing the name of the column for the number of
#'   exposures (traffic size) for fixed time period, e.g., 4 weeks.
#' @param table_type A string specifying the table type, whether
#'   \code{"duration"} for a table of minimum durations or \code{"mde"} for a
#'   table of MDEs or MMEs.
#' @param analysis The type of analysis conducted, either \code{"power"} for
#'   power analysis or \code{"precision"} for precision analysis. This currently
#'   applies only to the duration table, i.e., \code{table_type = "duration"},
#'   since \code{table_type = "mde"} and \code{table_type = "mme"} will
#'   specify the type of analysis by their argument values alone, i.e., a power
#'   analysis for \code{table_type = "mde"} and a precision analysis for
#'   \code{table_type = "mme"}.
#' @param effect A numeric vector of MDE values to determine the corresponding
#'   MRSS values in a power analysis or MME values for a precision analysis.
#' @param duration A numeric vector of durations (in weeks) to determine the
#'   corresponding MDE values in a power analysis or a precision analysis.
#' @param suffix A numeric vector of suffix values for the duration/MDE column
#'   names. Defaults to the duration/MDE sequence of values if \code{NULL}.
#' @param round_up Whether to round up the durations in a table of durations.
#' @param ancova_deff Whether the analysis of covariance (ANCOVA) design effect
#'   should be applied.
#' @param cluster_deff Whether the cluster design effect should be applied.
#' @param cluster_deff_type The type of cluster design effect to apply, whether
#'   \code{"min"} for the lower bound or \code{"max"} for the upper bound.\
#'   See documentation for the \code{ClusterDeff} function for details.
#' @param cs_mean Name of column with cluster size means.
#' @param cs_sd Name of column with cluster size standard deviations.
#' @param icc Name of column containing estimated intra-cluster correlations
#'   (ICCs).
#' @param arms The number of arms in the study.
#' @param ... Additional arguments to pass into the \code{MRSS} or \code{MDE}
#'   functions.
#' @return A data frame containing the duration or MDE table, with percent
#'   format in relevant columns.
#' @export

PowerPrecisionAnalysisTable <- function(data, exposure = "exposure",
                                        metric = "metric",
                                        resp_type = "resp_type",
                                        stat_type = "stat", n0 = "n",
                                        stat0 = "stat_value", sd0 = "sd",
                                        cor0 = "r", table_type = c(
                                          "duration", "mde", "mme"
                                        ), analysis = c("power", "precision"),
                                        effect = seq(0.02, 0.08, by = 0.02),
                                        duration = seq(4, 16, by = 4),
                                        suffix = NULL, round_up = TRUE,
                                        ancova_deff = FALSE,
                                        cluster_deff = FALSE,
                                        cluster_deff_type = c("min", "max"),
                                        cs_mean = "cs_mean", cs_sd = "cs_sd",
                                        icc = "icc", arms = 2, ...) {
  table_type <- match.arg(table_type)
  analysis <- match.arg(analysis)
  select_vars <- c(exposure, metric, resp_type, stat_type, n0, stat0, sd0)
  percent_vars <- c(stat0, sd0)
  if (ancova_deff) {
    select_vars <- c(select_vars, cor0)
    percent_vars <- c(percent_vars, cor0)
  }
  out <- data %>%
    select_at(select_vars)
  out[[resp_type]] <- as.character(out[[resp_type]])
  out[[stat_type]] <- as.character(out[[stat_type]])
  # Determine duration/MDE/MME sequence of values, column suffix, and column
  # names
  if (table_type == "duration") {
    sequence <- effect
  } else {
    sequence <- duration
    analysis <- switch(table_type,
      mde = "power",
      mme = "precision"
    )
  }
  if (is.null(suffix)) {
    suffix <- sequence
  }
  # Initialize columns and values
  col_names <- paste0(table_type, "_", suffix)
  for (j in seq_along(col_names)) {
    out[[col_names[j]]] <- 0
  }
  # Compute values for duration/MDE/MME columns
  for (i in seq(nrow(out))) {
    stat_type_i <- ifelse(
      startsWith(out[i, stat_type], "q"), "quantile", out[i, stat_type]
    )
    cor_i <- ifelse(ancova_deff, out[i, cor0], 0)
    if (cluster_deff) {
      cluster_stats <- list(
        mean = data[i, cs_mean], sd = data[i, cs_sd], icc = data[i, icc]
      )
    } else {
      cluster_stats <- list(mean = 1, sd = 0, icc = 0)
    }
    for (j in seq_along(sequence)) {
      if (table_type == "duration") {
        out[i, col_names[j]] <- MRSS(
          stat0 = out[i, stat0], sd0 = out[i, sd0], cor0 = cor_i,
          resp_type = out[i, resp_type], stat_type = stat_type_i,
          effect = sequence[j], analysis = analysis,
          cluster_stats = cluster_stats, cluster_deff_type = cluster_deff_type,
          ...
        ) %>%
          SampleSizeToDuration(
            n0 = out[i, n0], arms = arms, round_up = round_up
          )
      } else {
        out[i, col_names[j]] <- MDE(
          stat0 = out[i, stat0], sd0 = out[i, sd0], cor0 = cor_i,
          resp_type = out[i, resp_type], stat_type = stat_type_i,
          n = DurationToSampleSize(
            duration = sequence[j], n0 = out[i, n0], arms = arms
          ), analysis = analysis, cluster_stats = cluster_stats,
          cluster_deff_type = cluster_deff_type, ...
        )
      }
    }
  }
  # Use percent format for relevant table columns
  if (table_type != "duration") {
    percent_vars <- c(percent_vars, col_names)
  }
  out <- Percentize(out, x = percent_vars)
  return(out)
}
DurationTable <- function(...) {
  return(PowerPrecisionAnalysisTable(table_type = "duration", ...))
}
MdeTable <- function(...) {
  return(PowerPrecisionAnalysisTable(table_type = "mde", ...))
}
MaxErrorMarginTable <- function(...) {
  return(PowerPrecisionAnalysisTable(table_type = "mme", ...))
}
