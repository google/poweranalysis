% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{StartPars}
\alias{StartPars}
\title{Start Values of Log-likelihood Parameters}
\usage{
StartPars(
  y,
  dist = c("norm", "logis", "cauchy", "lnorm", "gamma"),
  sd_mle = TRUE
)
}
\arguments{
\item{y}{A numeric vector containing the observed data for estimating the
MLE of the parameters.}

\item{dist}{One of the distribution keyword strings specified above.}

\item{sd_mle}{Whether to use the MLE of the standard deviation and variance
for the sample standard deviation in the start values. If \code{TRUE},
the typical sample standard deviation and variance estimators are
multiplied by, respectively, the factors \code{sqrt((n - 1) / n} and
\code{(n - 1) / n}. If \code{FALSE}, the typical standard deviation and
variance estimators are used, where the latter is an unbiased estimator of
of the population variance.}
}
\value{
A vector containing the start values for estimating the MLE of the
  parameters under the assumed parametric model.
}
\description{
Given a numeric vector \code{y}, estimate the start values of the parameters
for maximizing the log-likelihood of \code{y} for one of the following
parametric distributions specified for the \code{LogLikFun} function:
(1) Normal (\code{"norm"}) and logistic (\code{"logis"}): The sample mean and
    standard deviation of \code{y}.
(2) Cauchy (\code{"cauchy"}): The sample median and interquartile range (IQR)
    of \code{y}.
(3) Lognormal (\code{"lnorm"}): The sample mean and standard deviation of
    \code{log(y)}.
(4) Gamma (\code{"gamma"}): The method-of-moments (MOM) estimators of the
    shape and rate parameters of \code{y}, based on the sample mean and
    standard deviation of \code{y}.
}
\details{
For the normal and lognormal distributions, the start value of the first
parameter, the sample mean, is exactly the maximum likelihood estimate (MLE),
while the start value of the second parameter, the sample standard deviation,
is exactly the MLE if \code{sd_mle = TRUE} and asymptotically the MLE if
\code{sd_mle = FALSE}. Hence, not much work is required in estimating the
MLE for these two distributions, since the MLEs has nice analytic closed
forms, which can be used as start values.
}
