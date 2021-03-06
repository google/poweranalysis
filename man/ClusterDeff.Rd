% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_analysis.R
\name{ClusterDeff}
\alias{ClusterDeff}
\title{Cluster Design Effect}
\usage{
ClusterDeff(
  cluster_size = NULL,
  mean = 1,
  sd = 0,
  k = NULL,
  icc = 0,
  type = c("min", "max")
)
}
\arguments{
\item{cluster_size}{A numeric vector of cluster sizes for the study.}

\item{mean}{The mean of the cluster sizes. Used if
\code{cluster_size} is \code{NULL}.}

\item{sd}{The standard deviation (SD) of the cluster sizes, mainly for
computing the CV, i.e., the SD divided by the mean. Used if
\code{cluster_size} is \code{NULL}.}

\item{k}{The number of clusters.}

\item{icc}{The intra-cluster correlation (ICC)}

\item{type}{Either the minimum design effect (\code{"min"}) or the maximum
(\code{"max"}). TODO(shyuemingloh@): Implement other design effect types
when the actual cluster sizes are provided.}
}
\value{
The cluster design effect for the study.
}
\description{
Estimate the design effect in clustered data from having within-cluster
correlation relative to a zero-correlation design, with the following two
options:
(1) Minimum design effect, given by \code{1 + (mean - 1) * icc}, where
    \code{mean} is the mean of the cluster sizes and \code{icc}, the
    intra-cluster correlation (ICC), i.e., the within-cluster correlation.
(2) Maximum design effect, given by
    \code{1 + ((cv^2 * (1 - 1 / k) + 1) * mean - 1) * icc}, where
    \code{cv} is the coefficient of variation (CV) for the cluster sizes and
    \code{k}, the number of clusters. If \code{k} is \code{NULL}, we assume
    that \code{k} is \code{Inf}, so that the design effect becomes
    \code{1 + (cv^2 * mean - 1) * icc}.
Eldridge et al. (2006) provided the above lower and upper bounds of the
cluster design effect, because the actual cluster sizes are typically not
known before study begins.
}
\references{
\itemize{
Eldridge, Sandra M., Deborah Ashby, and Sally Kerry (2006). "Sample Size for
  Cluster Randomized Trials: Effect of Coefficient of Variation of Cluster
  Size and Analysis Method." \emph{International Journal of Epidemiology,}
  35(5): 1292-1300.
}
}
