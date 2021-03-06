% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{Percentize}
\alias{Percentize}
\title{Apply Percent Printing Format to Columns in a Data Frame}
\usage{
Percentize(data, x = NULL, ...)
}
\arguments{
\item{data}{A data frame.}

\item{x}{A character vector containing names of columns to format.}

\item{...}{Additional arguments to pass to \code{formattable::percent}.}
}
\value{
The data frame passed into the function with the
 \code{formattable::percent} printing format applied to the selected columns.
}
\description{
Apply percent printing format from \code{formattable::percent} to a selection
of columns in a data frame.
}
