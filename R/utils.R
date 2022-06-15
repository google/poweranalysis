#' Create Formula Object from String Names
#'
#' Creates a formula object from string names of variables on left-hand side
#' (LHS) and right-hand side (RHS) of the formula.
#'
#' @param y A character vector containing names of the variables on the LHS of
#'   the formula.
#' @param x A character vector containing names of the variables on the RHS of
#'   the formula.
#' @param collapse A character for separating variables on either side of the
#'   the formula. Defaults to \code{+}.
#' @return A formula object with the LHS and RHS variables specified.

CreateFormula <- function(y = NULL, x = NULL, collapse = "+") {
  if (is.null(y) && is.null(x)) {
    stop("At least one of y and x needs to be non-null.")
  }
  if (is.null(y)) {
    y <- ""
  }
  y_formula <- paste(y, collapse = collapse)
  if (is.null(x)) {
    out <- as.formula(y_formula)
  } else {
    x_formula <- paste(x, collapse = collapse)
    out <- as.formula(paste(y_formula, x_formula, sep = "~"))
  }
  return(out)
}

#' Apply Percent Printing Format to Columns in a Data Frame
#'
#' Apply percent printing format from \code{formattable::percent} to a selection
#' of columns in a data frame.
#'
#' @param data A data frame.
#' @param x A character vector containing names of columns to format.
#' @param ... Additional arguments to pass to \code{formattable::percent}.
#' @return The data frame passed into the function with the
#'  \code{formattable::percent} printing format applied to the selected columns.

Percentize <- function(data, x = NULL, ...) {
  if (!is.null(x)) {
    data[x] <- data[x] %>%
      lapply(formattable::percent, ...)
  }
  return(data)
}
