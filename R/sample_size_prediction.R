#' Sample Size Data
#'
#' Create a data frame of sample counts and cumultative counts over time, which
#' be used for aggregating raw data for sample size plots using
#' \code{SampleSizePlot}.
#'
#' @param data A data frame containing the unit variable (e.g., user id),
#'   a grouping variable, and a date variable.
#' @param unit String name of unit table.
#' @param x String name of the variable for the x-axis of the plot, either
#'   a character or Date vector for dates, or a numeric vector, e.g., for
#'   number of days.
#' @param group String name of grouping variable.
#' @param compute_days Whether to compute the number of days from a reference
#'   date if the \code{x} argument is a date or character vector.
#' @param ref_date The reference date for computing the number of days. If
#'   \code{NULL}, then the minimum date in the \code{data} data frame is used.
#' @return A data frame containing the x-axis variable(s), grouping variables,
#'   and the sample size counts and cumulative counts for each combination of
#'   x variable(s) and grouping variables.
#' @export

SampleSizeData <- function(data, unit = "obscured_gaia_id",
                           group = "experiment_id",
                           x = "min_exposure_date", compute_days = TRUE,
                           ref_date = NULL) {
  # Check for x variable type and convert character vectors to Date vectors.
  stopifnot(
    is.numeric(data[[x]]) || is.character(data[[x]]) || is.Date(data[[x]])
  )
  if (is.character(data[[x]])) {
    data[[x]] <- as.Date(data[[x]])
  }
  data$unit <- data[[unit]]
  group_vars <- c(group, x)
  # Aggregate by unit and group variables, in case there are duplicate units,
  # before aggregating by the group variables, to compute aggregated counts
  # and cumulative counts.
  out <- data %>%
    group_by_at(c("unit", group_vars)) %>%
    group_by_at(group_vars) %>%
    summarize(count = length(unique(unit))) %>%
    mutate(cum_count = cumsum(count))
  # Compute the number of days from a reference date, if the x variable is a
  # date and compute_days = TRUE.
  if (is.Date(out[[x]]) && compute_days) {
    if (is.null(ref_date)) {
      ref_date <- min(out[[x]])
    }
    out$days <- as.numeric(out[[x]] - ref_date) + 1
    out <- out %>%
      dplyr::select_at(c(group, x, "days", "count", "cum_count"))
  }
  return(out)
}

#' Extend Data
#'
#' Extend a data frame containing time variables, grouping variables, and
#' other variables over a longer time horizon, so that it can be used for
#' forecasting, e.g., forecasting the growth of sample size counts.
#'
#' @param data A data frame containing the variables to extend and, optionally,
#'   the grouping variables and other variables.
#' @param extend_vars A character vector containing names of variables to extend
#'   over the specified time horizon, typically, the time variables used in the
#'   forecast exercise, e.g., the number of days, date, etc.
#' @param group_vars A character vector containing names of the grouping
#'   variables. These can be considered the keys for the data that we group by.
#' @param other_vars A character vector containing names of other measure
#'   variables to include in the data frame returned. These will have \code{NA}
#'   values for time points in the extension.
#' @param horizon A numeric vector of positive integers to add to the time
#'   variables to extend. Defaults to \code{1:60}.
#' @return The original data frame passed into \code{data} but with additional
#'   rows added for the forecast horizon.
#' @export

ExtendData <- function(data, extend_vars = "days", group_vars = NULL,
                       other_vars = NULL, horizon = 1:60) {
  stopifnot(length(extend_vars) > 0 & length(horizon) > 0)
  h <- length(horizon)
  # Create data frame of existing data with
  # columns appropriately selected and arranged
  existing_data <- data %>%
    ungroup() %>%
    select_at(c(group_vars, extend_vars, other_vars))
  # Create new data
  if (is.null(group_vars)) {
    new_data <- list()
    for (x in extend_vars) {
      new_data[[x]] <- max(data[[x]]) + horizon
    }
    # The for loop won't execute if other_vars is NULL.
    for (x in other_vars) {
      new_data[[x]] <- NA
    }
    new_data <- data.frame(new_data)
  } else {
    data_list <- split(data, data[, group_vars], drop = TRUE)
    new_data <- list()
    for (i in seq_along(data_list)) {
      lhs <- data_list[[i]] %>%
        select_at(group_vars) %>%
        slice(1)
      # Recursively calls ExtendData for keys = NULL
      data_selected <- data_list[[i]] %>%
        ungroup() %>%
        select_at(c(extend_vars, other_vars))
      rhs <- ExtendData(
        data_selected,
        extend_vars = extend_vars,
        other_vars = other_vars,
        horizon = horizon
      )
      new_data[[i]] <- bind_cols(lhs, rhs)
    }
    new_data <- do.call(bind_rows, new_data)
  }
  # Create output data from existing and new data
  out <- existing_data %>%
    bind_rows(new_data)
  return(out)
}

#' Sample Size Model
#'
#' Fits a forecast model for a sample size time series using \code{mgcv::gam}
#' for a spline regression model or \code{stats::lm} for a linear regression
#' model. Note the following:
#' (1) The spline regression model is recommended for forecasting sample size
#'     growth, while the linear regression model option is added as a benchmark.
#' (2) While the name of the function is \code{SampleSizeModel}, it can be used
#'     for forecasting any numeric variable \code{y} ordered by a numeric
#'     variable \code{x}.
#'
#' @param data A data frame containing the unit variable (e.g., user id),
#'   a grouping variable, and a date variable.
#' @param y String name for the y-axis variable, which has to be a numeric
#'   vector.
#' @param x String name for the x-axis variable of the plot, which has to be a
#'   numeric vector, e.g., the number of days or a time index.
#' @param group String name of grouping variable. Note the function is currently
#'   limited to using a single grouping variable. Hence, only the first entry
#'   is used if a character vector is passed in as the argument.
#' @param method Method used to fit a forecast model, either a spline regression
#'   model using \code{mgcv::gam} or a linear regression model using
#'   \code{stats::lm}.
#' @param k The number of knots for the spline regression, an argument passed
#'   into the {mgcv::s} smoothing function when \code{method = "spline"}.
#' @param bs The type of basis function for the spline regression, argument
#'   into the {mgcv::s} smoothing function when \code{method = "spline"}.
#'   Available types include \code{"tp"} for thin-plate splines and \code{"cr"}
#'   for cubic regression splines. See the following \code{mgcv} documentation
#'   page for other options:
#'   https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.terms.html
#' @return A \code{mgcv::gam} object or a \code{stats::lm} object for the fitted
#'   regression model for forecasting sample size growth.
#' @export

SampleSizeModel <- function(data, y, x, group = NULL,
                            method = c("spline", "linear"), k = -1, bs = "tp") {
  method <- match.arg(method)
  if (nzchar(group)) {
    group <- group[1]
    data[[group]] <- factor(data[[group]])
    formula_rhs <- switch(method,
      spline = glue("{group} + s({x}, by = {group}, k = {k}, bs = '{bs}')"),
      linear = glue("{group} + {x}:{group}")
    )
  } else {
    formula_rhs <- switch(method,
      spline = glue("s({x}, k = {k}, bs = '{bs}')"),
      linear = glue("{x}")
    )
  }
  model_formula <- Formula(y, formula_rhs)
  FUN <- switch(method,
    spline = mgcv::gam,
    linear = stats::lm
  )
  model <- FUN(model_formula, data = data)
  return(model)
}

#' Sample Size Prediction
#'
#' Add a prediction of sample size counts to a data frame containing sample
#' size data. Extend the data frame over a forecast horizon using
#' \code{ExtendData} as needed.
#'
#' @param data A data frame containing the unit variable (e.g., user id),
#'   a grouping variable, and a date variable.
#' @param y String name for the y-axis variable, which has to be a numeric
#'   vector.
#' @param y_pred String name for the prediction of the y-axis variable. Defaults
#'   to \code{y} followed by a suffix \code{"_pred"}.
#' @param x String name for the x-axis variable of the plot, which has to be a
#'   numeric vector, e.g., the number of days or a time index.
#' @param date String name for the date variable corresponding to the x-axis
#'   variable.
#' @param group_vars A character vector of names of grouping variables. If there
#'   are multiple grouping variables, a single grouping variable will be created
#'   as a concatenation of all the grouping variable for the
#'   \code{SampleSizeModel} function call.
#' @param other_vars A character vector of names of measure variables to retain
#'   in the prediction data. These names will be passed into the
#'   \code{other_vars} argument of the \code{ExtendData} function call.
#' @param train A numeric vector for values of the \code{x} variable to train
#'   on. Defaults to \code{NULL}, where all values of \code{x} are used. Note
#'   in addition to this filter, all observations with \code{NA} values for
#'   \code{y} will be removed as well.
#' @param horizon A numeric vector of positive integers to add to the time
#'   variables to extend. Defaults to \code{1:60}.
#' @param ... Additional arguments to pass into the \code{SampleSizeModel}
#'   function.
#' @return A data frame containing the variables specified in \code{y},
#'   \code{x}, \code{group_vars}, and \code{other_vars}, and the prediction
#'   variable \code{y_pred} and optionally, the combined group variable
#'   \code{"group"} if there are multiple grouping variables.
#' @export

SampleSizePrediction <- function(data, y = "cum_count",
                                 y_pred = paste0(y, "_pred"), x = "days",
                                 date = "min_exposure_date",
                                 group_vars = NULL, other_vars = "count",
                                 train = NULL, horizon = 1:60, ...) {
  # Define single grouping variable from the group_vars vector
  if (length(group_vars) == 0) {
    group <- NULL
  } else if (length(group_vars) == 1) {
    group <- group_vars
  } else {
    group <- "group"
    group_vars <- c(group_vars, group)
    data[[group]] <- data[[group_vars[1]]]
    for (i in seq_along(length(group_vars) - 1)) {
      data[[group]] <- paste(data[[group]], data[[group_vars[1 + i]]])
    }
  }
  # Train the model on non-missing values of y and the
  model_data_filter <- !is.na(data[[y]])
  if (!is.null(train)) {
    model_data_filter <- model_data_filter & (data[[x]] %in% train)
  }
  model_data <- data[model_data_filter, ]
  pred_model <- SampleSizeModel(model_data, y = y, x = x, group = group, ...)
  # Extend the data if the forecast horizon is specified
  if (is.null(horizon)) {
    pred_data <- data
  } else {
    pred_data <- ExtendData(
      data,
      extend_vars = c(date, x), other_vars = c(other_vars, y),
      group_vars = group_vars, horizon = horizon
    )
  }
  # Add prediction to data
  pred_data[[y_pred]] <- pred_model %>%
    predict(newdata = pred_data) %>%
    as.numeric()
  return(pred_data)
}

#' Sample Size Plot
#'
#' Create a plot of sample size counts over time using \code{ggplot2::ggplot},
#' optionally using \code{SampleSizeData} to create aggregated data for the
#' plot. We have the following dimensions for defining the plot:
#' (1) Data Aggregation: The function can aggregate the data using
#'     \code{SampleSizeData} if \code{preprocess == TRUE}. Otherwise, the
#'     function will take the data as is and assume that the data is already
#'     aggregated.
#' (2) Sample Size Counts: If the data needs to be aggregated first, then
#'     \code{y} would be set to one of two values, i.e., \code{y == "cum_count"}
#'     if \code{cumulative == TRUE} and \code{y == "count"} if
#'     if \code{cumulative == FALSE}. If the data is taken as is, then \code{y}
#'     should ideally be specified. If it is not, then \code{y == "cum_count"}
#'     or \code{y == "count"} depending on the \code{cumulative} argument.
#' (3) Estimated Mean: In addition to the line plot of the sample size counts,
#'     this function can include the line plot of estimated means of the counts,
#'     either provided in the data using the name provided in \code{y_pred}, or
#'     by adding a \code{geom_smooth()} to the plot with \code{method = "gam"}
#'     for a generalized additive model (GAM) smooth fit using \code{mgcv::gam}.
#'     The advantage of providing the estimated mean using \code{y_pred} is that
#'     we can plot \code{y_pred} for time points with \code{NA} values of
#'     \code{y}.
#'   See http://go/cloud-sample-size-prediction for details on sample size
#' growth.
#'
#' @param data A data frame containing the unit variable (e.g., user id),
#'   a grouping variable, and a date variable.
#' @param unit String name of unit variable, used when the data needs to be
#'   aggregated first by specifying \code{preprocess = TRUE}.
#' @param y String name of the y-axis variable to plot.
#' @param y_pred String name of the prediction of the y-axis variable to plot.
#' @param x String name of the variable for the x-axis of the plot, either
#'   a character or Date vector for dates, or a numeric vector, e.g., for
#'   number of days.
#' @param facet String name of the variable to facet the plot by.
#' @param color String name of the variable to color the plot by.
#' @param pred_range The range to show the prediction \code{y_pred}. If
#'   \code{"horizon"}, \code{y_pred} will be shown for values of \code{x} where
#'   the variable \code{y} has \code{NA} values. If \code{"all"}, \code{y_pred}
#'   will be shown for all values of \code{x}.
#' @param preprocess Whether to aggregate the data by date and group before
#'   passing the data to \code{ggplot2::ggplot}.
#' @param cumulative Whether the sample size counts are cumulative counts over
#'   time. Defaults to \code{TRUE}.
#' @param smooth Whether to add a smoothed mean of the data in the plot, in
#'   addition to the line geom of the \code{y} variable.
#' @param facet_ncol The number of columns for plot faceting.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @return A line plot for the sample size counts, whether individual time
#'   counts or cumulative counts over time.
#' @export

SampleSizePlot <- function(data, unit = "obscured_gaia_id", y = NULL,
                           y_pred = NULL, x = "min_exposure_date",
                           facet = NULL, color = NULL,
                           pred_range = c("horizon", "all"), preprocess = FALSE,
                           cumulative = TRUE, smooth = FALSE, facet_ncol = 2,
                           facet_scales = "free", xlab = NULL, ylab = NULL) {
  # Check for x variable type and convert character vectors to Date vectors.
  stopifnot(
    is.numeric(data[[x]]) || is.character(data[[x]]) || is.Date(data[[x]])
  )
  if (is.character(data[[x]])) {
    data[[x]] <- as.Date(data[[x]])
  }
  pred_range <- match.arg(pred_range)
  # Aggregate data if preprocess argument is TRUE.
  if (preprocess) {
    group_vars <- unique(c(facet, color))
    plot_data <- SampleSizeData(
      data = data, unit = unit, x = x, group = group_vars
    )
  } else {
    plot_data <- data
  }
  # If specified, show prediction only for NA values of the response predicted
  if (!is.null(y_pred) && pred_range == "horizon") {
    plot_data[[y_pred]] <- ifelse(
      is.na(plot_data[[y]]), plot_data[[y_pred]], NA
    )
  }
  # Specify y and ylab variables if not specified.
  if (is.null(y)) {
    y <- ifelse(cumulative, "cum_count", "count")
  }
  if (is.null(ylab)) {
    ylab <- paste0(ifelse(cumulative, "Cumulative Count", "Count"), "\n")
  }
  # Create plot of sample size
  p <- ggplot(plot_data, aes_string(x = x, color = color)) +
    geom_line(aes_string(y = y)) +
    xlab(xlab) +
    ylab(ylab)
  # Add faceting if specified
  if (!is.null(facet)) {
    facet_formula <- as.formula(paste("~", facet))
    p <- p +
      facet_wrap(facet_formula, ncol = facet_ncol, scales = facet_scales)
  }
  # Add plot of sample size prediction if specified
  if (!is.null(y_pred)) {
    p <- p + geom_line(aes_string(y = y_pred), color = "blue")
  }
  # Add smoothed mean of sample size if specified
  if (smooth) {
    p <- p + geom_smooth(aes_string(y = y), formula = y ~ s(x), method = "gam")
  }
  return(p)
}
