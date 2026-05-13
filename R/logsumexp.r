#' Numerically stable log-sum-exp
#'
#' Computes \code{log(sum(exp(log_values)))} in a numerically stable way by
#' factoring out the maximum value before exponentiation. Non-finite values
#' (\code{NA}, \code{NaN}, \code{-Inf}, \code{Inf}) are dropped prior to the
#' calculation. If no finite values remain, \code{-Inf} is returned. This
#' helper underlies the Bayesian normalization used by the local optimization
#' routines.
#'
#' @param log_values Numeric vector of values already on the log scale.
#' @return A single numeric value equal to \code{log(sum(exp(log_values)))},
#'   or \code{-Inf} when the input contains no finite elements.
#' @keywords internal
#' @noRd
logsumexp = function(log_values) {
    # Drop NA / NaN / +-Inf so the shift-by-max trick does not propagate them
    finite_log_values = log_values[is.finite(log_values)]

    if (length(finite_log_values) == 0) {
        return(-Inf)
    }

    # Shift by max to keep exp(.) in a safe range before summing
    max_log_value = max(finite_log_values)
    max_log_value + log(sum(exp(finite_log_values - max_log_value)))
}