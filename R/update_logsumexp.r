#' Incrementally update a log-sum-exp accumulator
#'
#' Adds a single new log-scale value to a running \code{log(sum(exp(...)))}
#' accumulator without having to recompute it from scratch. Non-finite
#' inputs are handled gracefully: a non-finite \code{new_log_value} leaves
#' the accumulator unchanged, and a non-finite \code{current_log_sum} is
#' replaced by \code{new_log_value}. Used by the iterative local
#' optimization routines to keep partial normalizers.
#'
#' @param current_log_sum Numeric. The current running log-sum-exp value.
#' @param new_log_value Numeric. A new log-scale value to fold in.
#' @return A numeric value equal to \code{log(exp(current_log_sum) + exp(new_log_value))}.
#' @keywords internal
#' @noRd
update_logsumexp = function(current_log_sum, new_log_value) {
    # If the new value is non-finite, nothing to add: keep the accumulator
    if (!is.finite(new_log_value)) {
        return(current_log_sum)
    }

    # If the accumulator has never been seeded, initialize it with the new value
    if (!is.finite(current_log_sum)) {
        return(new_log_value)
    }

    # Standard shift-by-max trick for numerically stable addition in log-space
    max_log_value = max(current_log_sum, new_log_value)
    max_log_value + log(exp(current_log_sum - max_log_value) + exp(new_log_value - max_log_value))
}