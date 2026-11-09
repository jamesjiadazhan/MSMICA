#' Build empirical evidence calibration for MSMICA likelihood terms
#'
#' The calibration object stores high-confidence training distributions
#' learned inside a dataset run. It is intentionally simple and
#' deterministic: RT evidence uses the empirical survival function of
#' absolute RT errors, while correlation evidence uses the empirical CDF
#' of Fisher-z transformed absolute correlations.
#'
#' @keywords internal
#' @noRd
build_msmica_empirical_calibration = function(rt_errors = NULL,
                                              corr_z = NULL,
                                              floor = 1e-6) {
    list(
        rt_errors = sort(abs(rt_errors[is.finite(rt_errors)])),
        corr_z = sort(corr_z[is.finite(corr_z)]),
        floor = floor
    )
}

#' Empirical RT log-likelihood from calibration anchors
#' @keywords internal
#' @noRd
empirical_rt_log_likelihood = function(time_difference, calibration) {
    if (is.null(calibration) || length(calibration$rt_errors) == 0 || !is.finite(time_difference)) {
        return(NA_real_)
    }
    p = mean(calibration$rt_errors >= abs(time_difference))
    log(max(p, calibration$floor))
}

#' Retention-time log likelihood with a safe no-calibration behavior
#' @keywords internal
#' @noRd
msmica_rt_log_likelihood = function(time_difference, calibration, rt_sigma) {
    log_rt = empirical_rt_log_likelihood(time_difference, calibration)
    if (is.finite(log_rt)) {
        return(log_rt)
    }
    if (length(rt_sigma) != 1 || !is.finite(rt_sigma) || rt_sigma <= 0 ||
        !is.finite(time_difference)) {
        # Missing calibration is not negative evidence and must not be replaced
        # by an invented uncertainty. Omit the RT term from the score.
        return(0)
    }
    log(max(dnorm(time_difference, mean = 0, sd = rt_sigma), 1e-300))
}

#' Correlation log likelihood with a safe no-calibration behavior
#' @keywords internal
#' @noRd
msmica_corr_log_likelihood = function(correlation, calibration, corr_mu, corr_sigma) {
    log_corr = empirical_corr_log_likelihood(correlation, calibration)
    if (is.finite(log_corr)) {
        return(log_corr)
    }
    if (!is.finite(correlation) || length(corr_mu) != 1 ||
        !is.finite(corr_mu) || length(corr_sigma) != 1 ||
        !is.finite(corr_sigma) || corr_sigma <= 0) {
        # No valid correlation calibration/prior is available. Omit the
        # correlation term instead of manufacturing a training distribution.
        return(0)
    }
    corr_value = min(max(correlation, -0.999999), 0.999999)
    z_obs = atanh(corr_value)
    log(max(dnorm(z_obs, mean = corr_mu, sd = corr_sigma), 1e-300))
}

#' Empirical correlation log-likelihood from calibration anchors
#' @keywords internal
#' @noRd
empirical_corr_log_likelihood = function(correlation, calibration) {
    if (is.null(calibration) || length(calibration$corr_z) == 0 || !is.finite(correlation)) {
        return(NA_real_)
    }
    corr_value = min(max(abs(correlation), 0), 0.999999)
    z_obs = atanh(corr_value)
    p = mean(calibration$corr_z <= z_obs)
    log(max(p, calibration$floor))
}

#' Write calibration summary and distributions when detail output is enabled
#' @keywords internal
#' @noRd
write_msmica_calibration_outputs = function(calibration, output_folder) {
    if (is.null(calibration) || is.null(output_folder)) {
        return(invisible(NULL))
    }

    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

    summary_df = data.frame(
        evidence = c("rt_error", "correlation_z"),
        n = c(length(calibration$rt_errors), length(calibration$corr_z)),
        mean = c(
            mean(calibration$rt_errors, na.rm = TRUE),
            mean(calibration$corr_z, na.rm = TRUE)
        ),
        sd = c(
            sd(calibration$rt_errors, na.rm = TRUE),
            sd(calibration$corr_z, na.rm = TRUE)
        )
    )

    readr::write_csv(summary_df, file.path(output_folder, "empirical_evidence_calibration_summary.csv"))

    max_len = max(length(calibration$rt_errors), length(calibration$corr_z), 1)
    distribution_df = data.frame(
        rt_error = c(calibration$rt_errors, rep(NA_real_, max_len - length(calibration$rt_errors))),
        correlation_z = c(calibration$corr_z, rep(NA_real_, max_len - length(calibration$corr_z)))
    )

    readr::write_csv(distribution_df, file.path(output_folder, "empirical_evidence_calibration_distributions.csv"))

    invisible(NULL)
}
