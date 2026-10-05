# Apply empirical clustering thresholds when estimable, otherwise return
# explicit fixed defaults and preserve the calibration failure reason.
msmica_threshold_estimate_with_fallback = function(estimator,
                                                   estimator_args,
                                                   type = c("adduct", "isotope"),
                                                   default_time,
                                                   default_correlation,
                                                   capture_fraction = 0.8,
                                                   correlation_floor = 0.4) {
    type = match.arg(type)
    time_name = if (type == "adduct") {
        "adduct_correlation_time_threshold"
    } else {
        "isotopic_correlation_time_threshold"
    }
    correlation_name = if (type == "adduct") {
        "adduct_correlation_r_threshold"
    } else {
        "isotopic_correlation_r_threshold"
    }
    result = tryCatch(
        do.call(estimator, estimator_args),
        error = function(e) structure(list(error = conditionMessage(e)), class = "msmica_threshold_error")
    )

    failure_reason = if (inherits(result, "msmica_threshold_error")) {
        result$error
    } else if (
        is.null(result[[time_name]]) || length(result[[time_name]]) != 1 ||
        !is.finite(result[[time_name]]) ||
        is.null(result[[correlation_name]]) || length(result[[correlation_name]]) != 1 ||
        !is.finite(result[[correlation_name]])
    ) {
        "empirical estimator returned a missing or non-finite threshold"
    } else {
        NULL
    }

    if (is.null(failure_reason)) {
        return(result)
    }

    message(
        "Empirical ", type,
        " threshold estimation unavailable; using fixed defaults (RT <= ",
        default_time, " sec, Spearman r >= ", default_correlation,
        "): ", failure_reason
    )

    pair_match = regexec("only ([0-9]+) calibration pairs", failure_reason)
    pair_capture = regmatches(failure_reason, pair_match)[[1]]
    observed_pairs = if (length(pair_capture) >= 2) as.integer(pair_capture[2]) else 0L

    summary = data.frame(
        threshold_method = "fixed_default_fallback",
        fallback_reason = failure_reason,
        capture_fraction = capture_fraction,
        n_pairs = observed_pairs,
        n_pairs_for_threshold = observed_pairs,
        adduct_correlation_floor = if (type == "adduct") correlation_floor else NA_real_,
        isotopic_correlation_floor = if (type == "isotope") correlation_floor else NA_real_,
        adduct_correlation_time_threshold = if (type == "adduct") default_time else NA_real_,
        adduct_correlation_r_threshold = if (type == "adduct") default_correlation else NA_real_,
        isotopic_correlation_time_threshold = if (type == "isotope") default_time else NA_real_,
        isotopic_correlation_r_threshold = if (type == "isotope") default_correlation else NA_real_,
        observed_time_capture = NA_real_,
        observed_correlation_capture = NA_real_,
        median_time_difference = NA_real_,
        median_correlation = NA_real_,
        stringsAsFactors = FALSE
    )
    result = list(
        pairs = data.frame(),
        summary = summary
    )
    result[[time_name]] = default_time
    result[[correlation_name]] = default_correlation
    result
}
