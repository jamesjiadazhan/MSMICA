#' Select the most likely feature for a single metabolite
#'
#' Given a vector of log-scale Bayesian scores, one per candidate
#' feature, this helper picks the feature with the highest posterior
#' probability and also returns the normalized posterior probability
#' across all candidates. Non-finite scores are treated as impossible
#' and excluded from the normalization. If no candidate has a finite
#' score (or the input is empty), the function returns \code{NULL} so
#' that the caller can skip that metabolite. Used by
#' \code{local_optimization_per_metabolite()}.
#'
#' @param feature_log_score Numeric vector of log-posterior scores, one
#'   per candidate feature.
#' @return Either \code{NULL} (no viable candidate) or a list with
#'   elements \code{best_feature_index}, \code{best_probability},
#'   \code{best_log_posterior}, and \code{all_probabilities} (a vector of
#'   posterior probabilities over all candidates, with \code{NA} at
#'   non-finite positions).
#' @keywords internal
#' @noRd
exact_bayesian_feature_selection_for_metabolite = function(feature_log_score) {
    if (length(feature_log_score) == 0) {
        return(NULL)
    }

    # Drop features with non-finite scores (e.g. incompatible mz/RT or missing data)
    finite_idx = is.finite(feature_log_score)

    if (!any(finite_idx)) {
        return(NULL)
    }

    candidate_idx = which(finite_idx)
    candidate_log_scores = feature_log_score[finite_idx]

    # Softmax normalization over candidate features in log-space
    normalized_log_probability = candidate_log_scores - logsumexp(candidate_log_scores)
    probability_values = exp(normalized_log_probability)

    # Pick the MAP feature (highest posterior)
    best_idx_local = which.max(candidate_log_scores)
    best_feature_idx = candidate_idx[best_idx_local]

    # Place probabilities back into a full-length vector so callers can index it
    # against the original feature list (non-candidates remain NA).
    all_probabilities = rep(NA_real_, length(feature_log_score))
    all_probabilities[candidate_idx] = probability_values

    list(
        best_feature_index = best_feature_idx,
        best_probability = probability_values[best_idx_local],
        best_log_posterior = candidate_log_scores[best_idx_local],
        all_probabilities = all_probabilities
    )

}