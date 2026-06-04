#' Greedy metabolite-first Bayesian feature selection
#'
#' Assigns at most one feature to each metabolite, processing metabolites
#' from highest to lowest concentration prior. Once a feature is selected
#' for one metabolite, it is removed from the candidate pool for later
#' metabolites. The direct concentration prior is not added here; the
#' concentration-to-intensity relationship enters only through the same
#' rank-matching bonus used by \code{exact_bayesian_assignment()}.
#'
#' @keywords internal
#' @noRd
exact_bayesian_metabolite_feature_assignment = function(local_log_score_matrix,
                                                        metabolite_log_concentration,
                                                        feature_log_intensity,
                                                        feature_rank_prior_strength = 1) {
    n_features = nrow(local_log_score_matrix)
    n_metabolites = ncol(local_log_score_matrix)

    if (n_features == 0 || n_metabolites == 0) {
        return(NULL)
    }

    feature_intensity_scaled =
        (feature_log_intensity - min(feature_log_intensity)) /
        max(1e-12, max(feature_log_intensity) - min(feature_log_intensity))

    metabolite_abundance_scaled =
        (metabolite_log_concentration - min(metabolite_log_concentration)) /
        max(1e-12, max(metabolite_log_concentration) - min(metabolite_log_concentration))

    metabolite_order = order(metabolite_log_concentration, decreasing = TRUE, na.last = TRUE)
    available_features = rep(TRUE, n_features)

    selected_feature = integer(0)
    selected_metabolite = integer(0)
    selected_probability = numeric(0)
    selected_log_posterior = numeric(0)

    for (metabolite_idx in metabolite_order) {
        row_idx = which(available_features & is.finite(local_log_score_matrix[, metabolite_idx]))

        if (length(row_idx) == 0) {
            next
        }

        candidate_log_scores = local_log_score_matrix[row_idx, metabolite_idx]
        candidate_rank_bonus = -feature_rank_prior_strength *
            abs(feature_intensity_scaled[row_idx] - metabolite_abundance_scaled[metabolite_idx])
        adjusted_log_scores = candidate_log_scores + candidate_rank_bonus

        normalized_log_probability = adjusted_log_scores - logsumexp(adjusted_log_scores)
        probability_values = exp(normalized_log_probability)

        best_idx_local = which.max(adjusted_log_scores)
        best_feature_idx = row_idx[best_idx_local]

        selected_feature = c(selected_feature, best_feature_idx)
        selected_metabolite = c(selected_metabolite, metabolite_idx)
        selected_probability = c(selected_probability, probability_values[best_idx_local])
        selected_log_posterior = c(selected_log_posterior, adjusted_log_scores[best_idx_local])

        available_features[best_feature_idx] = FALSE

        if (!any(available_features)) {
            break
        }
    }

    if (length(selected_feature) == 0) {
        return(NULL)
    }

    list(
        feature_index = selected_feature,
        metabolite_index = selected_metabolite,
        probability = selected_probability,
        log_posterior = selected_log_posterior
    )
}