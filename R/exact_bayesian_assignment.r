#' Feature-wise Bayesian assignment of metabolites
#'
#' Given a pre-computed log-posterior score matrix where rows are features
#' and columns are candidate metabolites, this helper assigns each feature
#' independently to the metabolite that maximises its posterior score. A
#' feature-rank prior links intense features to abundant metabolites: the
#' stronger the feature (by log mean intensity), the closer the expected
#' log concentration of its assigned metabolite. This gives MSMICA a soft
#' "rank-matching" preference without forcing a hard one-to-one assignment.
#' Used internally by \code{local_optimization_per_feature()}.
#'
#' @param local_log_score_matrix Numeric matrix of log-posterior scores
#'   with features in rows and metabolites in columns.
#' @param feature_rank_order Integer vector giving the order in which
#'   features should be scored (typically decreasing log intensity).
#' @param metabolite_log_concentration Numeric vector of log concentration
#'   priors for the candidate metabolites (column order).
#' @param feature_log_intensity Numeric vector of log mean intensities for
#'   the features (row order).
#' @param feature_rank_prior_strength Numeric scalar; how strongly the
#'   rank-matching bonus pulls strong features toward abundant
#'   metabolites. Default is 1.
#' @return \code{NULL} if the matrix has no finite rows/columns;
#'   otherwise a list with \code{assignment} (metabolite index per
#'   feature in original order), \code{probability} (posterior
#'   probability of the chosen metabolite), and \code{log_posterior}
#'   (score of the chosen pair).
#' @keywords internal
#' @noRd
exact_bayesian_assignment = function(local_log_score_matrix,
                                    feature_rank_order,
                                    metabolite_log_concentration,
                                    feature_log_intensity,
                                    feature_rank_prior_strength = 1) {
    # Process features in the caller-supplied order so the rank bonus can
    # reference a "strongest feature first" ordering.
    ranked_local_log_score_matrix = local_log_score_matrix[feature_rank_order, , drop = FALSE]
    ranked_feature_log_intensity = feature_log_intensity[feature_rank_order]

    n_features = nrow(ranked_local_log_score_matrix)
    n_metabolites = ncol(ranked_local_log_score_matrix)

    if (n_features == 0 || n_metabolites == 0) {
        return(NULL)
    }

    best_assignment_ranked = integer(n_features)
    pair_probability_ranked = numeric(n_features)
    pair_log_posterior_ranked = numeric(n_features)

    # Rescale feature intensities and metabolite concentrations to a common
    # [0, 1] interval so the two are directly comparable for the rank bonus.
    # The 1e-12 floor guards against zero-variance inputs.
    feature_intensity_scaled =
        (ranked_feature_log_intensity - min(ranked_feature_log_intensity)) /
        max(1e-12, max(ranked_feature_log_intensity) - min(ranked_feature_log_intensity))

    metabolite_abundance_scaled =
        (metabolite_log_concentration - min(metabolite_log_concentration)) /
        max(1e-12, max(metabolite_log_concentration) - min(metabolite_log_concentration))

    for (rank_idx in seq_len(n_features)) {
        row_vals = ranked_local_log_score_matrix[rank_idx, ]
        finite_idx = is.finite(row_vals)

        if (!any(finite_idx)) {
            return(NULL)
        }

        candidate_cols = which(finite_idx)
        candidate_log_scores = row_vals[finite_idx]

        # Rank bonus: reward metabolites whose (scaled) log concentration is
        # close to this feature's (scaled) log intensity. The larger the
        # absolute gap, the more negative the bonus, so the score shrinks.
        candidate_rank_bonus = -feature_rank_prior_strength *
            abs(metabolite_abundance_scaled[candidate_cols] - feature_intensity_scaled[rank_idx])

        adjusted_log_scores = candidate_log_scores + candidate_rank_bonus

        # Softmax over candidate metabolites for this feature (log-stable)
        normalized_log_probability = adjusted_log_scores - logsumexp(adjusted_log_scores)
        probability_values = exp(normalized_log_probability)

        # MAP assignment for this feature
        best_idx_local = which.max(adjusted_log_scores)
        best_metabolite_idx = candidate_cols[best_idx_local]

        best_assignment_ranked[rank_idx] = best_metabolite_idx
        pair_probability_ranked[rank_idx] = probability_values[best_idx_local]
        pair_log_posterior_ranked[rank_idx] = adjusted_log_scores[best_idx_local]
    }

    # Re-map results back into the caller's original (unranked) feature order
    # so downstream joins by feature position remain valid.
    best_assignment_original_order = integer(n_features)
    pair_probability_original_order = numeric(n_features)
    pair_log_posterior_original_order = numeric(n_features)

    for (rank_idx in seq_len(n_features)) {
        original_feature_idx = feature_rank_order[rank_idx]

        best_assignment_original_order[original_feature_idx] = best_assignment_ranked[rank_idx]
        pair_probability_original_order[original_feature_idx] = pair_probability_ranked[rank_idx]
        pair_log_posterior_original_order[original_feature_idx] = pair_log_posterior_ranked[rank_idx]
    }

    list(
        assignment = best_assignment_original_order,
        probability = pair_probability_original_order,
        log_posterior = pair_log_posterior_original_order
    )
}
