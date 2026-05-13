#' Plackett-Luce rank log-likelihood
#'
#' Computes the log-likelihood of an observed ranking under the
#' Plackett-Luce model given log-scale preference scores (here, log
#' concentration priors). The input is assumed to be supplied in rank
#' order: the most-preferred item first. The log-likelihood is the sum
#' over positions of the log-softmax of the remaining items. Used by the
#' Bayesian local optimization to score how consistent a candidate
#' feature-to-metabolite assignment is with the relative abundance of the
#' candidate metabolites.
#'
#' @param log_concentration_values Numeric vector of log-scale scores
#'   (for MSMICA, log concentration priors) arranged in rank order.
#' @return A list with two components: \code{total}, the summed
#'   log-likelihood, and \code{per_rank}, the vector of per-position
#'   log-likelihoods. For vectors of length 0 or 1, both are returned as
#'   zero.
#' @keywords internal
#' @noRd
plackett_luce_rank_loglik = function(log_concentration_values) {

    scaled_log_concentration = log_concentration_values
    n_ranked = length(scaled_log_concentration)

    # Degenerate case: ranking of 0 or 1 items has zero log-likelihood
    if (n_ranked <= 1) {
        return(list(total = 0, per_rank = rep(0, n_ranked)))
    }

    per_rank_loglik = numeric(n_ranked)

    # At each rank, the PL probability of picking item rank_idx is
    # exp(score_rank_idx) / sum(exp(scores_of_remaining_items)).
    # In log space this reduces to score - logsumexp(remaining scores).
    for (rank_idx in seq_len(n_ranked)) {
        per_rank_loglik[rank_idx] = scaled_log_concentration[rank_idx] - logsumexp(scaled_log_concentration[rank_idx:n_ranked])
    }

    list(
        total = sum(per_rank_loglik),
        per_rank = per_rank_loglik
    )
}
