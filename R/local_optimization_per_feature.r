#' Feature-wise Bayesian local optimization for a single monoisotopic mass
#'
#' Scores every (feature, metabolite) pair in a single-monoisotopic-mass
#' group using three sources of evidence -- retention-time agreement
#' (Gaussian on \code{time_difference}), biological-network correlation
#' (Gaussian on the Fisher z-transformed correlation), and the
#' biospecimen-specific metabolite concentration prior. The
#' log-likelihood matrix is then handed to
#' \code{exact_bayesian_assignment()} to pick one metabolite per feature.
#' When only one candidate feature has matched m/z, the concentration
#' prior is boosted and RT / correlation weights are halved to reflect
#' the reduced discriminatory power.
#'
#' @param mass_group_data Data frame with candidate rows for one
#'   monoisotopic mass group, including \code{mz_time}, \code{InChIKey},
#'   \code{time_difference}, \code{correlation},
#'   \code{log_concentration_prior}, \code{mean_intensity} and
#'   \code{log_mean_intensity}.
#' @param rt_sigma Standard deviation of the RT Gaussian (defaults to
#'   the global \code{rt_mapping_sigma_predret}).
#' @param corr_mu,corr_sigma Mean and SD of the Fisher-z correlation
#'   prior (defaults to globals \code{pp_mu}, \code{pp_sigma}).
#' @param w_rt,w_corr,w_prior Weights on the three evidence terms.
#' @param single_feature_prior_multiplier Boost factor applied to
#'   \code{w_prior} when only one candidate feature is present.
#' @param feature_rank_prior_strength Strength of the feature-rank
#'   prior passed through to \code{exact_bayesian_assignment()}.
#' @return \code{NULL} if no candidate assignment is possible;
#'   otherwise a data frame with \code{mz_time}, assigned
#'   \code{InChIKey}, \code{log_posterior}, and posterior
#'   \code{probability}.
#' @keywords internal
#' @noRd
local_optimization_per_feature = function(mass_group_data,
                                            rt_sigma = rt_mapping_sigma_predret,
                                            corr_mu = pp_mu,
                                            corr_sigma = pp_sigma,
                                            w_rt = 1,
                                            w_corr = 1,
                                            w_prior = 1,
                                            single_feature_prior_multiplier = 10,
                                            feature_rank_prior_strength = 1
                                            ) {

    # All rows in mass_group_data share one monoisotopic mass
    current_monomass = unique(mass_group_data$Mono_mass)

    # group by InChIKey, mz, and time and then filter the row with the max mean_intensity: this is to avoid the bug caused by the features with exactly the mz and time after rounding
    mass_group_data = mass_group_data %>%
        mutate(
            InChIKey_mz_time = paste0(InChIKey, "_", mz, "_", time)
        )  %>%
        group_by(InChIKey_mz_time) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup() %>%
        ## remove the InChIKey_mz_time column
        dplyr::select(-InChIKey_mz_time)

    features = unique(mass_group_data$mz_time)
    metabolites = unique(mass_group_data$InChIKey)

    N = length(features)
    M = length(metabolites)

    if (N == 0 || M == 0) {
        return(NULL)
    }

    # boost metabolite concentration prior by multiplying it with the single_feature_prior_multiplier and reduce the weight of the retention time likelihood by 100% when there is exactly 1 observed feature with matched m/z
    if (N == 1) {
        w_prior_effective = w_prior * single_feature_prior_multiplier
        w_rt_effective = w_rt * 0.5
        w_corr_effective = w_corr * 0.5
    } else {
        w_prior_effective = w_prior
        w_rt_effective = w_rt
        w_corr_effective = w_corr
    }

    # ----------------------------
    # 1. Initialize score matrix
    # ----------------------------
    log_lik_matrix = matrix(NA_real_, nrow = N, ncol = M)
    rownames(log_lik_matrix) = features
    colnames(log_lik_matrix) = metabolites

    feature_metadata = mass_group_data %>%
        distinct(mz_time, mean_intensity, log_mean_intensity) %>%
        mutate(feature_order = match(mz_time, features)) %>%
        arrange(feature_order)

    metabolite_metadata = mass_group_data %>%
        distinct(InChIKey, log_concentration_prior) %>%
        mutate(metabolite_order = match(InChIKey, metabolites)) %>%
        arrange(metabolite_order)

    # ----------------------------
    # 2. Score each feature-metabolite pair
    # ----------------------------
    for (i in seq_along(features)) {
        feat = features[i]

        feature_data = mass_group_data %>%
            dplyr::filter(mz_time == feat)

        for (j in seq_along(metabolites)) {
            met = metabolites[j]

            pair_data = feature_data %>%
                dplyr::filter(InChIKey == met)

            if (nrow(pair_data) == 0) {
                next
            }

            td = pair_data$time_difference[1]

            # A. RT likelihood
            p_rt = dnorm(td, mean = 0, sd = rt_sigma)
            p_rt = max(p_rt, 1e-300)
            log_rt = log(p_rt)

            # B. Correlation likelihood
            if (is.na(pair_data$correlation[1])) {
                log_corr = 0
            } else {
                corr_value = min(max(pair_data$correlation[1], -0.999999), 0.999999)
                z_obs = atanh(corr_value)
                p_corr = dnorm(z_obs, mean = corr_mu, sd = corr_sigma)
                p_corr = max(p_corr, 1e-300)
                log_corr = log(p_corr)
            }

            # C. Metabolite concentration prior
            log_prior = pair_data$log_concentration_prior[1]

            log_lik_matrix[i, j] =
                w_rt_effective * log_rt +
                w_corr_effective * log_corr +
                w_prior_effective * log_prior
        }
    }

    # ----------------------------
    # 3. Feature-wise assignment
    # ----------------------------
    feature_rank_order = feature_metadata %>%
        arrange(desc(log_mean_intensity), mz_time) %>%
        pull(feature_order)

    assignment_result = exact_bayesian_assignment(
        local_log_score_matrix = log_lik_matrix,
        feature_rank_order = feature_rank_order,
        feature_log_intensity = feature_metadata$log_mean_intensity,
        metabolite_log_concentration = metabolite_metadata$log_concentration_prior,
        feature_rank_prior_strength = feature_rank_prior_strength
    )

    if (is.null(assignment_result)) {
        return(NULL)
    }

    results = data.frame(
        mz_time = features,
        InChIKey = metabolites[assignment_result$assignment],
        log_posterior = assignment_result$log_posterior,
        probability = assignment_result$probability,
        stringsAsFactors = FALSE
    )

    return(results)
}
