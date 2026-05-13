#' Select the best feature for each metabolite across a mass group
#'
#' When the same metabolite-adduct combination has been matched to more
#' than one observed feature (e.g. two near-co-eluting peaks of the same
#' m/z), this function picks the single best feature for each metabolite
#' using a log-space combination of three independent evidence terms:
#' retention-time agreement (Gaussian on \code{time_difference}),
#' precursor-product/transporter correlation (Gaussian on the
#' Fisher-z-transformed Spearman r), and a softmax prior over observed
#' intensities. It calls
#' \code{exact_bayesian_feature_selection_for_metabolite()} to convert
#' the combined log-score vector into a posterior probability and the
#' MAP feature choice per metabolite.
#'
#' @param mass_group_data Data frame of candidate (feature, metabolite)
#'   rows for one or more metabolites, typically the subset where the
#'   same metabolite appears on multiple features. Required columns:
#'   \code{InChIKey}, \code{mz_time}, \code{time_difference},
#'   \code{correlation}, \code{mean_intensity} (or
#'   \code{log_mean_intensity}), \code{mz}, \code{time}.
#' @param rt_sigma Standard deviation of the RT-difference Gaussian
#'   likelihood; defaults to the global \code{rt_mapping_sigma_predret}.
#' @param corr_mu,corr_sigma Mean and SD of the Fisher-z correlation
#'   prior; default to globals \code{pp_mu} and \code{pp_sigma}.
#' @param w_rt,w_corr,w_intensity Weights on the three evidence terms.
#'   All default to 1 (equal weighting).
#' @param intensity_temperature Softmax temperature for the
#'   intensity-based within-metabolite prior. Values > 1 flatten
#'   intensity differences; values < 1 sharpen them. Default is 1.
#' @return A data frame in the same format as \code{mass_group_data}
#'   with additional columns \code{log_rt}, \code{log_corr},
#'   \code{log_intensity_prior}, \code{log_posterior},
#'   \code{probability}, and \code{best_feature}. Returns \code{NULL}
#'   if no viable candidates remain.
#' @keywords internal
#' @noRd
local_optimization_per_metabolite = function(mass_group_data,
                                            rt_sigma = rt_mapping_sigma_predret,
                                            corr_mu = pp_mu,
                                            corr_sigma = pp_sigma,
                                            w_rt = 1,
                                            w_corr = 1,
                                            w_intensity = 1,
                                            intensity_temperature = 1) {

    # remove exact duplicated rows caused by rounded mz/time collisions
    mass_group_data = mass_group_data %>%
        mutate(InChIKey_mz_time = paste0(InChIKey, "_", mz, "_", time)) %>%
        group_by(InChIKey_mz_time) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup() %>%
        dplyr::select(-InChIKey_mz_time)

    metabolites = unique(mass_group_data$InChIKey)

    if (length(metabolites) == 0) {
        return(NULL)
    }

    result_list = vector("list", length(metabolites))

    for (j in seq_along(metabolites)) {
        met = metabolites[j]

        metabolite_data = mass_group_data %>%
            dplyr::filter(InChIKey == met) %>%
            dplyr::distinct(mz_time, .keep_all = TRUE)

        features = metabolite_data$mz_time
        n_features = length(features)

        if (n_features == 0) {
            result_list[[j]] = NULL
            next
        }

        # initialize
        log_rt_vec = rep(NA_real_, n_features)
        log_corr_vec = rep(NA_real_, n_features)
        log_intensity_prior_vec = rep(NA_real_, n_features)

        # ----------------------------
        # 1. RT and correlation evidence
        # ----------------------------
        for (i in seq_len(n_features)) {
            pair_data = metabolite_data[i, ]

            td = pair_data$time_difference[1]

            # A. RT likelihood
            p_rt = dnorm(td, mean = 0, sd = rt_sigma)
            p_rt = max(p_rt, 1e-300)
            log_rt_vec[i] = log(p_rt)

            # B. Correlation likelihood
            if (is.na(pair_data$correlation[1])) {
                log_corr_vec[i] = 0
            } else {
                corr_value = min(max(pair_data$correlation[1], -0.999999), 0.999999)
                z_obs = atanh(corr_value)
                p_corr = dnorm(z_obs, mean = corr_mu, sd = corr_sigma)
                p_corr = max(p_corr, 1e-300)
                log_corr_vec[i] = log(p_corr)
            }
        }

        # ----------------------------
        # 2. Intensity-based within-metabolite prior
        # ----------------------------
        if ("log_mean_intensity" %in% colnames(metabolite_data)) {
            intensity_score = metabolite_data$log_mean_intensity
        } else {
            intensity_score = log(pmax(metabolite_data$mean_intensity, 1))
        }

        # soften or sharpen the intensity effect
        intensity_score = intensity_score / intensity_temperature

        # convert to within-metabolite prior by softmax
        log_intensity_prior_vec = intensity_score - logsumexp(intensity_score)

        # ----------------------------
        # 3. Combine all evidence
        # ----------------------------
        feature_log_score =
            w_rt * log_rt_vec +
            w_corr * log_corr_vec +
            w_intensity * log_intensity_prior_vec

        selection_result = exact_bayesian_feature_selection_for_metabolite(feature_log_score)

        if (is.null(selection_result)) {
            result_list[[j]] = NULL
            next
        }

        metabolite_result = metabolite_data %>%
            mutate(
                log_rt = log_rt_vec,
                log_corr = log_corr_vec,
                log_intensity_prior = log_intensity_prior_vec,
                log_posterior = feature_log_score,
                probability = selection_result$all_probabilities,
                best_feature = mz_time == features[selection_result$best_feature_index]
            )

        result_list[[j]] = metabolite_result
    }

    results = dplyr::bind_rows(result_list)

    if (nrow(results) == 0) {
        return(NULL)
    }

    results
}