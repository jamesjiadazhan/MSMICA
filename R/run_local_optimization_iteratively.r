#' Iteratively run the feature-wise local optimization until convergence
#'
#' Repeats the Bayesian local optimization loop
#' (\code{local_optimization_per_mass()} → precursor-product filtering
#' → concentration-priority pruning) until no further assignments can be
#' made or \code{max_iter} is reached. At each iteration, newly assigned
#' (feature, metabolite) pairs are removed from the candidate pool before
#' the next round, so the pool strictly shrinks. The function stops early
#' if: the input pool is empty; no assignment survives the scoring and
#' filtering steps; the pool size did not decrease; or the iteration
#' count exceeds \code{max_iter}.
#'
#' This iterative design lets confident, unambiguous assignments be made
#' first and then used as anchors for subsequent rounds, gradually
#' resolving more ambiguous cases as the candidate space narrows.
#'
#' @param input_data Data frame of candidate (feature, metabolite) rows
#'   in the MSMICA scoring schema, split by \code{Mono_mass}.
#' @param rt_sigma Numeric. SD of the RT Gaussian passed to
#'   \code{local_optimization_per_mass()}.
#' @param pp_mu,pp_sigma Mean and SD of the Fisher-z correlation prior
#'   passed to \code{local_optimization_per_mass()}.
#' @param detail Logical. If \code{TRUE} and \code{output_folder} is set,
#'   intermediate input and result tables are written as CSV files.
#' @param output_folder Character path. Directory for intermediate CSV
#'   output when \code{detail = TRUE}.
#' @param max_iter Integer. Maximum number of optimization rounds before
#'   a hard stop. Default is 20.
#' @return A named list with \code{result} (bound data frame of all
#'   accepted assignments across all iterations, sorted by \code{mz}
#'   and \code{time}) and \code{remaining_data} (candidate rows that
#'   were never assigned).
#' @keywords internal
#' @noRd
run_local_optimization_iteratively = function(input_data, rt_sigma, pp_mu, pp_sigma,
                                              rxn_connection, cor_input, adduct_corr_time_thresh,
                                              detail = FALSE, output_folder = NULL, max_iter = 20,
                                              concentration_prior_weight = 1,
                                              feature_rank_prior_strength = 1,
                                              single_feature_prior_multiplier = 10,
                                              cross_mass_concentration_tiebreak = TRUE,
                                              use_concentration_tag = TRUE,
                                              local_optimization_assignment = "metabolite_greedy",
                                              local_rt_weight = 1,
                                              local_corr_weight = 1,
                                              metabolite_rt_weight = 1,
                                              metabolite_corr_weight = 1,
                                              metabolite_intensity_weight = 1,
                                              calibration_method = "current",
                                              evidence_calibration = NULL) {
    current_data = input_data
    all_results = list()
    iter = 1

    repeat {
        message("*******************************************************************************")
        message(paste0(
            "Performing local optimization per feature, iteration ", iter, ":"
        ))

        if (nrow(current_data) == 0) {
            message("No remaining input data. Stop.")
            break
        }

        # optional save of input
        if (detail == TRUE && !is.null(output_folder)) {
            write_csv(
                current_data,
                paste0(
                    output_folder, "/",
                    "MSMICA_local_optimization_per_mass_main_adduct_input_data_",
                    iter, ".csv"
                )
            )
        }

        # split by Mono_mass
        current_split = split(current_data, current_data$Mono_mass)

        # progress bar
        total_iterations = length(current_split)
        pb_local_optimization_per_mass = txtProgressBar(
            min = 0,
            max = total_iterations,
            style = 3
        )

        i = 0
        last_progress_updete = 1

        current_result_split = lapply(seq_along(current_split), function(j) {
            res = local_optimization_per_mass(
                mass_group_data = current_split[[j]],
                rt_sigma = rt_sigma,
                corr_mu = pp_mu,
                corr_sigma = pp_sigma,
                w_rt = local_rt_weight,
                w_corr = local_corr_weight,
                w_prior = concentration_prior_weight,
                single_feature_prior_multiplier = single_feature_prior_multiplier,
                feature_rank_prior_strength = feature_rank_prior_strength,
                assignment_strategy = local_optimization_assignment,
                calibration_method = calibration_method,
                evidence_calibration = evidence_calibration
            )
            setTxtProgressBar(pb_local_optimization_per_mass, j)
            res
        })

        close(pb_local_optimization_per_mass)

        # combine
        current_result = do.call(rbind, current_result_split)

        if (is.null(current_result) || nrow(current_result) == 0) {
            message("No new local optimization result found. Stop.")
            break
        }

        current_result = tibble(current_result)

        # join back
        current_joined = current_data %>%
            inner_join(
                current_result,
                by = c("mz_time", "InChIKey"),
                relationship = "many-to-many"
            )

        # process identification_method
        current_processed = process_identification_method(
            current_joined,
            use_concentration_tag = use_concentration_tag
        )

        message("Here is the current local optimization result:")
        print(current_processed)

        # group by mz_time and select the metabolites with the highest Concentration_average: this is because the previous steps use monoisotopic mass, but some metabolites have very close but not the same monoisotopic mass (like Butyrobetaine and Acetylcholine, but Butyrobetaine is much more abundant should be prioritized). 
        if (cross_mass_concentration_tiebreak) {
            current_filtered = current_processed %>%
                group_by(mz_time) %>%
                filter(Concentration_average == max(Concentration_average) | is.na(Concentration_average)) %>%
                ungroup()
        }
    
        # standardize output
        current_filtered = current_filtered %>%
            dplyr::select(-log_posterior) %>%
            rename(Probability = probability) %>%
            mutate(Probability = Probability * 100)

        # stop if nothing survives this iteration
        if (nrow(current_filtered) == 0) {
            message("No new filtered result found. Stop.")
            break
        }

        # save result
        if (detail == TRUE && !is.null(output_folder)) {
            write_csv(
                current_filtered,
                paste0(
                    output_folder, "/",
                    "MSMICA_local_optimization_per_mass_main_adduct_result_",
                    iter, ".csv"
                )
            )
        }

        # store result
        all_results[[iter]] = current_filtered

        # keep your exact logic here
        next_data = current_data %>%
            filter(
                !(mz_time %in% current_filtered$mz_time) &
                !(InChIKey %in% current_filtered$InChIKey)
            )

        # stop if no shrinkage
        if (nrow(next_data) == nrow(current_data)) {
            message("No further reduction in input data. Stop to avoid infinite loop.")
            break
        }

        current_data = next_data
        iter = iter + 1

        if (iter > max_iter) {
            message("Reached max_iter. Stop.")
            break
        }
    }

    combined_result = dplyr::bind_rows(all_results) %>%
        arrange(mz, time)

    return(list(
        result = combined_result,
        remaining_data = current_data
    ))
}
