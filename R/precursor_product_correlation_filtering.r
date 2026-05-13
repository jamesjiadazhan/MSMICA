#' Filter duplicated candidate identifications using precursor-product evidence
#'
#' Takes the output of \code{local_optimization_per_feature()} and
#' resolves cases where a single metabolite-adduct combination is still
#' associated with multiple candidate features. For each duplicated
#' group the function re-scores candidates by combining:
#' the number of significant precursor-product correlations, the spread
#' of retention times relative to the adduct-correlation time window,
#' and a per-metabolite Bayesian re-scoring via
#' \code{local_optimization_per_metabolite()}. A second pass handles
#' residual duplicates introduced when the filtered set is merged back
#' with the original identification table.
#'
#' @param MSMICA_local_optimization_per_feature_result Data frame of
#'   candidate identifications from the feature-wise local optimization
#'   step.
#' @return A filtered version of the input with at most one candidate
#'   per metabolite-adduct combination, ordered by \code{mz} and
#'   \code{time}.
#' @keywords internal
#' @noRd
precursor_product_correlation_filtering = function(MSMICA_local_optimization_per_feature_result,
                                                   rxn_connection,
                                                   cor_input,
                                                   adduct_corr_time_thresh,
                                                   rt_sigma,
                                                   pp_mu,
                                                   pp_sigma) {
    # Compute precursor-product correlations for every still-ambiguous pair
    MSMICA_local_optimization_per_feature_result_pp = precursor_product_correlation_processing(
        MSMICA_local_optimization_per_feature_result,
        rxn_connection = rxn_connection,
        cor_input = cor_input
    )

    # add correlation values
    ## extract the $MSMICA_col_names_connection_12 from the MSMICA_local_optimization_per_feature_result_pp and remove the correlation column
    MSMICA_local_optimization_per_feature_result_pp_12 = MSMICA_local_optimization_per_feature_result_pp$MSMICA_col_names_connection_12
    # remove the correlation column from MSMICA_local_optimization_per_feature_result_pp_12
    MSMICA_local_optimization_per_feature_result_pp_12 = MSMICA_local_optimization_per_feature_result_pp_12 %>%
        dplyr::select(-correlation)

    # add the correlation values to the MSMICA_local_optimization_per_feature_result_pp_12
    MSMICA_local_optimization_per_feature_result_pp_backpropagation = inner_join(MSMICA_local_optimization_per_feature_result_pp_12, MSMICA_local_optimization_per_feature_result_pp$cor_MSMICA, by = c("mz_time_1", "mz_time_2"), relationship = "many-to-many")

    # Summarize the correlation to median correlation and then count the total number of correlations and significant correlations
    MSMICA_local_optimization_per_feature_result_pp_backpropagation_2 = MSMICA_local_optimization_per_feature_result_pp_backpropagation %>%
        # select connection_2's records
        dplyr::select(KEGG_ID, HMDB_ID, connection_2_InChIKey, identified_Name_2, identification_type_1, ion_mode, Mono_mass, Adduct_2, mz_2, time_2, mz_time_2, identification_method_2, time_predicted, time_difference, Concentration_average, correlation, p_value, mean_intensity, probability) %>%
        # group by Mono_mass, Confirmed_Name, Adduct_annotated, and then create match_category = "single" if there is only 1 row for each Mono_mass, Confirmed_Name, Adduct_annotated
        group_by(Mono_mass, identified_Name_2, Adduct_2) %>%
        mutate(match_category = if_else(n() == 1, "single", "multiple")) %>%
        ungroup() %>%
        rename(InChIKey = connection_2_InChIKey, identified_Name = identified_Name_2, identification_type = identification_type_1, Adduct = Adduct_2, mz = mz_2, time = time_2, mz_time = mz_time_2, identification_method = identification_method_2) %>%
        group_by(KEGG_ID, HMDB_ID, InChIKey, identified_Name, Mono_mass, identification_type, ion_mode, Adduct, mz, time, mz_time, identification_method, match_category, time_predicted, time_difference, Concentration_average, mean_intensity, probability) %>%
        summarize(
            ## calculate the median of the significant correlations (p_value < 0.05)
            correlation = median(correlation[p_value < 0.05], na.rm = TRUE), 
            ## count the total number of correlations
            correlation_count = n(),
            ## count the total number of significant correlations
            significant_correlation_count = sum(p_value < 0.05, na.rm = TRUE),
            .groups = "keep") %>%
        ungroup()
    
    # group by identified_Name and keep the rows with the highest significant_correlation_count
    MSMICA_local_optimization_per_feature_result_pp_backpropagation_3 = MSMICA_local_optimization_per_feature_result_pp_backpropagation_2 %>%
        group_by(identified_Name) %>%
        # calculate the max retention time difference for the same identified_Name
        mutate(
            max_time_difference = if (all(is.na(time_difference))) NA_real_ else max(time_difference, na.rm = TRUE)
        ) %>%
        filter(
            if (all(is.na(significant_correlation_count))) {
                TRUE
            } else if (!is.na(first(max_time_difference)) && first(max_time_difference) < adduct_corr_time_thresh && any(grepl("clustering of adducts and isotopes", identification_method))) {
                TRUE
            } else {
                significant_correlation_count == max(significant_correlation_count, na.rm = TRUE)
            }
        ) %>%
        ungroup()

    # group by Mono_mass and filter out the rows with time_difference > 200 seconds
    MSMICA_local_optimization_per_feature_result_pp_backpropagation_3 = MSMICA_local_optimization_per_feature_result_pp_backpropagation_3 %>%
        group_by(Mono_mass) %>%
        mutate(all_time_difference_200 = all(time_difference >= 200)) %>%
        ungroup() %>%
        filter(all_time_difference_200 | time_difference < 200) %>%
        dplyr::select(-all_time_difference_200)
    
    # make the correlation value absolute
    MSMICA_local_optimization_per_feature_result_pp_backpropagation_3$correlation = abs(MSMICA_local_optimization_per_feature_result_pp_backpropagation_3$correlation)

    # find where identified_Name is duplicated (multiple adducts)
    identified_Name_duplicated = MSMICA_local_optimization_per_feature_result_pp_backpropagation_3 %>%
        group_by(identified_Name, Adduct) %>%
        filter(n() > 1) %>%
        ungroup()
    
    # if identified_Name_duplicated is not empty, then select the best feature for each metabolite based on time_difference, precursor-product/transporter correlation, mean_intensity
    if (nrow(identified_Name_duplicated) > 0) {
        # select the best feature for each metabolite based on time_difference, precursor-product/transporter correlation, mean_intensity
        best_feature_results = local_optimization_per_metabolite(
            identified_Name_duplicated,
            rt_sigma = rt_sigma,
            corr_mu = pp_mu,
            corr_sigma = pp_sigma,
            w_rt = 1,
            w_corr = 1,
            w_intensity = 1
        )

        # within best_feature_results, select where best_feature == TRUE
        best_feature_results_TRUE = best_feature_results %>%
            filter(best_feature == TRUE)

        # within MSMICA_local_optimization_per_feature_result_pp_backpropagation_3, exclude InChIKey and mz_time that are in best_feature_results_TRUE
        MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered = MSMICA_local_optimization_per_feature_result_pp_backpropagation_3 %>%
            filter(!(InChIKey %in% best_feature_results_TRUE$InChIKey | mz_time %in% best_feature_results_TRUE$mz_time))

        # remove log_rt, log_corr, log_intensity_prior columns, best_feature, and log_posterior columns in the best_feature_results_TRUE
        best_feature_results_TRUE = best_feature_results_TRUE %>%
            dplyr::select(-log_rt, -log_corr, -log_intensity_prior, -best_feature, -log_posterior)

        # combine best_feature_results_TRUE and MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered
        MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered = rbind(best_feature_results_TRUE, MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered)
    } else {
        # if identified_Name_duplicated is empty, then MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered is the same as MSMICA_local_optimization_per_feature_result_pp_backpropagation_3
        MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered = MSMICA_local_optimization_per_feature_result_pp_backpropagation_3
    }

    # remove duplicates
    MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered = MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered %>%
        distinct() %>%
        arrange(mz, time)

    # within MSMICA_local_optimization_per_feature_result, keep only the rows with mz_time and InChIKey that are in MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered
    MSMICA_local_optimization_per_feature_result_pp_filtered = MSMICA_local_optimization_per_feature_result %>%
        filter(mz_time %in% MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered$mz_time & InChIKey %in% MSMICA_local_optimization_per_feature_result_pp_backpropagation_3_filtered$InChIKey)

    # within MSMICA_local_optimization_per_feature_result, remove the rows with mz_time and InChIKey that are in MSMICA_local_optimization_per_feature_result_pp_filtered
    MSMICA_local_optimization_per_feature_result_other_filtered = MSMICA_local_optimization_per_feature_result %>%
        filter(!(mz_time %in% MSMICA_local_optimization_per_feature_result_pp_filtered$mz_time | InChIKey %in% MSMICA_local_optimization_per_feature_result_pp_filtered$InChIKey))

    # combine MSMICA_local_optimization_per_feature_result_pp_filtered and MSMICA_local_optimization_per_feature_result_other_filtered
    MSMICA_local_optimization_per_feature_result_filtered = rbind(MSMICA_local_optimization_per_feature_result_pp_filtered, MSMICA_local_optimization_per_feature_result_other_filtered) %>% 
        arrange(mz, time)

    # find where identified_Name is duplicated (multiple adducts)
    identified_Name_duplicated_2 = MSMICA_local_optimization_per_feature_result_filtered %>%
        group_by(identified_Name, Adduct) %>%
        filter(n() > 1) %>%
        ungroup()

    # if identified_Name_duplicated_2 is not empty, then select the best feature for each metabolite based on time_difference, precursor-product/transporter correlation, mean_intensity
    if (nrow(identified_Name_duplicated_2) > 0) {

        # select the best feature for each metabolite based on time_difference, precursor-product/transporter correlation, mean_intensity
        best_feature_results_2 = local_optimization_per_metabolite(
            identified_Name_duplicated_2,
            rt_sigma = rt_sigma,
            corr_mu = pp_mu,
            corr_sigma = pp_sigma,
            w_rt = 1,
            w_corr = 1,
            w_intensity = 1
        )

        # within best_feature_results_2, select where best_feature == TRUE
        best_feature_results_2_TRUE = best_feature_results_2 %>%
            filter(best_feature == TRUE)

        # within MSMICA_local_optimization_per_feature_result_filtered, exclude InChIKey and mz_time that are in best_feature_results_2_TRUE
        MSMICA_local_optimization_per_feature_result_filtered_2 = MSMICA_local_optimization_per_feature_result_filtered %>%
            filter(!(InChIKey %in% best_feature_results_2_TRUE$InChIKey | mz_time %in% best_feature_results_2_TRUE$mz_time))

        # remove log_rt, log_corr, log_intensity_prior columns, best_feature, and log_posterior columns in the best_feature_results_2_TRUE
        best_feature_results_2_TRUE = best_feature_results_2_TRUE %>%
            dplyr::select(-log_rt, -log_corr, -log_intensity_prior, -best_feature)
        
        # combine best_feature_results_2_TRUE and MSMICA_local_optimization_per_feature_result_filtered_2
        MSMICA_local_optimization_per_feature_result_filtered_2 = rbind(best_feature_results_2_TRUE, MSMICA_local_optimization_per_feature_result_filtered_2)

        # remove duplicates
        MSMICA_local_optimization_per_feature_result_filtered_2 = MSMICA_local_optimization_per_feature_result_filtered_2 %>%
            distinct() %>%
            arrange(mz, time)
    } else {
        # if identified_Name_duplicated_2 is empty, then MSMICA_local_optimization_per_feature_result_filtered_2 is the same as MSMICA_local_optimization_per_feature_result_filtered
        MSMICA_local_optimization_per_feature_result_filtered_2 = MSMICA_local_optimization_per_feature_result_filtered
    }

    return(MSMICA_local_optimization_per_feature_result_filtered_2)
}