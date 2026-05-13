#' Backpropagation correlation analysis for MSMICA
#'
#' Merges the expanded precursor-product candidate pairs with their
#' computed correlations and summarizes the biological-network evidence
#' supporting each candidate product. Optionally applies
#' Benjamini-Hochberg FDR correction within each product InChIKey,
#' collapses bidirectional pair records by median significant
#' correlation, and counts total vs. significant correlations.
#' \code{duplicate_removal} controls whether precursors with a mix of
#' unique and redundant product links are pruned to the unique ones
#' (while keeping entirely-redundant precursors intact).
#'
#' @param MSMICA_col_names_connection_12 Data frame of expanded
#'   precursor-product candidate pairs produced by
#'   \code{precursor_product_correlation_processing()}.
#' @param cor_MSMICA Tibble of per-pair correlations to join in.
#' @param backpropagation_correlation_direction Character, one of
#'   \code{"positive"} (keep only positive correlations) or
#'   \code{"both"} (use absolute value).
#' @param FDR_correction Logical. Apply BH FDR correction of p-values
#'   within each product InChIKey if \code{TRUE}.
#' @param duplicate_removal Logical. If \code{TRUE}, prune a precursor's
#'   redundant product links while keeping its unique ones.
#' @return A named list with \code{MSMICA_col_names_connection_identified_final_4}
#'   (the summarized per-product evidence table) and
#'   \code{Backpropagation_input_clean_detailed} (the row-level joined
#'   table, useful for diagnostics).
#' @keywords internal
#' @noRd
backpropagation_correlation_analysis = function(MSMICA_col_names_connection_12, cor_MSMICA, backpropagation_correlation_direction, FDR_correction, duplicate_removal) {
    # Join per-pair correlations onto the candidate precursor-product table
    MSMICA_col_names_connection_12_2 = inner_join(MSMICA_col_names_connection_12, cor_MSMICA, by = c("mz_time_1", "mz_time_2"), relationship = "many-to-many")

    # make correlation values absolute if backpropagation_correlation_direction is "both"
    if (backpropagation_correlation_direction == "both") {
        MSMICA_col_names_connection_12_2$correlation = abs(MSMICA_col_names_connection_12_2$correlation)
        print("backpropagation correlation direction is both: positive and negative correlations are considered")
    } else if (backpropagation_correlation_direction == "positive") {
        MSMICA_col_names_connection_12_2 = MSMICA_col_names_connection_12_2 %>%
            filter(correlation > 0)
        print("backpropagation correlation direction is positive: only positive correlations are considered")
    }

    # perform the FDR correction (Benjamini-Hochberg correction)
    if (FDR_correction == TRUE) {
        MSMICA_col_names_connection_12_2 = MSMICA_col_names_connection_12_2 %>%
            group_by(connection_2_InChIKey) %>%
            mutate(p_value = p.adjust(p_value, method = "BH")) %>%
            ungroup()
    }

    # arrange the correlation from high to low by each mz_1 and connection_2_InChIKey
    MSMICA_col_names_connection_12_2 = MSMICA_col_names_connection_12_2 %>%
        arrange(mz_1, connection_2_InChIKey, desc(correlation))

    # filter out those with NA in correlation
    MSMICA_col_names_connection_12_2 = MSMICA_col_names_connection_12_2 %>%
        filter(!is.na(correlation))

    # arrange by identified_Name_1 and mz_1
    MSMICA_col_names_connection_12_2 = MSMICA_col_names_connection_12_2 %>%
        arrange(identified_Name_1) %>%
        arrange(mz_1)

    # note, here the MSMICA_col_names_connection_12_2 is bidirectional: connection_1 and connection_2 have precursor-product/transporter relationship. Also, connection_1's records are repeated in connection_2's records, so there is no need to do additional things.

    # group by mz and time and remove all rows with duplicated mz_time_identified_final values 
    # mz_time_1 = precursor feature
    # mz_time_2 = product feature
    # duplicate_count == 1 means that precursor-product pair is unique
    # duplicate_count > 1 means that pair appears multiple times
    Backpropagation_input_clean_detailed = MSMICA_col_names_connection_12_2 %>%
        mutate(mz_time_1 = paste0(mz_1, "_", time_1)) %>%
        # relocate mz_time right after time
        relocate(mz_time_1, .after = time_1) %>%
        group_by(mz_time_1, mz_time_2) %>%
        mutate(duplicate_count = n()) %>%
        ungroup()

    # duplicate_removal controls whether the function keeps all precursor-product matches, or prunes repeated product assignments for each precursor before the downstream correlation summary
    if (duplicate_removal == FALSE) {
        message("Duplicated features used as different precursors or products for each feature's backpropagation are kept")
    } else if (duplicate_removal == TRUE) {
        message("Duplicated features used as different precursors or products for each feature's backpropagation are removed")

        # For a given precursor mz_time_1, if it has any unique product links, keep only those unique links (duplicate_count == 1).
        # If that precursor has no unique links at all, keep all of its duplicated links.
        # In other words, its function is to prefer cleaner, non-redundant precursor-product relationships when available, and avoid deleting an entire precursor’s backpropagation set if everything in that set is duplicated.
        Backpropagation_input_clean_detailed = Backpropagation_input_clean_detailed %>%
            group_by(mz_time_1, mz_time_2) %>%
            mutate(any_unique = any(duplicate_count == 1), all_duplicate = all(duplicate_count > 1)) %>%
            ungroup()

        # group by mz_time_1, and if there are any rows with any_unique == TRUE, then select the rows with any_unique == TRUE and duplicate_count == 1
        Backpropagation_input_clean_detailed = Backpropagation_input_clean_detailed %>%
            # 1. Group ONLY by the connection_1 (mz_time_1) to assess the whole set for that ion
            group_by(mz_time_1) %>%
            mutate(
                # Check if this entire precursor group contains ANY row with duplicate_count == 1
                group_has_unique = any(duplicate_count == 1)
            ) %>%
            # 2. Apply the selection logic
            filter(
                # Case A: Keep the row if it is a unique connection
                duplicate_count == 1 | 
                # Case B: OR keep the row if the group has NO unique connections at all
                !group_has_unique
            ) %>%
            # Clean up helper columns and grouping
            ungroup() %>%
            select(-group_has_unique)
    }

    # Summarize the correlation to median correlation and then count the total number of correlations and significant correlations
    Backpropagation_input_clean_detailed_2 = Backpropagation_input_clean_detailed %>%
        # select connection_2's records
        dplyr::select(KEGG_ID, HMDB_ID, connection_2_InChIKey, identified_Name_2, identification_type_1, ion_mode, Mono_mass, Adduct_2, mz_2, time_2, mz_time_2, identification_method_2, time_predicted, time_difference, Concentration_average, correlation, p_value) %>%
        # group by Mono_mass, Confirmed_Name, Adduct_annotated, and then create match_category = "single" if there is only 1 row for each Mono_mass, Confirmed_Name, Adduct_annotated
        group_by(Mono_mass, identified_Name_2, Adduct_2) %>%
        mutate(match_category = if_else(n() == 1, "single", "multiple")) %>%
        ungroup() %>%
        rename(InChIKey = connection_2_InChIKey, identified_Name = identified_Name_2, identification_type = identification_type_1, Adduct = Adduct_2, mz = mz_2, time = time_2, mz_time = mz_time_2, identification_method = identification_method_2) %>%
        group_by(KEGG_ID, HMDB_ID, InChIKey, identified_Name, Mono_mass, identification_type, ion_mode, Adduct, mz, time, mz_time, identification_method, match_category, time_predicted, time_difference, Concentration_average) %>%
        summarize(
            ## calculate the median of the significant correlations (p_value < 0.05)
            correlation = median(correlation[p_value < 0.05], na.rm = TRUE), 
            ## count the total number of correlations
            correlation_count = n(),
            ## count the total number of significant correlations
            significant_correlation_count = sum(p_value < 0.05, na.rm = TRUE),
            .groups = "keep") %>%
        ungroup()

    # select specific columns from final_results
    Backpropagation_input_clean_detailed_3 = Backpropagation_input_clean_detailed_2 %>%
        dplyr::select(InChIKey, KEGG_ID, HMDB_ID, ion_mode, identification_type, identified_Name, Adduct, mz, time, correlation, correlation_count, significant_correlation_count, identification_method, match_category, time_predicted, time_difference, Concentration_average)

    # arrange by mz and time
    Backpropagation_input_clean_detailed_3 = Backpropagation_input_clean_detailed_3 %>%
        arrange(mz, time)

    # create a new column, mz_time, to combine mz and time (thus, each mz_time is a unique feature)
    Backpropagation_input_clean_detailed_3$mz_time = paste0(Backpropagation_input_clean_detailed_3$mz, "_", Backpropagation_input_clean_detailed_3$time)

    # move mz_time right after time
    Backpropagation_input_clean_detailed_3 = Backpropagation_input_clean_detailed_3 %>%
        relocate(mz_time, .after = time)

    return(list(MSMICA_col_names_connection_identified_final_4 = Backpropagation_input_clean_detailed_3, Backpropagation_input_clean_detailed = Backpropagation_input_clean_detailed))
}