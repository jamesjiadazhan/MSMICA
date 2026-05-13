#' Compute precursor-product correlations for backpropagation
#'
#' Second stage of the precursor-product correlation pipeline. Given a
#' table of already-identified metabolites, this helper expands each one
#' into its candidate reactive partners using \code{reaction_connection},
#' builds the full set of unique \code{(mz_time_1, mz_time_2)} feature
#' pairs, and computes a Spearman correlation (with p-value) for each
#' pair against the intensity matrix \code{MSMICA_cor_input}. The
#' returned list feeds the backpropagation step, which uses these
#' correlations to adjust identification confidence for connected
#' features.
#'
#' Relies on \code{reaction_connection} and \code{MSMICA_cor_input}
#' being available in the calling scope.
#'
#' @param identified_metabolite_final Data frame of already-identified
#'   metabolites with the MSMICA output schema (including columns
#'   \code{InChIKey}, \code{mz}, \code{time}, \code{Adduct}, etc.).
#' @return A named list with \code{cor_MSMICA} (a tibble of
#'   \code{mz_time_1}, \code{mz_time_2}, \code{correlation},
#'   \code{p_value}) and \code{MSMICA_col_names_connection_12} (the
#'   expanded candidate-pair table used to produce it).
#' @keywords internal
#' @noRd
precursor_product_correlation_processing = function(identified_metabolite_final,
                                                    rxn_connection,
                                                    cor_input) {
    # Defensive numeric coercion: upstream sources can yield character columns
    identified_metabolite_final$mz = as.numeric(identified_metabolite_final$mz)
    identified_metabolite_final$time = as.numeric(identified_metabolite_final$time)

    # arrange the identified_metabolite_final_clustering by mz and time
    identified_metabolite_final = identified_metabolite_final %>%
        arrange(mz, time)

    # inner join the MSMICA_TCA_col_names_connection with the reaction_connection
    MSMICA_col_names_connection = identified_metabolite_final %>%
        inner_join(rxn_connection, by = c("InChIKey" = "connection_1_InChIKey"), relationship = "many-to-many")

    # Arrange the data by the mz column in ascending order
    MSMICA_col_names_connection_simplifed_2 = MSMICA_col_names_connection %>%
        arrange(mz)

    # select only the metabolites with the precursor-product relationships with the metabolites in the identified_metabolite_final
    # MSMICA_col_names_connection_simplifed_2 is the metabolites with relationships with the other metabolites that have precursor-product reactions
    MSMICA_col_names_connection_simplifed_2 = MSMICA_col_names_connection_simplifed_2[, c("InChIKey", "identification_type", "identified_Name", "Adduct", "mz", "time","identification_method", "connection_2_InChIKey", "react_id", "enzyme_transporter", "source")]

    # remove the duplicates based on "identification_type", "identified_Name", "Adduct", "mz", "time", "identification_method", "connection_2_InChIKey",
    MSMICA_col_names_connection_simplifed_2 = MSMICA_col_names_connection_simplifed_2 %>%
        distinct(identification_type, identified_Name, Adduct, mz, time, identification_method, connection_2_InChIKey, .keep_all = TRUE)

    # MSMICA_col_names_connection_simplifed_2
    MSMICA_col_names_connection_simplifed = MSMICA_col_names_connection_simplifed_2

    # round mz to 4 decimal places and time to integers
    MSMICA_col_names_connection_simplifed$mz = round(MSMICA_col_names_connection_simplifed$mz, 4)
    MSMICA_col_names_connection_simplifed$time = round(MSMICA_col_names_connection_simplifed$time, 0)
    # add mz_time to the MSMICA_col_names_connection_simplifed
    MSMICA_col_names_connection_simplifed$mz_time = paste0(MSMICA_col_names_connection_simplifed$mz, "_", MSMICA_col_names_connection_simplifed$time)

    # round mz to 4 decimal places and time to integers
    identified_metabolite_final$mz = round(identified_metabolite_final$mz, 4)
    identified_metabolite_final$time = round(identified_metabolite_final$time, 0)
    # add mz_time to the identified_metabolite_final
    identified_metabolite_final$mz_time = paste0(identified_metabolite_final$mz, "_", identified_metabolite_final$time)

    # inner join the MSMICA_col_names_connection_simplifed with the identified_metabolite_final
    MSMICA_col_names_connection_12 = inner_join(
        MSMICA_col_names_connection_simplifed, 
        identified_metabolite_final, 
        by = c("connection_2_InChIKey" = "InChIKey"),
        suffix = c("_1", "_2"),
        relationship = "many-to-many"
        )

    # remove those precursor-product/transporter correlation values that are self-correlation values
    ## remove where InChIKey == Connection_2_InChIKey (self-correlation values) and feature_1 == feature_1_final
    MSMICA_col_names_connection_12 = MSMICA_col_names_connection_12 %>%
        filter(InChIKey != connection_2_InChIKey) %>%
        mutate(feature_1 = paste0(mz_1, "_", time_1)) %>%
        mutate(feature_2 = paste0(mz_2, "_", time_2)) %>%
        filter(feature_1 != feature_2) %>%
        dplyr::select(-feature_1, -feature_2)

    # Filter out rows with NA in mz_time_1 or mz_time_2 early
    valid_rows = MSMICA_col_names_connection_12 %>%
        filter(!is.na(mz_time_1) & !is.na(mz_time_2))
    
    # remove duplicated rows with the same mz_time_1 and mz_time_2
    valid_rows = valid_rows %>%
        distinct(mz_time_1, mz_time_2, .keep_all = TRUE)

    # Extract the columns of interest from valid_rows
    identified_metabolites = valid_rows$mz_time_1
    annotated_metabolites = valid_rows$mz_time_2

    # Preallocate vectors for results
    n = length(identified_metabolites)
    cor_values = numeric(n)
    p_values = numeric(n)

    # Loop through each pair of metabolites
    # Here, we just use typical correlation analysis between 2 features. 
    for (i in seq_len(n)) {

        # extract the intensity value from the MSMICA_cor_input for the precursor and product metabolites
        precursor_vector = cor_input[[identified_metabolites[i]]]
        product_vector = cor_input[[annotated_metabolites[i]]]

        # Conduct Spearman correlation
        cor_test_result = cor.test(
            precursor_vector,
            product_vector,
            method = "spearman",
            alternative = "two.sided",
            exact = FALSE
        )
        cor_values[i] = cor_test_result$estimate
        p_values[i] = cor_test_result$p.value
    }

    # Create the final data frame
    cor_MSMICA_data = data.frame(
        mz_time_1 = identified_metabolites,
        mz_time_2 = annotated_metabolites,
        correlation = cor_values,
        p_value = p_values
    )

    cor_MSMICA = tibble(cor_MSMICA_data)

    return(list(cor_MSMICA = cor_MSMICA, MSMICA_col_names_connection_12 = MSMICA_col_names_connection_12))

}