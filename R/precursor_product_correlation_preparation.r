#' Assemble precursor-product feature pairs for correlation scoring
#'
#' Expands each identified/annotated metabolite into the set of candidate
#' precursor-product feature pairs implied by the reaction network
#' (\code{reaction_connection}) and by the candidate features already in
#' the input table. The resulting table is the input to the actual
#' correlation computation step. Self-correlations are dropped by
#' comparing canonical \code{mz_time} strings instead of raw InChIKeys,
#' which is more robust to rounding-induced collisions between close
#' adducts.
#'
#' Relies on \code{reaction_connection} being available in the calling
#' scope (it is defined by \code{MSMICA_algorithm()} depending on the
#' user's \code{reaction_database} choice).
#'
#' @param precursor_product_correlation_preparation_input_data Data
#'   frame of candidate identifications with (at least) columns
#'   \code{InChIKey}, \code{identification_type},
#'   \code{identified_Name}, \code{Adduct}, \code{mz}, and \code{time}.
#' @return A data frame of unique precursor-product candidate pairs
#'   ready for correlation scoring.
#' @keywords internal
#' @noRd
precursor_product_correlation_preparation = function(precursor_product_correlation_preparation_input_data,
                                                     rxn_connection) {
    # Join each identified metabolite to its reactive partners from the network
    MSMICA_col_names_connection = precursor_product_correlation_preparation_input_data %>%
        inner_join(rxn_connection, by = c("InChIKey" = "connection_1_InChIKey"), relationship = "many-to-many")

    # Arrange the data by the mz column in ascending order
    MSMICA_col_names_connection = MSMICA_col_names_connection %>%
        arrange(mz)

    # select only the metabolites with the precursor-product relationships with the metabolites in the identified_metabolite_final
    # MSMICA_col_names_connection_simplifed_2 is the metabolites with relationships with the other metabolites that have precursor-product reactions
    MSMICA_col_names_connection_simplifed = MSMICA_col_names_connection[, c("InChIKey", "identification_type", "identified_Name", "Adduct", "mz", "time",  "connection_2_InChIKey", "react_id", "enzyme_transporter", "source")]

    # remove the duplicates based on "identification_type", "identified_Name", "Adduct", "mz", "time",  "connection_2_InChIKey",
    MSMICA_col_names_connection_simplifed = MSMICA_col_names_connection_simplifed %>%
        distinct(identification_type, identified_Name, Adduct, mz, time, connection_2_InChIKey, .keep_all = TRUE)

    # round mz to 4 decimal places and time to integers
    MSMICA_col_names_connection_simplifed$mz = round(MSMICA_col_names_connection_simplifed$mz, 4)
    MSMICA_col_names_connection_simplifed$time = round(MSMICA_col_names_connection_simplifed$time, 0)
    # add mz_time to the MSMICA_col_names_connection_simplifed
    MSMICA_col_names_connection_simplifed$mz_time = paste0(MSMICA_col_names_connection_simplifed$mz, "_", MSMICA_col_names_connection_simplifed$time)

    # round mz to 4 decimal places and time to integers
    precursor_product_correlation_preparation_input_data$mz = round(precursor_product_correlation_preparation_input_data$mz, 4)
    precursor_product_correlation_preparation_input_data$time = round(precursor_product_correlation_preparation_input_data$time, 0)
    # add mz_time to the precursor_product_correlation_preparation_input_data
    precursor_product_correlation_preparation_input_data$mz_time = paste0(precursor_product_correlation_preparation_input_data$mz, "_", precursor_product_correlation_preparation_input_data$time)

    # inner join the MSMICA_col_names_connection_simplifed with the precursor_product_correlation_preparation_input_data
    MSMICA_col_names_connection_identified_final = inner_join(
        MSMICA_col_names_connection_simplifed, 
        precursor_product_correlation_preparation_input_data, 
        by = c("connection_2_InChIKey" = "InChIKey"),
        suffix = c("_identified", "_identified_final"),
        relationship = "many-to-many"
        )
    
    # remove duplicate rows based on InChIKey and connection_2_InChIKey
    MSMICA_col_names_connection_identified_final = MSMICA_col_names_connection_identified_final %>%
        distinct(InChIKey, connection_2_InChIKey, .keep_all = TRUE)

    # remove those precursor-product correlation values that are self-correlation values (if feature_identified == feature_identified_final)
    ## in the past, the self-correlation values were removed if they are equal to 1, but now, we are removing them if they are the same feature based on mz and time. This is more robust because it tolerates extremely high correlation values (equal to 1 after rounding to 4 decimal places) that may exist between M+H and M+H[+1].
    MSMICA_col_names_connection_identified_final = MSMICA_col_names_connection_identified_final %>%
        mutate(feature_identified = paste0(mz_identified, "_", time_identified)) %>%
        mutate(feature_identified_final = paste0(mz_identified_final, "_", time_identified_final)) %>%
        filter(feature_identified != feature_identified_final) %>%
        dplyr::select(-feature_identified, -feature_identified_final)

    # filter out the rows with NA in InChIKey and connection_2_InChIKey
    MSMICA_col_names_connection_identified_final = MSMICA_col_names_connection_identified_final %>%
        filter(!is.na(InChIKey) & !is.na(connection_2_InChIKey))
    
    # filter out the rows with duplicated InChIKey and connection_2_InChIKey
    MSMICA_col_names_connection_identified_final = MSMICA_col_names_connection_identified_final %>%
        distinct(InChIKey, connection_2_InChIKey, .keep_all = TRUE)

    return(MSMICA_col_names_connection_identified_final)
}