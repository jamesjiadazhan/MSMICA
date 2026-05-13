#' target_metabolite_search
#' 
#' This function searches for target metabolites in a reference library and a feature table.
#' @param reference_library A reference library table. The reference mz should be called "reference_mz" and the reference retention time should be called "reference_time".
#' @param feature_table A feature table. The mz should be called "mz" and the time should be called "time".
#' @param use_retention_time A logical value indicating whether to use retention time for matching. Default is TRUE.
#' @param mz_threshold The m/z threshold for matching features between two data sets. Default is 10 ppm.
#' @param time_threshold The retention time threshold for matching features between two data sets. Default is 30 seconds. If you use minutes, please convert the unit appropriately.
#' @return A data frame with the overlapping features, including
#'   `mz_ppm_difference` and, when retention-time matching is used,
#'   `time_difference`.
#' @export target_metabolite_search
target_metabolite_search = function(reference_library, feature_table, use_retention_time = TRUE, mz_threshold = 10, time_threshold = 30){
    library(data.table)

    # define the function to find overlapping mz and time between two data sets
    find_overlapping_mz_time <- function(dataA, dataB, mz.thresh = 5, time.thresh = NA) {
        data_a <- data.table::as.data.table(dataA)
        data_b <- data.table::as.data.table(dataB)

        data.table::setnames(data_a, old = names(data_a)[1], new = "mz_data_A")
        data.table::setnames(data_b, old = names(data_b)[1], new = "mz.data_B")
        data.table::set(data_b, j = "mz.data_B_actual", value = data_b[["mz.data_B"]])

        use_time <- !is.na(time.thresh)

        if (use_time) {
            data.table::setnames(data_a, old = names(data_a)[2], new = "time_data_A")
            data.table::setnames(data_b, old = names(data_b)[2], new = "time_data_B")
            data.table::set(data_b, j = "time_data_B_actual", value = data_b[["time_data_B"]])
            message("Using the 1st column as 'mz' and 2nd column as 'retention time'")
        } else {
            message("Using the 1st column as 'mz'")
        }

        data.table::set(data_a, j = "index_A", value = seq_len(nrow(data_a)))
        data.table::set(data_b, j = "index_B", value = seq_len(nrow(data_b)))
        data.table::set(data_a, j = "mz_tol", value = mz.thresh * data_a[["mz_data_A"]] / 1e6)
        data.table::set(data_a, j = "mz_lower", value = data_a[["mz_data_A"]] - data_a[["mz_tol"]])
        data.table::set(data_a, j = "mz_upper", value = data_a[["mz_data_A"]] + data_a[["mz_tol"]])

        matches <- data_b[
            data_a,
            on = c("mz.data_B>=mz_lower", "mz.data_B<=mz_upper"),
            allow.cartesian = TRUE,
            nomatch = 0
        ]

        if (nrow(matches) == 0) {
            if (use_time) {
                return(
                    base::data.frame(
                        index_A = integer(),
                        mz_data_A = numeric(),
                        time_data_A = numeric(),
                        index_B = integer(),
                        mz.data_B = numeric(),
                        time_data_B = numeric(),
                        mz_difference_ppm = numeric(),
                        time_difference_sec = numeric(),
                        check.names = FALSE
                    )
                )
            }

            return(
                base::data.frame(
                    index_A = integer(),
                    mz_data_A = numeric(),
                    index_B = integer(),
                    mz.data_B = numeric(),
                    mz_difference_ppm = numeric(),
                    check.names = FALSE
                )
            )
        }

        data.table::set(
            matches,
            j = "mz_difference_ppm",
            value = (abs(matches[["mz_data_A"]] - matches[["mz.data_B_actual"]]) / matches[["mz_data_A"]]) * 1e6
        )
        data.table::set(matches, j = "mz.data_B", value = matches[["mz.data_B_actual"]])

        if (use_time) {
            data.table::set(
                matches,
                j = "time_difference_sec",
                value = abs(matches[["time_data_A"]] - matches[["time_data_B_actual"]])
            )
            data.table::set(matches, j = "time_data_B", value = matches[["time_data_B_actual"]])
            matches <- matches[matches[["time_difference_sec"]] <= time.thresh]

            if (nrow(matches) == 0) {
                return(
                    base::data.frame(
                        index_A = integer(),
                        mz_data_A = numeric(),
                        time_data_A = numeric(),
                        index_B = integer(),
                        mz.data_B = numeric(),
                        time_data_B = numeric(),
                        mz_difference_ppm = numeric(),
                        time_difference_sec = numeric(),
                        check.names = FALSE
                    )
                )
            }

            selected_columns <- c(
                "index_A",
                "mz_data_A",
                "time_data_A",
                "index_B",
                "mz.data_B",
                "time_data_B",
                "mz_difference_ppm",
                "time_difference_sec"
            )

            return(
                matches[, selected_columns, with = FALSE] |>
                    base::as.data.frame(check.names = FALSE)
            )
        }

        selected_columns <- c("index_A", "mz_data_A", "index_B", "mz.data_B", "mz_difference_ppm")

        matches[, selected_columns, with = FALSE] |>
            base::as.data.frame(check.names = FALSE)
    }

    if (use_retention_time == TRUE){
        print("Using m/z and retention time for matching features between the reference library and the feature table.")
        print(paste("m/z threshold:", mz_threshold, "ppm"))
        print(paste("Retention time threshold:", time_threshold, "seconds"))
        # select only mz and time columns
        feature_table_1 = feature_table[,c("mz", "time")]
        # select only reference_mz and reference_time columns
        reference_library_1 = reference_library[,c("reference_mz", "reference_time")]
        # find overlapping mz within 5 ppm and time within 30 seconds
        masteroverlap.raw_reference_library = find_overlapping_mz_time(reference_library_1, feature_table_1, mz.thresh = mz_threshold, time.thresh = time_threshold)
    }
    else{
        print("Using m/z for matching features between the reference library and the feature table.")
        print(paste("m/z threshold:", mz_threshold, "ppm"))
        print(paste("Retention time threshold:", time_threshold, "seconds"))
        # select only mz columns
        feature_table_1 = feature_table[,c("mz")]
        # select only reference_mz columns
        reference_library_1 = reference_library[,c("reference_mz")]
        # find overlapping mz within 5 ppm
        masteroverlap.raw_reference_library = find_overlapping_mz_time(reference_library_1, feature_table_1, mz.thresh = mz_threshold)
    }

    if (nrow(masteroverlap.raw_reference_library) == 0) {
        return(reference_library[0, , drop = FALSE])
    }

    # select the overlapping rows from reference_library
    reference_library_2 = dplyr::slice(reference_library, masteroverlap.raw_reference_library$index_A)

    # select the overlapping rows from feature_table
    feature_table_2 = dplyr::slice(feature_table, masteroverlap.raw_reference_library$index_B)

    match_metrics = data.frame(
        mz_ppm_difference = masteroverlap.raw_reference_library$mz_difference_ppm
    )

    if ("time_difference_sec" %in% colnames(masteroverlap.raw_reference_library)) {
        match_metrics$time_difference = masteroverlap.raw_reference_library$time_difference_sec
    } else {
        match_metrics$time_difference = NA_real_
    }

    # combine the two data frames by columns: now we have the metabolite, KEGGID, adduct, mz (sample), time (sample), and sample intensity columns
    target_metabolite_search = cbind(reference_library_2, match_metrics, feature_table_2)

    # return the data as tibble
    target_metabolite_search = tibble::as_tibble(target_metabolite_search)

    # return the data
    return(target_metabolite_search)
}
