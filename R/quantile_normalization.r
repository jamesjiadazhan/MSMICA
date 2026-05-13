#' Quantile-normalize a metabolomics feature table
#'
#' Applies sample-wise quantile normalization to the intensity columns of
#' a wide-format metabolomics feature table using
#' \code{preprocessCore::normalize.quantiles}. The first two columns are
#' preserved as m/z and retention time metadata. Quantile normalization
#' is used here to remove systematic differences in total signal across
#' samples before downstream correlation-based analyses.
#'
#' @param met_raw_wide A wide-format metabolomics feature table. The
#'   first two columns must be m/z and retention time; the remaining
#'   columns are per-sample intensities.
#' @return A tibble with the same layout as the input, where the
#'   intensity columns have been quantile-normalized.
#' @keywords internal
#' @noRd
quantile_normalization = function(met_raw_wide){
    message("Quantile normalization is performed.")
    ## remove mz and time columns and convert to matrix
    met_raw_wide_quantile_normalized = met_raw_wide %>%
        dplyr::select(-c(mz, time)) %>%
        as.matrix()
    ## save the sample column name
    sample_column_name = colnames(met_raw_wide_quantile_normalized)
    ## perform quantile normalization
    met_raw_wide_quantile_normalized = preprocessCore::normalize.quantiles(met_raw_wide_quantile_normalized)
    ## give the sample column name
    colnames(met_raw_wide_quantile_normalized) = sample_column_name
    ## convert to data frame
    met_raw_wide_quantile_normalized = as.data.frame(met_raw_wide_quantile_normalized)
    ## bring mz and time columns back
    met_raw_wide_quantile_normalized = cbind(met_raw_wide$mz, met_raw_wide$time, met_raw_wide_quantile_normalized)
    ## convert to tibble
    met_raw_wide_quantile_normalized = tibble(met_raw_wide_quantile_normalized)
    ## update the column name
    colnames(met_raw_wide_quantile_normalized) = c("mz", "time", sample_column_name)
    ## update the met_raw_wide using the met_raw_wide_quantile_normalized
    met_raw_wide = met_raw_wide_quantile_normalized
    # return the metabolomics data
    return(met_raw_wide)
}