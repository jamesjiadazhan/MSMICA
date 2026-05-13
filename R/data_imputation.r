#' Impute missing intensities in a metabolomics feature table
#'
#' Wrapper that dispatches to one of the supported imputation strategies
#' for a wide-format metabolomics feature table. Two strategies are
#' available: \code{"QRILC"}, the Quantile Regression Imputation of
#' Left-Censored data implemented in \pkg{imputeLCMD} (recommended when
#' replicate samples are available), and \code{"half_min"}, which
#' replaces missing values in each feature row with half of that row's
#' observed minimum. Any other value disables imputation.
#'
#' @param met_raw_wide A wide-format metabolomics feature table. Columns
#'   one and two are treated as m/z and retention time metadata; all
#'   remaining columns are per-sample intensities.
#' @param imputation_method Character string. One of \code{"QRILC"},
#'   \code{"half_min"}, or any other value (including \code{NA}) to
#'   disable imputation.
#' @return The input table with intensity columns imputed according to
#'   the chosen strategy.
#' @keywords internal
#' @noRd
data_imputation = function(met_raw_wide, imputation_method){
    # Branch 1: QRILC (quantile regression imputation of left-censored data)
    if (imputation_method == "QRILC"){
        met_raw_wide_meta = as.data.frame(met_raw_wide[, 1:2])
        met_raw_wide_intensity = as.data.frame(met_raw_wide[,-c(1,2)])
        met_raw_wide_intensity_imputed = imputeLCMD::impute.QRILC(met_raw_wide_intensity)
        met_raw_wide_intensity = as.data.frame(met_raw_wide_intensity_imputed[1])
        met_raw_wide = bind_cols(met_raw_wide_meta, met_raw_wide_intensity)

        message("Missing data imputation completed. The imputation method is: QRILC")
        # Remove 'X' from column names if it is the first character
        colnames(met_raw_wide) = gsub("^X", "", colnames(met_raw_wide))
    } else if (imputation_method == "half_min"){

        message("Missing data imputation completed. The imputation method is: half minimum")

        # Function to impute missing values with half of the minimum value in the row
        impute_missing_values = function(values) {
            # Calculate half of the minimum value in the row (excluding NA)
            half_min_value = min(values, na.rm = TRUE) / 2
            # Impute NA values with half of the minimum value
            values[is.na(values)] = half_min_value
            # Combine the first two columns with the imputed values
            return(values)
        }

        # separate the met_raw_wide into meta data and intensity data
        met_raw_wide_meta = met_raw_wide[, 1:2]
        met_raw_wide_intensity = met_raw_wide[, -c(1,2)]

        # Apply the imputation function to each row
        imputed_data = t(apply(met_raw_wide_intensity, 1, impute_missing_values))

        # Convert the result back to a data frame and ensure the column names are preserved
        met_raw_wide_intensity_imputed = as.data.frame(imputed_data, stringsAsFactors = FALSE)
        colnames(met_raw_wide_intensity_imputed) = colnames(met_raw_wide_intensity)
        met_raw_wide_intensity_imputed = tibble(met_raw_wide_intensity_imputed)

        # combine the meta data and the intensity data
        met_raw_wide = bind_cols(met_raw_wide_meta, met_raw_wide_intensity_imputed)
    } else {
        message("No data imputation is performed. All 0 values are replaced with NA.")
    }
    return(met_raw_wide)
}
