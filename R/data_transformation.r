#' Log2-transform metabolomics intensities
#'
#' Applies a \code{log2} transform to the intensity columns of a
#' metabolomics feature table in wide format. The first two columns are
#' treated as metadata (m/z and retention time) and are renamed to
#' \code{"mz"} and \code{"time"} before transformation. Any \code{NaN}
#' values that arise from taking the log of zero or negative intensities
#' are converted to \code{NA}.
#'
#' @param met_raw_wide A data frame (or tibble) in wide format. The first
#'   column must hold m/z values, the second must hold retention times,
#'   and the remaining columns must hold per-sample intensities.
#' @return The input table with the intensity columns log2-transformed
#'   and the first two columns renamed to \code{mz} and \code{time}.
#' @keywords internal
#' @noRd
data_transformation = function(met_raw_wide){
    # Standardize the metadata column names so downstream code can assume them
    colnames(met_raw_wide)[1] = "mz"
    # update the second column name as time
    colnames(met_raw_wide)[2] = "time"

    met_raw_wide = as.data.frame(met_raw_wide)

    # log2 transformation for all intensity values
    met_raw_wide[,-c(1,2)] = log2(met_raw_wide[,-c(1,2)])

    message("log2 transformation is performed.")

    suppressWarnings({
        # Convert all NaN values to NA
        is.nan.data.frame = function(x)
            do.call(cbind, lapply(x, is.nan))

        met_raw_wide[is.nan(met_raw_wide)] = NA
    })

    return(met_raw_wide)
}