#' Annotate the \code{identification_method} column with supporting evidence
#'
#' Appends a human-readable tag to each row's \code{identification_method}
#' string describing which of the three downstream pieces of evidence
#' actually contributed to that identification. Tags are only appended
#' once per category, so repeated calls are idempotent. A sentinel
#' \code{Concentration_average} value of \code{1e-6} (set upstream to
#' represent "no biospecimen-specific concentration available") is
#' replaced with \code{NA} so it does not get tagged.
#'
#' @param df Data frame of MSMICA results containing (at least) the
#'   columns \code{identification_method}, \code{time_difference},
#'   \code{correlation} and \code{Concentration_average}.
#' @return The input data frame with its \code{identification_method}
#'   column augmented in place.
#' @keywords internal
#' @noRd
process_identification_method = function(df) {
    # Tag: retention time prediction (anywhere RT prediction gave a time_difference)
    df$identification_method = ifelse(
        !is.na(df$time_difference) & !grepl("retention time prediction", df$identification_method),
        paste0(df$identification_method, "; retention time prediction"),
        df$identification_method
    )

    # Tag: precursor-product/transporter correlation (biological network evidence)
    df$identification_method = ifelse(
        !is.na(df$correlation) & !grepl("precursor-product/transporter correlation", df$identification_method),
        paste0(df$identification_method, "; precursor-product/transporter correlation"),
        df$identification_method
    )
    # Sentinel 1e-6 was used upstream to mean "no HMDB concentration"; restore to NA
    df$Concentration_average[df$Concentration_average == 1e-6] = NA
    # Tag: biospecimen-specific concentration prior was used in scoring
    df$identification_method = ifelse(
        !is.na(df$Concentration_average) & !grepl("biospecimen-specific concentration", df$identification_method),
        paste0(df$identification_method, "; biospecimen-specific concentration"),
        df$identification_method
    )
    return(df)
}