#' Get MSMICA adduct presets
#'
#' Returns beginner-friendly adduct presets for common LC-MS ion mode and
#' sample-type combinations. The returned vector can be passed directly to
#' the \code{All_Adduct} argument of \code{MSMICA_algorithm()}.
#'
#' @param mode Character. Ionization mode, either \code{"positive"} or
#'   \code{"negative"}.
#' @param sample_type Character. Sample type preset, either \code{"fluid"} or
#'   \code{"tissue"}.
#' @return A character vector of adduct names.
#' @examples
#' msmica_adducts(mode = "positive", sample_type = "fluid")
#' msmica_adducts(mode = "negative", sample_type = "tissue")
#' @export
msmica_adducts = function(mode = c("positive", "negative"),
                          sample_type = c("fluid", "tissue")) {
    mode = match.arg(tolower(mode), c("positive", "negative"))
    sample_type = match.arg(tolower(sample_type), c("fluid", "tissue"))

    preset = paste(mode, sample_type, sep = "_")

    switch(
        preset,
        positive_fluid = c(
            "M+H", "M+Na", "M+2Na-H", "M+H-H2O", "M+H-NH3",
            "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O"
        ),
        negative_fluid = c(
            "M-H", "M+Cl", "M+FA-H", "M+Hac-H", "M-H+HCOONa",
            "M+Na-2H", "M-2H", "2M-H", "M+ACN-H"
        ),
        positive_tissue = c(
            "M+H", "M+K", "M+2K-H", "M+H-H2O", "M+H-NH3",
            "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O"
        ),
        negative_tissue = c(
            "M-H", "M+Cl", "M+FA-H", "M+Hac-H", "M-H+HCOOK",
            "M-2H", "2M-H", "M+ACN-H", "M+K-2H"
        )
    )
}
