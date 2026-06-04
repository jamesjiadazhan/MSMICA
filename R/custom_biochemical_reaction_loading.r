#' Load the bundled custom biochemical reaction dataset
#'
#' Reads the curated biochemical reaction table shipped with MSMICA as a
#' tibble for use inside precursor-product correlation preparation.
#'
#' @return A tibble of custom biochemical reactions with simplified InChIKeys.
#' @keywords internal
#' @noRd
custom_biochemical_reaction_loading <- function() {
  custom_biochemical_reaction_path <- system.file("MSMICA_collected_reactions_inchikey_simple.csv", package = "MSMICA")
  message("Custom biochemical reaction dataset is loading...")
  custom_biochemical_reaction = readr::read_csv(custom_biochemical_reaction_path, show_col_types = FALSE)
  return(custom_biochemical_reaction)
}
