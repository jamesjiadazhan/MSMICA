#' @title Open the custom biochemical reaction dataset 
#' @description
#'
#' `custom_biochemical_reaction_loading` opens the custom biochemical reaction dataset (collected by James Zhan) as a tibble and return it. User can add their own biochemical reaction dataset in the same format and load it with this function.
#' 
#' @return a tibble
#' @export
custom_biochemical_reaction_loading <- function() {
  custom_biochemical_reaction_path <- system.file("biochemical_reaction", "MSMICA_collected_reactions_inchikey_simple.csv", package = "MSMICA")
  custom_biochemical_reaction_path = read_csv('/Users/james/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/Publication/Abstract/MSMICA_package/MSMICA/inst/biochemical_reaction/MSMICA_collected_reactions_inchikey_simple.csv')
  message("Custom biochemical reaction dataset is loading...")
  setwd('/Users/james/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/Publication/Abstract/MSMICA_package/MSMICA/data')
  save(custom_biochemical_reaction_path, file = "custom_biochemical_reaction_path.rda")
  custom_biochemical_reaction = readr::read_csv(custom_biochemical_reaction_path)
  return(custom_biochemical_reaction)
}
