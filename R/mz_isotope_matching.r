#' Match feature table m/z values to isotopologue-shifted adduct library
#'
#' Sibling of \code{mz_matching()} that matches observed feature m/z to
#' the isotopologue-shifted theoretical m/z (\code{mz_isotope}) of each
#' metabolite-adduct entry instead of the monoisotopic m/z. This is what
#' allows MSMICA to annotate \code{[M+H]+[+1]}-style isotopic partner
#' peaks, which are then used as evidence during the adduct/isotope
#' clustering step of the Bayesian identification. The \code{Mass_diff}
#' value (e.g. \code{"+1"}, \code{"+2"}) is preserved as a bracketed
#' \code{isotope} tag on each matched row.
#'
#' Expects the globals \code{met_raw_wide_original_mean_intensity},
#' \code{ion_mode}, and \code{mz_matching_folder_name} to be defined in
#' the calling scope (as they are inside \code{MSMICA_algorithm()}).
#'
#' @param met_raw_wide A wide-format metabolomics feature table with
#'   \code{mz} and \code{time} as the first two columns.
#' @param metabolite_mz_library A tibble of metabolite-adduct entries
#'   including the columns \code{mz_isotope}, \code{Mass_diff} and the
#'   usual metadata.
#' @param mz_threshold Numeric. m/z matching tolerance in ppm.
#' @param detail Logical. If \code{TRUE}, writes the cleaned match
#'   table to \code{mz_matching_folder_name} as a CSV.
#' @return A tibble of annotated feature-adduct-isotopologue matches.
#' @keywords internal
#' @noRd
mz_isotope_matching = function(met_raw_wide, metabolite_mz_library, mz_threshold, detail,
                               mean_intensity_ref, ion_mode, output_folder = NULL){
    # Observed m/z values from the feature table
    met_raw_wide_isotope_1 = c(met_raw_wide$mz)
    # Theoretical isotopologue m/z (not monoisotopic) from the adduct library
    metabolite_mz_library_isotope_1 = c(metabolite_mz_library$mz_isotope)
    # make met_raw_wide_isotope_1 a data frame
    met_raw_wide_isotope_1 = as.data.frame(met_raw_wide_isotope_1)
    # make metabolite_mz_library_isotope_1 a data frame
    metabolite_mz_library_isotope_1 = as.data.frame(metabolite_mz_library_isotope_1)
    # find the overlapping mz between met_raw_wide and kegg library using the specified mz tolerance (in ppm)
    masteroverlap.met_raw_wide_kegg_isotope = find.Overlapping.mzs(met_raw_wide_isotope_1, metabolite_mz_library_isotope_1, mz.thresh = mz_threshold)
    # select the matched mz and kegg id columns from the kegg library
    metabolite_mz_library_isotope_2 = slice(metabolite_mz_library, masteroverlap.met_raw_wide_kegg_isotope$index.B)
    # rename the mz_isotope column as mz_annotated_isotope
    metabolite_mz_library_isotope_2 = metabolite_mz_library_isotope_2 %>%
        dplyr::rename(mz_annotated_isotope = mz_isotope) %>%
        ## remove the mz column
        select(-mz)
    # select the matched mz and retention time columns from the met_raw_wide feature table
    met_raw_wide_isotope_2 = slice(met_raw_wide, masteroverlap.met_raw_wide_kegg_isotope$index.A)
    # rename the mz column as mz_sample_isotope and retention time column as time_sample_isotope
    colnames(met_raw_wide_isotope_2)[1] = "mz_sample_isotope"
    colnames(met_raw_wide_isotope_2)[2] = "time_sample_isotope"
    # add mean_intensity column to the met_raw_wide_isotope_2 by merging with met_raw_wide_original_mean_intensity: mz_sample_isotope == mz and time_sample_isotope == time
    met_raw_wide_isotope_3 = met_raw_wide_isotope_2 %>%
        left_join(mean_intensity_ref, by = c("mz_sample_isotope" = "mz", "time_sample_isotope" = "time")) %>%
        ## move the mean_intensity column to the column after the time_sample_isotope column
        relocate(mean_intensity, .after = time_sample_isotope) %>%
        ## round mean_intensity to integers
        mutate(mean_intensity = round(mean_intensity))
    # combine the met_raw_wide feature table and the kegg library
    met_raw_wide_isotope_metabolite_annotated = cbind(metabolite_mz_library_isotope_2, met_raw_wide_isotope_3)
    # convert annotated results from data frame to tibble
    met_raw_wide_isotope_metabolite_annotated_tb = tibble(met_raw_wide_isotope_metabolite_annotated)
    # add the metabolomics column (HILIC or C18) and mz_time columns to each of them
    met_raw_wide_isotope_metabolite_annotated_tb_1 = met_raw_wide_isotope_metabolite_annotated_tb %>%
        mutate(
            ion_mode = ion_mode,
            mz_sample_isotope = round(mz_sample_isotope, 4),
            time_sample_isotope = round(time_sample_isotope, 0),
            mz_time_sample_isotope = paste0(mz_sample_isotope, "_", time_sample_isotope)
            ) %>%
        # reorder the columns
        dplyr::select(ion_mode, KEGG_ID, mz_annotated_isotope, Name, Formula, Most_abundant_isotopologue_mass, Adduct, mz_sample_isotope, time_sample_isotope, mz_time_sample_isotope, everything()) %>%
        arrange(KEGG_ID)
    # create a new column, isotope, and set it to the pasted combiation of [, Mass_diff, and ]: this is for labeling the mass difference for the isotopic adducts
    met_raw_wide_isotope_metabolite_annotated_tb_1$isotope = paste0("[", met_raw_wide_isotope_metabolite_annotated_tb_1$Mass_diff, "]")
    # remove the Mass_diff column
    met_raw_wide_isotope_metabolite_annotated_tb_1 = met_raw_wide_isotope_metabolite_annotated_tb_1 %>%
        dplyr::select(-c(Mass_diff)) %>%
        ## relocate the isotope column right after the Adduct column
        relocate(isotope, .after = Adduct)
    # round Mono_mass to 4 decimal places
    met_raw_wide_isotope_metabolite_annotated_tb_1$Mono_mass = round(met_raw_wide_isotope_metabolite_annotated_tb_1$Mono_mass, 4)
    # remove duplicates
    met_raw_wide_isotope_metabolite_annotated_tb_1 = distinct(met_raw_wide_isotope_metabolite_annotated_tb_1)
    # if detail is TRUE, save the data in the mz_matching_folder_name folder
    if (detail == TRUE){
        # keep only the first 28 columns
        met_raw_wide_isotope_metabolite_annotated_tb_1_output = met_raw_wide_isotope_metabolite_annotated_tb_1[, 1:28]
        # rename mz_annotated_isotope as theoretical_mz, mz_sample_isotope as mz, time_sample_isotope as time, and mz_time_sample_isotope as mz_time
        met_raw_wide_isotope_metabolite_annotated_tb_1_output = met_raw_wide_isotope_metabolite_annotated_tb_1_output %>%
            dplyr::rename(theoretical_mz = mz_annotated_isotope, mz = mz_sample_isotope, time = time_sample_isotope, mz_time = mz_time_sample_isotope)
        # add mz_matching_ppm between theoretical_mz and mz
        met_raw_wide_isotope_metabolite_annotated_tb_1_output$mz_matching_ppm = abs(((met_raw_wide_isotope_metabolite_annotated_tb_1_output$mz - met_raw_wide_isotope_metabolite_annotated_tb_1_output$theoretical_mz)/met_raw_wide_isotope_metabolite_annotated_tb_1_output$theoretical_mz)*1000000)
        # relocate the mz_matching_ppm right after mz_time column, and relocate theoretical_mz before mz_matching_ppm column
        met_raw_wide_isotope_metabolite_annotated_tb_1_output = met_raw_wide_isotope_metabolite_annotated_tb_1_output %>%
            relocate(mz_matching_ppm, .after = mz_time) %>%
            relocate(theoretical_mz, .before = mz_matching_ppm)
        # round theoretical_mz and mz to 4 decimal places
        met_raw_wide_isotope_metabolite_annotated_tb_1_output$theoretical_mz = round(met_raw_wide_isotope_metabolite_annotated_tb_1_output$theoretical_mz, 4)
        met_raw_wide_isotope_metabolite_annotated_tb_1_output$mz = round(met_raw_wide_isotope_metabolite_annotated_tb_1_output$mz, 4)
        # round time to 0 decimal place and mz_matching_ppm to 1 decimal place
        met_raw_wide_isotope_metabolite_annotated_tb_1_output$time = round(met_raw_wide_isotope_metabolite_annotated_tb_1_output$time, 0)
        met_raw_wide_isotope_metabolite_annotated_tb_1_output$mz_matching_ppm = round(met_raw_wide_isotope_metabolite_annotated_tb_1_output$mz_matching_ppm, 1)
        # round Most_abundant_isotopologue_mass to 4 decimal place
        met_raw_wide_isotope_metabolite_annotated_tb_1_output$Most_abundant_isotopologue_mass = round(met_raw_wide_isotope_metabolite_annotated_tb_1_output$Most_abundant_isotopologue_mass, 4)
        # clean the Name by selecting only the values before the first ":"
        met_raw_wide_isotope_metabolite_annotated_tb_1_output = met_raw_wide_isotope_metabolite_annotated_tb_1_output %>%
            ## select only the values in the Name column before the first ";"
            mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1))
        write_csv(met_raw_wide_isotope_metabolite_annotated_tb_1_output, paste0(output_folder, "/", ion_mode, "_metabolite_annotated_", mz_threshold, "ppm_mzonly_isotope.csv"), progress = FALSE)
        rm(met_raw_wide_isotope_metabolite_annotated_tb_1_output)
    }
    
    # return the data
    return(met_raw_wide_isotope_metabolite_annotated_tb_1)
}