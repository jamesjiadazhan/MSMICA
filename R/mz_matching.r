#' Match feature table m/z values to an adduct library
#'
#' Finds all rows of a metabolomics feature table whose m/z values fall
#' within \code{mz_threshold} (ppm) of a theoretical adduct m/z in a
#' metabolite library. Matching is delegated to
#' \code{find.Overlapping.mzs()}. The result is an annotated tibble
#' where each row pairs an observed feature (its sample m/z, retention
#' time, and mean intensity) with a candidate metabolite-adduct entry
#' (KEGG_ID, Name, Formula, Mono_mass, Adduct). Optionally writes a
#' CSV of the match to disk when \code{detail = TRUE}.
#'
#' Downstream MSMICA steps use these m/z matches as the starting
#' candidate set for Bayesian identification. Expects the globals
#' \code{met_raw_wide_original_mean_intensity}, \code{ion_mode}, and
#' \code{mz_matching_folder_name} to be defined in the calling scope
#' (as they are inside \code{MSMICA_algorithm()}).
#'
#' @param met_raw_wide A wide-format metabolomics feature table with
#'   \code{mz} and \code{time} as the first two columns.
#' @param metabolite_mz_library A tibble of metabolite-adduct entries
#'   including \code{mz}, \code{mz_isotope}, \code{KEGG_ID},
#'   \code{Name}, \code{Formula}, \code{Mono_mass},
#'   \code{Most_abundant_isotopologue_mass}, \code{Adduct}, and
#'   \code{Mass_diff}.
#' @param mz_threshold Numeric. m/z matching tolerance in ppm.
#' @param detail Logical. If \code{TRUE}, writes the cleaned match
#'   table to \code{mz_matching_folder_name} as a CSV.
#' @return A tibble of annotated feature-adduct matches.
#' @keywords internal
#' @noRd
mz_matching = function(met_raw_wide, metabolite_mz_library, mz_threshold, detail,
                       mean_intensity_ref, ion_mode, output_folder = NULL){
    # Pull the observed and theoretical m/z vectors as single-column data
    # frames so find.Overlapping.mzs() can return positional indices.
    met_raw_wide_1 = c(met_raw_wide$mz)
    # select the mz column
    metabolite_mz_library_1 = c(metabolite_mz_library$mz)
    # make met_raw_wide_1 a data frame
    met_raw_wide_1 = as.data.frame(met_raw_wide_1)
    # make metabolite_mz_library_1 a data frame
    metabolite_mz_library_1 = as.data.frame(metabolite_mz_library_1)
    # find the overlapping mz between met_raw_wide and metabolite database using the specified mz tolerance (in ppm)
    masteroverlap.met_raw_wide_kegg = find.Overlapping.mzs(met_raw_wide_1, metabolite_mz_library_1, mz.thresh = mz_threshold)
    # select the matched mz and kegg id columns from the  metabolite database
    metabolite_mz_library_2 = slice(metabolite_mz_library, masteroverlap.met_raw_wide_kegg$index.B)
    # rename the mz column as mz_annotated
    metabolite_mz_library_2 = metabolite_mz_library_2 %>%
        dplyr::rename(mz_annotated = mz)
    # select the matched mz and retention time columns from the met_raw_wide feature table
    met_raw_wide_2 = slice(met_raw_wide, masteroverlap.met_raw_wide_kegg$index.A)
    # rename the mz column as mz_sample and retention time column as time_sample
    colnames(met_raw_wide_2)[1] = "mz_sample"
    colnames(met_raw_wide_2)[2] = "time_sample"
    # add mean_intensity column to the met_raw_wide_2 by merging with met_raw_wide_original_mean_intensity: mz_sample == mz and time_sample == time
    met_raw_wide_3 = met_raw_wide_2 %>%
        left_join(mean_intensity_ref, by = c("mz_sample" = "mz", "time_sample" = "time")) %>%
        ## move the mean_intensity column to the column after the time_sample column
        relocate(mean_intensity, .after = time_sample) %>%
        ## round mean_intensity to integers
        mutate(mean_intensity = round(mean_intensity))
    # combine the met_raw_wide feature table and the kegg library
    met_raw_wide_metabolite_annotated = cbind(metabolite_mz_library_2, met_raw_wide_3)
    # convert annotated results from data frame to tibble
    met_raw_wide_metabolite_annotated_tb = tibble(met_raw_wide_metabolite_annotated)
    # add the metabolomics column (HILIC or C18) and mz_time columns to each of them
    met_raw_wide_metabolite_annotated_tb_1 = met_raw_wide_metabolite_annotated_tb %>%
        mutate(
            ion_mode = ion_mode,
            mz_sample = round(mz_sample, 4),
            time_sample = round(time_sample, 0),
            mz_time_sample = paste0(mz_sample, "_", time_sample)
            ) %>%
        # remove Most_abundant_isotopologue_mass, mz_isotope
        dplyr::select(-c(Most_abundant_isotopologue_mass, mz_isotope)) %>%
        # reorder the columns
        dplyr::select(ion_mode, KEGG_ID, mz_annotated, Name, Formula, Mono_mass, Adduct, mz_sample, time_sample, mz_time_sample, everything()) %>%
        arrange(KEGG_ID)
    # create a new column, isotope, and set it to NA since this is not the isotope data
    met_raw_wide_metabolite_annotated_tb_1$isotope = NA
    # remove the Mass_diff column
    met_raw_wide_metabolite_annotated_tb_1 = met_raw_wide_metabolite_annotated_tb_1 %>%
        dplyr::select(-c(Mass_diff)) %>%
        ## relocate the isotope column right after the Adduct column
        relocate(isotope, .after = Adduct)
    # round Mono_mass to 4 decimal places
    met_raw_wide_metabolite_annotated_tb_1$Mono_mass = round(met_raw_wide_metabolite_annotated_tb_1$Mono_mass, 4)
    # remove duplicates
    met_raw_wide_metabolite_annotated_tb_1 = distinct(met_raw_wide_metabolite_annotated_tb_1)
    # if detail is TRUE, save the data in the mz_matching_folder_name folder
    if (detail == TRUE){
        # keep only the first 27 columns
        met_raw_wide_metabolite_annotated_tb_1_output = met_raw_wide_metabolite_annotated_tb_1[, 1:27]
        # rename mz_annotated as theoretical_mz, mz_sample as mz, time_sample as time, and mz_time_sample as mz_time
        met_raw_wide_metabolite_annotated_tb_1_output = met_raw_wide_metabolite_annotated_tb_1_output %>%
            dplyr::rename(theoretical_mz = mz_annotated, mz = mz_sample, time = time_sample, mz_time = mz_time_sample)
        # add mz_matching_ppm between theoretical_mz and mz
        met_raw_wide_metabolite_annotated_tb_1_output$mz_matching_ppm = abs(((met_raw_wide_metabolite_annotated_tb_1_output$mz - met_raw_wide_metabolite_annotated_tb_1_output$theoretical_mz)/met_raw_wide_metabolite_annotated_tb_1_output$theoretical_mz)*1000000)
        # relocate the mz_matching_ppm right after mz_time column, and relocate theoretical_mz before mz_matching_ppm column
        met_raw_wide_metabolite_annotated_tb_1_output = met_raw_wide_metabolite_annotated_tb_1_output %>%
            relocate(mz_matching_ppm, .after = mz_time) %>%
            relocate(theoretical_mz, .before = mz_matching_ppm)
        # round theoretical_mz and mz to 4 decimal places
        met_raw_wide_metabolite_annotated_tb_1_output$theoretical_mz = round(met_raw_wide_metabolite_annotated_tb_1_output$theoretical_mz, 4)
        met_raw_wide_metabolite_annotated_tb_1_output$mz = round(met_raw_wide_metabolite_annotated_tb_1_output$mz, 4)
        # round time to 0 decimal place and mz_matching_ppm to 1 decimal place
        met_raw_wide_metabolite_annotated_tb_1_output$time = round(met_raw_wide_metabolite_annotated_tb_1_output$time, 0)
        met_raw_wide_metabolite_annotated_tb_1_output$mz_matching_ppm = round(met_raw_wide_metabolite_annotated_tb_1_output$mz_matching_ppm, 1)
        # clean the Name by selecting only the values before the first ":"
        met_raw_wide_metabolite_annotated_tb_1_output = met_raw_wide_metabolite_annotated_tb_1_output %>%
            ## select only the values in the Name column before the first ";"
            mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1))
        # save the file
        write_csv(met_raw_wide_metabolite_annotated_tb_1_output, paste0(output_folder, "/", ion_mode, "_metabolite_annotated_", mz_threshold, "ppm_mzonly.csv"), progress = FALSE)
        rm(met_raw_wide_metabolite_annotated_tb_1_output)
    }

    # return the data
    return(met_raw_wide_metabolite_annotated_tb_1)
}
