#' MSMICA algorithm
#' 
#' This function is used to perform the MSMICA algorithm for metabolite identification using the metabolomics feature table and the KEGG database. The function takes the following inputs:
#' @import dplyr
#' @import readr
#' @import tidyr
#' @import progress
#' @param met_raw_wide a metabolomics feature table in wide format with mz as the first column, time as the second column, and intensity values as the remaining columns.
#' @param metabolite_identified a data frame containing the metabolites that have been identified in the metabolomics feature table. The first column is ion_mode (positive or negative). The second column is metabolite (the name of metabolite). The third column is KEGGID (if no KEGGID, use NA). The fourth column is adduct (the adduct form of a metabolite). The fifth column is mz and the sixth column is time. The rest of the columns are the intensity values. Default is NULL.
#' @param identified_target: a vector containing the KEGGID that are used as additional starting points for the MSMICA algorithm. This can be used when users just want to use specific identified metabolites as the starting point for the MSMICA algorithm. Default is NULL.
#' @param Adduct a vector containing the adduct forms of the metabolites. Default is c("M+H","M+2Na-H","M+Na","M-H2O+H","M+K","M+2H") for the positive mode.
#' @param prefix a prefix to be added to the output files. Default is "".
#' @param ion_mode the ionization mode of the metabolomics data. Default is "positive". Other options are "negative".
#' @param detail a logical value indicating whether to save the intermediate results as csv files. Default is FALSE. WARNING, this can create thousands of files with a lot of space. Use with caution.
#' @param adduct_correlation_r_threshold the correlation threshold for adduct correlation analysis. Default is 0.4 (spearman correlation).
#' @param adduct_correlation_time_threshold the retention time threshold for adduct correlation analysis. Default is 5 (seconds).
#' @param isotopic_correlation_r_threshold the correlation threshold for isotopic correlation analysis. Default is 0.5 (spearman correlation).
#' @param isotopic_correlation_time_threshold the retention time threshold for isotopic correlation analysis. Default is 5 (seconds).
#' @param mean_absolute_isotope_ratio_deviation_threshold the threshold for the mean absolute isotope ratio deviation between predicted isotopic relative abundance and actual isotopic relative abundance compared to the primary adduct without heavy isotope. Default is 5 (percent).
#' @export MSMICA_algorithm

MSMICA_algorithm <- function(met_raw_wide, metabolite_identified = NULL, identified_target = NULL, Adduct = c("M+H","M+2Na-H","M+Na","M-H2O+H","M+K","M+2H"), prefix="", ion_mode="positive", detail=FALSE, adduct_correlation_r_threshold=0.4, adduct_correlation_time_threshold=5, isotopic_correlation_r_threshold=0.5, isotopic_correlation_time_threshold=5, mean_absolute_isotope_ratio_deviation_threshold=5) {
    # import the KEGG compound database with predicted isotopic mass
    data(KEGG_database)

    # import the precursor-product enzyme-based reaction KEGG database
    data(kegg_connection)

    # Cross join Adduct with KEGG_database to form a adduct list for each unique KEGGID 
    KEGG_database_2 <- tidyr::crossing(KEGG_database, Adduct)

    # calculate mz using Exact_mass and Adduct
    KEGG_database_2_mz = MetaboCoreUtilsAdduct::mass2mz_df(mass=KEGG_database_2$Exact_mass, adduct=KEGG_database_2$Adduct)
    # add mz to the KEGG_database_2
    KEGG_database_2$mz = KEGG_database_2_mz$mz
    # remove KEGG_database_2_mz to save memory
    rm(KEGG_database_2_mz)

    # calculate the heavy isotope mz using isotopic_mass and and Adduct
    KEGG_database_2_mz_heavy = MetaboCoreUtilsAdduct::mass2mz_df(mass=KEGG_database_2$isotopic_mass, adduct=KEGG_database_2$Adduct)
    # add heavy isotope mz to the KEGG_database_2
    KEGG_database_2$mz_isotope = KEGG_database_2_mz_heavy$mz
    # remove KEGG_database_2_mz_heavy to save memory
    rm(KEGG_database_2_mz_heavy)

    # remove mz less than 0
    kegg_library = KEGG_database_2 %>%
        filter(mz > 0)
    
    # update the first column name as mz
    colnames(met_raw_wide)[1] = "mz"
    # update the second column name as time
    colnames(met_raw_wide)[2] = "time"

    ########### mz matching for all primary and secondary adducts
    # select the mz column
    met_raw_wide_1 = met_raw_wide[,"mz"]
    # select the mz column
    kegg_library_1 = kegg_library[,"mz"]
    # find the overlapping mz between met_raw_wide and kegg library using 5 ppm mz tolerance
    masteroverlap.met_raw_wide_kegg = find.Overlapping.mzs(met_raw_wide_1, kegg_library_1, mz.thresh = 5, time.thresh = NA, alignment.tool = NA)
    # select the matched mz and kegg id columns from the kegg library
    kegg_library_2 = slice(kegg_library, masteroverlap.met_raw_wide_kegg$index.B)
    # rename the mz column as mz_annotated
    colnames(kegg_library_2)[9] = "mz_annotated"
    # dplyr::select the matched mz and retention time columns from the met_raw_wide feature table
    met_raw_wide_2 = slice(met_raw_wide, masteroverlap.met_raw_wide_kegg$index.A)
    # rename the mz column as mz_sample and retention time column as time_sample
    colnames(met_raw_wide_2)[1] = "mz_sample"
    colnames(met_raw_wide_2)[2] = "time_sample"
    met_raw_wide_3 = met_raw_wide_2
    # calculate the average itensity of each feature
    met_raw_wide_4 = met_raw_wide_3 %>%
        group_by(mz_sample, time_sample) %>%
        # calculate the mean of the intensity by doing the average for each row
        mutate(
            intensity_mean = rowMeans(across(where(is.numeric)))
        ) %>%
        # move the intensity_mean column to the column after the time_sample column
        relocate(intensity_mean, .after = time_sample) %>%
        ungroup()
    # combine the met_raw_wide feature table and the kegg library
    met_raw_wide_kegg_annotated = cbind(kegg_library_2, met_raw_wide_4)
    # convert annotated results from data frame to tibble
    met_raw_wide_kegg_annotated_tb = tibble(met_raw_wide_kegg_annotated)
    # add the metabolomics column (HILIC or C18) and mz_time columns to each of them
    met_raw_wide_kegg_annotated_tb_1 = met_raw_wide_kegg_annotated_tb %>%
        mutate(
            ion_mode = ion_mode,
            mz_time_sample = paste0(round(mz_sample, 4), "_", round(time_sample, 1))
            ) %>%
        # remove isotopic_mass, mz_isotope
        dplyr::select(-c(isotopic_mass, mz_isotope)) %>%
        # reorder the columns
        dplyr::select(ion_mode, KEGGID, mz_annotated, Name, Formula, Exact_mass, Adduct, mz_sample, time_sample, mz_time_sample, everything()) %>%
        arrange(KEGGID)
    # create a new column, isotope, and set it to NA since this is not the isotope data
    met_raw_wide_kegg_annotated_tb_1$isotope = NA
    # remove the mass_diff column
    met_raw_wide_kegg_annotated_tb_1 = met_raw_wide_kegg_annotated_tb_1 %>%
        dplyr::select(-c(mass_diff)) %>%
        ## relocate the isotope column right after the Adduct column
        relocate(isotope, .after = Adduct)
    # remove duplicates
    met_raw_wide_kegg_annotated_tb_1 = distinct(met_raw_wide_kegg_annotated_tb_1)
    # save the data if detail = TRUE
    if (detail == TRUE) {
        write_csv(met_raw_wide_kegg_annotated_tb_1, paste0(ion_mode, "_kegg_annotated_5ppm_mzonly.csv"))
    }

    ########### mz matching for all primary and secondary adducts
    # select the mz column
    met_raw_wide_isotope_1 = met_raw_wide[,"mz"]
    # select the mz column
    kegg_library_isotope_1 = kegg_library[,"mz_isotope"]
    # find the overlapping mz between met_raw_wide and kegg library using 5 ppm mz tolerance
    masteroverlap.met_raw_wide_kegg_isotope = find.Overlapping.mzs(met_raw_wide_isotope_1, kegg_library_isotope_1, mz.thresh = 5, time.thresh = NA, alignment.tool = NA)
    # select the matched mz and kegg id columns from the kegg library
    kegg_library_isotope_2 = slice(kegg_library, masteroverlap.met_raw_wide_kegg_isotope$index.B)
    # rename the mz column as mz_annotated
    colnames(kegg_library_isotope_2)[10] = "mz_annotated_isotope"
    # dplyr::select the matched mz and retention time columns from the met_raw_wide feature table
    met_raw_wide_isotope_2 = slice(met_raw_wide, masteroverlap.met_raw_wide_kegg_isotope$index.A)
    # rename the mz column as mz_sample_isotope and retention time column as time_sample_isotope
    colnames(met_raw_wide_isotope_2)[1] = "mz_sample_isotope"
    colnames(met_raw_wide_isotope_2)[2] = "time_sample_isotope"
    met_raw_wide_isotope_3 = met_raw_wide_isotope_2
    # calculate the average itensity of each feature
    met_raw_wide_isotope_4 = met_raw_wide_isotope_3 %>%
        group_by(mz_sample_isotope, time_sample_isotope) %>%
        # calculate the mean of the intensity by doing the average for each row
        mutate(
            intensity_mean = rowMeans(across(where(is.numeric)))
        ) %>%
        # move the intensity_mean column to the column after the time_sample_isotope column
        relocate(intensity_mean, .after = time_sample_isotope) %>%
        ungroup()
    # combine the met_raw_wide feature table and the kegg library
    met_raw_wide_isotope_kegg_annotated = cbind(kegg_library_isotope_2, met_raw_wide_isotope_4)
    # convert annotated results from data frame to tibble
    met_raw_wide_isotope_kegg_annotated_tb = tibble(met_raw_wide_isotope_kegg_annotated)
    # add the metabolomics column (HILIC or C18) and mz_time columns to each of them
    met_raw_wide_isotope_kegg_annotated_tb_1 = met_raw_wide_isotope_kegg_annotated_tb %>%
        mutate(
            ion_mode = ion_mode,
            mz_time_sample_isotope = paste0(round(mz_sample_isotope, 4), "_", round(time_sample_isotope, 1))
            ) %>%
        # remove Exact_mass, mz
        dplyr::select(-c(Exact_mass, mz)) %>%
        # reorder the columns
        dplyr::select(ion_mode, KEGGID, mz_annotated_isotope, Name, Formula, isotopic_mass, Adduct, mz_sample_isotope, time_sample_isotope, mz_time_sample_isotope, everything()) %>%
        arrange(KEGGID)
    # create a new column, isotope, and set it to the pasted combiation of [, mass_diff, and ]
    met_raw_wide_isotope_kegg_annotated_tb_1$isotope = paste0("[", met_raw_wide_isotope_kegg_annotated_tb_1$mass_diff, "]")
    # remove the mass_diff column
    met_raw_wide_isotope_kegg_annotated_tb_1 = met_raw_wide_isotope_kegg_annotated_tb_1 %>%
        dplyr::select(-c(mass_diff)) %>%
        ## relocate the isotope column right after the Adduct column
        relocate(isotope, .after = Adduct)
    # remove duplicates
    met_raw_wide_isotope_kegg_annotated_tb_1 = distinct(met_raw_wide_isotope_kegg_annotated_tb_1)
    # save the data if detail = TRUE
    if (detail == TRUE) {
        write_csv(met_raw_wide_isotope_kegg_annotated_tb_1, paste0(ion_mode, "_kegg_annotated_5ppm_mzonly_isotope.csv"))
    }

    # Update the met_raw_wide as the met_raw_wide_kegg_annotated_tb_1 for later use
    met_raw_wide_final = met_raw_wide_kegg_annotated_tb_1
    # remove the met_raw_wide_kegg_annotated_tb_1 to save memory
    rm(met_raw_wide_kegg_annotated_tb_1)

    # Update the met_raw_wide_isotope as the met_raw_wide_isotope_kegg_annotated_tb_1 for later use
    met_raw_wide_final_isotope = met_raw_wide_isotope_kegg_annotated_tb_1
    # remove the met_raw_wide_isotope_kegg_annotated_tb_1 to save memory
    rm(met_raw_wide_isotope_kegg_annotated_tb_1)

    ######################### prepare the sample data in the long format
    # Write a function to combine two data frames (fill NA for the data frame with fewer rows) into a single data frame
    combine_data_frames <- function(df1, df2) {
        # Find the number of rows in both data frames
        nrows1 <- nrow(df1)
        nrows2 <- nrow(df2)

        # Identify which data frame has fewer rows and add NA rows
        ## If df1 has fewer rows than df2, add NA rows to df1
        if (nrows1 < nrows2) {
            na_rows <- data.frame(matrix(NA, ncol = ncol(df1), nrow = nrows2 - nrows1))
            colnames(na_rows) <- colnames(df1)
            df1 <- rbind(df1, na_rows)
        } 
        ## If df2 has fewer rows than df1, add NA rows to df2
        else if (nrows2 < nrows1) {
            na_rows <- data.frame(matrix(NA, ncol = ncol(df2), nrow = nrows1 - nrows2))
            colnames(na_rows) <- colnames(df2)
            df2 <- rbind(df2, na_rows)
        }

        # Combine the data frames into a single data frame
        result <- cbind(df1, df2)

        # convert to tb
        result = tibble(result)

        return(result)
    }

    # raw metabolomics data
    ## add column_mode_KEGGID_adduct_mz_time_sample column by combining KEGGID, mz (rounded to 4 decimals) and time (rounded to 1 decimal), separated by "_"
    met_raw_wide_final_1 = met_raw_wide_final %>%
        mutate(
            column_mode_KEGGID_adduct_mz_time_sample = paste0("annotated", "_", ion_mode, "_", " ", "_", KEGGID, "_", Adduct, "_", round(mz_sample, 4), "_", round(time_sample, 1))
        ) 

    # group by column_mode_KEGGID_adduct_mz_time_sample and select the feature with highest intensity_mean if there are duplicates. These are the duplicates after the mz and time rounding
    met_raw_wide_final_1 = met_raw_wide_final_1 %>%
        group_by(column_mode_KEGGID_adduct_mz_time_sample) %>%
        filter(intensity_mean == max(intensity_mean)) %>%
        ungroup()

    # only do the following if identified_target is NULL and metabolite_identified is not NULL: this is for the all metabolites analysis (using all identified metabolites as starting points)
    if (is.null(identified_target) & !is.null(metabolite_identified)) {
        # exclude those annotated metabolites that are identified already
        met_raw_wide_final_1 = met_raw_wide_final_1 %>%
            ## include only those rows in the KEGGID column that are not in the metabolite_identified
            filter(!(KEGGID %in% metabolite_identified$KEGGID))
    } 
    # if identified_target is not NULL and metabolite_identified is not NULL, only remove the identified_target using the KEGGID column: this is for using specific identified metabolites as the starting point
    else if (!is.null(identified_target) & !is.null(metabolite_identified)) {
        met_raw_wide_final_1 = met_raw_wide_final_1 %>%
            ## include only those rows in the KEGGID column that are not in the identified_target
            filter(!(KEGGID %in% identified_target))
    }

    # dplyr::select only the column_mode_KEGGID_adduct_mz_time_sample and columns starting at 14th column
    met_raw_wide_final_2 = met_raw_wide_final_1[, c(14:ncol(met_raw_wide_final_1))]

    # reorder the columns by moving column_mode_KEGGID_adduct_mz_time_sample to the first column
    met_raw_wide_final_2 = met_raw_wide_final_2 %>%
        dplyr::select(column_mode_KEGGID_adduct_mz_time_sample, everything())

    # transform the data from wide to long format by setting column_mode_KEGGID_adduct_mz_time_sample as the id column and all other columns as value columns
    met_raw_long = t(met_raw_wide_final_2[,2:ncol(met_raw_wide_final_2)])
    met_raw_long = as.data.frame(met_raw_long)
    colnames(met_raw_long) = met_raw_wide_final_2$column_mode_KEGGID_adduct_mz_time_sample
    met_raw_long_tb = tibble(met_raw_long)

    MSMICA_cor_input = met_raw_long_tb

    # if detail = TRUE, then save the MSMICA_cor_input as a csv file temporarily
    if (detail == TRUE) {
        MSMICA_cor_input_file_name <- paste0(prefix, "MSMICA_cor_input.csv")
        write_csv(MSMICA_cor_input, MSMICA_cor_input_file_name)
    }

    # if detail = TRUE, then create the following folders
    if (detail == TRUE){
        # Create a folder to store all the new_identified_metabolite files
        new_identified_folder_name <- "new_identified_metabolite"
        dir.create(new_identified_folder_name)

        # Create a folder to store all simplified connection files
        simplified_connection_annotation_folder_name <- "simplified_connection_annotation"
        dir.create(simplified_connection_annotation_folder_name)

        # Create a folder to store all MSMICA_decision_input CSV files
        MSMICA_decision_input_folder_name <- "MSMICA_decision_input_all"
        dir.create(MSMICA_decision_input_folder_name)
    }

    # if detail = TRUE, then create a sub-folder to store the MSMICA_decision_input CSV files for each run
    if (detail == TRUE) {
        # Create a sub-folder to store the MSMICA_decision_input CSV files for each run
        MSMICA_decision_input_run_folder_name <- paste0(MSMICA_decision_input_folder_name, "/", "MSMICA_decision_input_run", "_", 0)
        dir.create(MSMICA_decision_input_run_folder_name)
    }

    # select the necessary columns for bayesian probability calculation for metabolite adduct clusters
    met_raw_wide_final_3 = met_raw_wide_final_2 %>%
        # separate the column_mode_KEGGID_adduct_mz_time_sample column by "_" into 7 columns while keeping the original column: identification_type_annotated, ion_mode, Confirmed_Name, KEGGID, Adduct_annotated, mz_annotated, time_annotated
        separate_wider_delim(column_mode_KEGGID_adduct_mz_time_sample, delim = "_", names = c("identification_type_annotated", "ion_mode", "Confirmed_Name", "KEGGID", "Adduct_annotated", "mz_annotated", "time_annotated"), cols_remove = FALSE) %>%
        # rename column_mode_KEGGID_adduct_mz_time_sample as col_names_annotated
        rename(col_names_annotated = column_mode_KEGGID_adduct_mz_time_sample) %>%
        # make the splited mz_annotated and time_annotated columns as numeric
        mutate(
            mz_annotated = as.numeric(mz_annotated),
            time_annotated = as.numeric(time_annotated)
        ) %>%
        group_by(KEGGID) %>%
        filter(n() >= 2) %>%
        ungroup() %>%
        # only select the following columns: KEGGID, col_names_annotated, identification_type_annotated, ion_mode, Confirmed_Name, Adduct_annotated, mz_annotated, time_annotated
        dplyr::select(KEGGID, col_names_annotated, identification_type_annotated, ion_mode, Confirmed_Name, Adduct_annotated, mz_annotated, time_annotated)


    # Splitting the data by KEGGID
    data_split_cluster <- split(met_raw_wide_final_3, met_raw_wide_final_3$KEGGID)

    # Before the lapply function, create a new progress bar for bayesian probability calculation for clustering patterns:
    total_iterations_cluster <- length(data_split_cluster)
    pb_cluster <- txtProgressBar(min = 0, max = total_iterations_cluster, style = 3)

    i = 0
    KEGG_completion = c()
    last_progress_updete = 1

    # create the a function to calculate the bayesian probability for each metabolite adduct cluster
    bayesian_probability_calculation_cluster = function(data, MSMICA_cor_input) {

        # add current KEGGID to the KEGG_completion
        KEGG_completion <<- c(KEGG_completion, data$KEGGID[1])

        # print(data$KEGGID[1])

        # calculation the total iterations left
        i <<- i + 1
        total_iterations_left = total_iterations_cluster - i

        # update the progress bar when the i is divisible by 50
        if (i %% 50 == 0) {
            ## update the progress bar
            setTxtProgressBar(pb_cluster, i)
            # print the KEGG_completion from the last progress update to the current i
            print(KEGG_completion[last_progress_updete:i])
            ## save the current i
            last_progress_updete <<- i
        }

        ####################### adduct correlation analysis #######################
        # Compute the adduct correlation matrix only if there are more than 1 rows in the data data frame
        if(nrow(data) > 1) {
            MSMICA_cor_input_annotated <- MSMICA_cor_input[, data$col_names_annotated]
            full_cor_matrix_annotated <- suppressWarnings(cor(MSMICA_cor_input_annotated, MSMICA_cor_input_annotated, method="spearman", use="complete.obs"))

            # Function to find high correlation pairs
            find_high_correlation_pairs <- function(cor_matrix, threshold) {
                # Filter the correlation matrix by setting the lower triangle and diagonal to NA
                cor_matrix[lower.tri(cor_matrix)] <- NA
                diag(cor_matrix) <- NA

                high_cor_indices <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)

                # if there are no high correlation pairs, return a empty list
                if (length(high_cor_indices) == 0) {
                    return(list())
                }

                # Extract and sort the column names based on indices to get unique pairs
                ## this works because the indices are in the same order as the column names

                #                                                    annotated_positive_ _C00062_M-H2O+H_157.1088_107.3 annotated_positive_ _C00062_M-H2O+H_157.109_208.2 annotated_positive_ _C00062_M+H_175.1195_210
                # annotated_positive_ _C00062_M-H2O+H_157.1088_107.3                                                 NA                                        0.08051876                                    0.0230049
                # annotated_positive_ _C00062_M-H2O+H_157.109_208.2                                                  NA                                                NA                                    0.5183700
                # annotated_positive_ _C00062_M+H_175.1195_210                                                       NA                                                NA                                           NA
                # annotated_positive_ _C00062_M+Na_197.1016_204.1                                                    NA                                                NA                                           NA
                # annotated_positive_ _C00062_M+K_213.0753_105.4                                                     NA                                                NA                                           NA
                # annotated_positive_ _C00062_M+2Na-H_219.0836_203.4                                                 NA                                                NA                                           NA
                #                                                    annotated_positive_ _C00062_M+Na_197.1016_204.1 annotated_positive_ _C00062_M+K_213.0753_105.4 annotated_positive_ _C00062_M+2Na-H_219.0836_203.4
                # annotated_positive_ _C00062_M-H2O+H_157.1088_107.3                                    0.0008610155                                     -0.1324536                                        -0.09210378
                # annotated_positive_ _C00062_M-H2O+H_157.109_208.2                                     0.5145237072                                     -0.4565664                                         0.44300350
                # annotated_positive_ _C00062_M+H_175.1195_210                                          0.3992441537                                     -0.2486592                                         0.47886772
                # annotated_positive_ _C00062_M+Na_197.1016_204.1                                                 NA                                     -0.4200910                                         0.47138997
                # annotated_positive_ _C00062_M+K_213.0753_105.4                                                  NA                                             NA                                        -0.21157729
                # annotated_positive_ _C00062_M+2Na-H_219.0836_203.4                                              NA                                             NA                                                 NA

                #                                                     row col
                # annotated_positive_ _C00062_M-H2O+H_157.109_208.2   2   3

                ## for example, here the correlation between M-H2O+H_157.109_208.2 and C00062_M+H_175.1195_210 IS 0.51, so they are a high correlation pair. Thus, when we extract the column name 2 and 3, this is equal to extract the row name 2 and column name 3 in the original cor_matrix, as the row and column names are the same in the cor_matrix.

                high_cor_pairs <- apply(high_cor_indices, 1, function(idx) {
                    colnames(cor_matrix)[idx]
                })

                # Convert to unique pairs
                unique_pairs <- unique(apply(high_cor_pairs, 2, function(pair) {
                    sort(pair)
                }, simplify = FALSE))

                # Convert the list of unique pairs to a list of character vectors
                unique_pairs_list <- lapply(unique_pairs, function(pair) {
                    return(pair)
                })

                return(unique_pairs_list)
            }

            # Find all high correlation pairs
            intra_correlation_all <- find_high_correlation_pairs(full_cor_matrix_annotated, adduct_correlation_r_threshold)

            # Define a function to check rt_annotated max and min difference among all the intra correlation groups in the intra_correlation_all list
            check_time_range_group <- function(group) {
                # Filter the data frame for rows matching any value in the group
                filtered_df <- data %>% 
                    filter(col_names_annotated %in% group)
                
                # Calculate max and min time_annotated and their difference
                max_time <- max(filtered_df$time_annotated, na.rm = TRUE)
                min_time <- min(filtered_df$time_annotated, na.rm = TRUE)
                time_difference <- max_time - min_time
                
                # Return a data frame with the result for the group
                return(data.frame(group = paste(group, collapse = ", "), max_time = max_time, min_time = min_time, time_difference = time_difference))
            }


            # Iterate over the groups in intra_correlation_all list and apply the function. Then bind the results into a single data frame
            result_group <- do.call(rbind, lapply(intra_correlation_all, check_time_range_group))

            # do the following if the result_group data frame is not NULL: this is to filter the result_group data frame to retain only rows where time_difference is less than adduct_correlation_time_threshold (default is 5s) and then convert the group strings back to character vectors
            if (!is.null(result_group)){
                # Filter the result_group data frame to retain only rows where time_difference is less than adduct_correlation_time_threshold (default is 5s)
                filtered_groups <- result_group %>% 
                    filter(time_difference < adduct_correlation_time_threshold)

                # Convert the group strings back to character vectors
                filtered_group_vectors <- strsplit(as.character(filtered_groups$group), split = ", ")

                # Recursive function to compare each filtered_group_vector against elements in intra_correlation_all
                get_matching_groups <- function(intra_correlation_list, filtered_group_vector) {

                    matching_groups <- list()

                    # Loop through each element in intra_correlation_list
                    for (i in seq_along(intra_correlation_list)) {
                        element <- intra_correlation_list[[i]]

                        if (is.list(element)) {
                            matching_groups <- c(matching_groups, get_matching_groups(element, filtered_group_vector))
                        } 
                        else if (all(filtered_group_vector %in% element) && all(element %in% filtered_group_vector)) {
                            matching_groups <- c(matching_groups, list(element))
                        }
                    }
                    return(matching_groups)
                }

                # Apply the function to each vector in filtered_group_vectors
                filtered_original_groups <- lapply(filtered_group_vectors, function(group_vec) {
                    get_matching_groups(intra_correlation_all, group_vec)
                })

            } else {
                filtered_original_groups = list()
            }

            intra_correlation_all_clean = filtered_original_groups

            # Flatten the nested list
            flat_list <- do.call(c, intra_correlation_all_clean)

            intra_correlation_all_clean_df = NULL 

            # check if any intra connection > adduct_correlation_r_threshold found. If yes, then perform the following. If no, then return NULL
            if (length(flat_list) > 0) {
                # Transform to data frame
                intra_correlation_all_clean_df <- data.frame(cluster = rep(1:length(flat_list), sapply(flat_list, length)),
                    variable = unlist(flat_list))

                # rename the cluster column as adduct_corr_cluster
                intra_correlation_all_clean_df = intra_correlation_all_clean_df %>%
                    rename(adduct_corr_cluster = cluster)

                # Merge the data with the intra_correlation_all_clean_df to add the adduct_corr_cluster column and assign values for those within the adduct correlation cluster
                data_adduct_corr = data %>%
                    full_join(intra_correlation_all_clean_df, by = c("col_names_annotated" = "variable"), relationship = "many-to-many")
            } else {
                # if there is no adduct correlation cluster, then add NA to the adduct_corr_cluster column
                data_adduct_corr = data %>%
                    mutate(adduct_corr_cluster = NA)
            }

        } else {
            # if there is only 1 row in the annotation data, then add NA to the adduct_corr_cluster column (because 1 annotation cannot form a cluster, which has at least 2 adducts)
            data_adduct_corr = data %>%
                mutate(adduct_corr_cluster = NA)
        }

        # select the current annotated metabolite's KEGGID
        current_kegg_ids = unique(data_adduct_corr$KEGGID)

        # remove duplicate rows in the data_adduct_corr data frame based on KEGGID, Adduct, mz, and time
        data_adduct_corr_2 = data_adduct_corr %>%
            distinct(KEGGID, Adduct_annotated, mz_annotated, time_annotated, .keep_all = TRUE) %>%
            # replace the numeric values in the adduct_corr_cluster with TRUE
            mutate(adduct_corr_cluster = ifelse(!is.na(adduct_corr_cluster), TRUE, NA))

        ####################### isotopic adduct correlation analysis #######################

        # Initialize an empty data frame to store the ratios
        abundance_ratios <- data.frame()

        # Filter both datasets for the current KEGG ID
        primary_df = met_raw_wide_final %>%
            filter(KEGGID %in% current_kegg_ids)

        isotope_df = met_raw_wide_final_isotope %>%
            filter(KEGGID %in% current_kegg_ids)

        # if there is no isotopic adduct, then add NA to the mz_isotope and time_isotope columns
        if (nrow(isotope_df) == 0) {
            data_adduct_corr_3 = data_adduct_corr_2 %>%
                mutate(
                    mz_isotope = NA,
                    time_isotope = NA
                    )
        } 
        # if the isotopic adduct annotation at least has 1 same Adduct in the primary adduct, then perform the following
        else if (any(primary_df$Adduct %in% isotope_df$Adduct) == TRUE) {
            ## loop over each row in primary_df to calculate isotopic abundance ratio and correlation if there is any isotopic adduct
            for (i in 1:nrow(primary_df)) {

                current_primary_df = primary_df[i,]
                adduct_primary = current_primary_df$Adduct
                time_primary = current_primary_df$time_sample

                ## loop over each row in isotope_df
                for (j in 1:nrow(isotope_df)) {

                    current_isotope_df = isotope_df[j,]
                    adduct_isotope = current_isotope_df$Adduct
                    time_isotope = current_isotope_df$time_sample_isotope

                    time_difference = abs(time_isotope - time_primary)

                    ## if the adduct and isotope are the same and time_difference <= isotopic_correlation_time_threshold (default 5s), then calculate the ratio
                    if (adduct_primary == adduct_isotope & time_difference <= isotopic_correlation_time_threshold) {

                        # Assume that column names for sample intensities match between the two dataframes
                        sample_columns <- colnames(primary_df)[-(1:13)] # Assuming the first 13 columns are not sample intensities

                        primary_adduct_intensity = t(current_primary_df[1, sample_columns])

                        isotope_adduct_intensity = t(current_isotope_df[1, sample_columns])

                        ## calculate the ratio of the isotopic adduct to the primary adduct
                        ratios <- isotope_adduct_intensity / primary_adduct_intensity

                        ## create a dataframe to store the results
                        result_df <- data.frame(KEGGID = current_kegg_ids, mz_primary = current_primary_df$mz_sample, time_primary = current_primary_df$time_sample, mz_isotope = current_isotope_df$mz_sample_isotope, time_isotope = current_isotope_df$time_sample_isotope, Adduct = adduct_primary, primary = primary_adduct_intensity, isotope = isotope_adduct_intensity, predicted_ratio = current_isotope_df$predicted_isotopic_abundance_ratio, ratios)
                        
                        ## add the result to the abundance_ratios dataframe
                        abundance_ratios <- rbind(abundance_ratios, result_df)
                    }
                }
            }

            # if there is no isotopic adduct, then add NA to the mz_isotope and time_isotope columns
            if (nrow(abundance_ratios) == 0) {
                data_adduct_corr_3 = data_adduct_corr_2 %>%
                    mutate(
                        mz_isotope = NA,
                        time_isotope = NA
                        )
            } 
            else {
                # Convert it back to a tibble for nicer printing, and set row names as a new column
                abundance_ratios_2 <- as_tibble(abundance_ratios)

                # Convert the ratios to numeric, and replace any non-finite values with 0
                abundance_ratios_2$ratios[!is.finite(abundance_ratios_2$ratios)] <- NA

                # multiply the ratios by 100
                abundance_ratios_2$ratios = abundance_ratios_2$ratios * 100

                ## calculate the pearson correlation between the isotopic adduct and the primary adduct
                abundance_ratios_2 = abundance_ratios_2 %>%
                    mutate(
                        correlation = cor(primary, isotope, method = "pearson", use = "complete.obs")
                        )

                # Convert all 0 to NA
                abundance_ratios_2$ratios[abundance_ratios_2$ratios == 0] <- NA

                # summarize the abundance ratio using mean by KEGGID
                abundance_ratios_3 = abundance_ratios_2 %>%
                    group_by(KEGGID, mz_primary, time_primary, mz_isotope, time_isotope, Adduct) %>%
                    summarize(
                        mean_abundance_ratio = mean(ratios, na.rm = TRUE),
                        mean_predicted_abundance_ratio = mean(predicted_ratio, na.rm = TRUE),
                        mean_correlation = mean(correlation, na.rm = TRUE), 
                        .groups = "keep"
                    ) %>%
                    ungroup()

                # calculate mean_absolute_isotope_ratio_deviation (absolute difference between the mean of observed and predicted ratio
                abundance_ratios_3 = abundance_ratios_3 %>%
                    mutate(
                        mean_absolute_isotope_ratio_deviation = abs(mean_abundance_ratio - mean_predicted_abundance_ratio)
                    ) %>%
                    ungroup() %>%
                    ## remove the mean_absolute_isotope_ratio_deviation > mean_absolute_isotope_ratio_deviation (default is 5)
                    filter(mean_absolute_isotope_ratio_deviation <= mean_absolute_isotope_ratio_deviation_threshold)

                # remove the abundance_ratios_3 with mean_abundance_ratio > 100 if any (since the isotopic adduct should have lower intensity than the primary adduct overall)
                abundance_ratios_4 = abundance_ratios_3 %>%
                    filter(mean_abundance_ratio <= 100) %>%
                    ## remove mean_correlation < isotopic_correlation_r_threshold (default is 0.5)
                    filter(mean_correlation >= isotopic_correlation_r_threshold)

                if (nrow(abundance_ratios_4) > 0) {
                    # select only these columns: mz_primary, time_primary, KEGGID, mz_isotope, time_isotope, Adduct, mean_correlation
                    abundance_ratios_4 = abundance_ratios_4[, c("mz_primary", "time_primary", "KEGGID", "mz_isotope", "time_isotope", "Adduct", "mean_correlation")]

                    # rename the mean_correlation as correlation
                    abundance_ratios_4 = abundance_ratios_4 %>%
                        rename(correlation = mean_correlation)

                    # round the correlation to 2 decimal places
                    abundance_ratios_4$correlation = round(abundance_ratios_4$correlation, 2)

                    # round mz_primary to 4 decimal places
                    abundance_ratios_4$mz_primary = round(abundance_ratios_4$mz_primary, 4)

                    # round time_primary to 1 decimal place
                    abundance_ratios_4$time_primary = round(abundance_ratios_4$time_primary, 1)

                    # inner join the abundance_ratios_4 with the data_adduct_corr_2 to add the abundance_ratios_4 columns to the data_adduct_corr_2 by KEGGID and Adduct_annotated
                    data_adduct_corr_3 = data_adduct_corr_2 %>%
                        left_join(abundance_ratios_4, by = c("mz_annotated"="mz_primary", "time_annotated"="time_primary", "KEGGID", "Adduct_annotated"="Adduct"), relationship = "many-to-many")

                    # filter out those rows where the mz_annotated is NA
                    data_adduct_corr_3 = data_adduct_corr_3 %>%
                        filter(!is.na(mz_annotated))

                } else {
                    # if there is no isotopic adduct correlation, then add NA to the mz_isotope and time_isotope columns
                    data_adduct_corr_3 = data_adduct_corr_2 %>%
                        mutate(
                            mz_isotope = NA,
                            time_isotope = NA
                            )
                }
            }
        }
        # else, if there is no isotopic adduct that has the same adduct in the primary adduct, then add NA to the mz_isotope and time_isotope columns
        else {
            data_adduct_corr_3 = data_adduct_corr_2 %>%
                mutate(
                    mz_isotope = NA,
                    time_isotope = NA
                    )
        }

        # remove duplicates based on KEGGID, Adduct_annotated, mz_annotated, and time_annotated
        data_adduct_corr_3 = data_adduct_corr_3 %>%
            distinct(KEGGID, Adduct_annotated, mz_annotated, time_annotated, .keep_all = TRUE)

        # if there is no adduct correlation cluster and isotope correlation cluster (all adduct_corr_cluster column's value is NA or all mz_isotope's value is NA), then just return an empty data frame without even calculating the Bayesian probability because all we care about this function is to calculate the Bayesian probability for those metabolite adducts with adduct correlation cluster and isotope correlation cluster
        if (all(is.na(data_adduct_corr_3$adduct_corr_cluster)) & all(is.na(data_adduct_corr_3$mz_isotope))) {
            return(data.frame())
        }

        if (detail == TRUE){
            # create a folder for the current KEGGID if at least one of the adduct_corr_cluster or mz_isotope columns is not NA
            dir.create(paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids), showWarnings = FALSE)
        }
        
        ### check any adduct correlation results available. If intra_correlation_all_clean_df is not NULL (any adduct connection clusters), then perform the following. If NULL (no adduct connection clusters), then follow the next step
        if (!all(is.na(data_adduct_corr_3$adduct_corr_cluster))){
            # if detail is TRUE, then save the data frame with correlation coefficient
            if (detail == TRUE){
                ## Simplify the data using intra_correlation_all_clean_df to only include the adducts that are in the intra_correlation_all_clean_df
                data_final = data %>%
                    filter(col_names_annotated %in% intra_correlation_all_clean_df$variable)
                MSMICA_cor_input_annotated_final <- MSMICA_cor_input[, data_final$col_names_annotated]
                full_cor_matrix_annotated_final <- suppressWarnings(cor(MSMICA_cor_input_annotated_final, MSMICA_cor_input_annotated_final, method="spearman", use="complete.obs"))
                full_cor_matrix_annotated_final = as.data.frame(full_cor_matrix_annotated_final)
                ## save the correlation matrix as a csv file
                write.csv(full_cor_matrix_annotated_final, paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids, "/", current_kegg_ids, "_adduct_cor_matrix.csv"))
            }
        }

        # if detail is TRUE and there is isotope correlation, then save the data frame with correlation coefficient
        if (detail == TRUE & !all(is.na(data_adduct_corr_3$mz_isotope))){
            abundance_ratios_4 = as.data.frame(abundance_ratios_4)
            ## save the abundance_ratios_4 as a csv file
            write_csv(abundance_ratios_4, paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids, "/", current_kegg_ids, "_isotope_correlation.csv"))
        }

        # calculate the Bayesian posterior probability using the following evidence:
        ## evidence 1: adduct formation
        ## evidence 2: adduct correlation cluster
        ## evidence 3: isotopic adduct correlation
        bayesian_probability_calculation_cluster_child <- function(data) {

            # initiate the data frame to store the probability combination results with 0 and 1 as the values
            ## Generate all combinations of 0 and 1 for each column
            combinations_probability <- expand.grid(rep(list(0:1), nrow(data)))
            colnames(combinations_probability) <- data$col_names_annotated                

            # remove the rows with all 0s
            combinations_probability <- combinations_probability[rowSums(combinations_probability) > 0, ]

            ################################# adduct formation probability calculation #################################

            # Get the adduct for the current result
            adduct <- data$Adduct_annotated

            # if the ion_mode is the positive mode, then perform the following:
            if (unique(data$ion_mode) == "positive") {
                p_M_H <- 440/524
                p_M_Other <- 84/524
            }
            # else, if the ion_mode is the negative mode, then perform the following:
            else if (unique(data$ion_mode) == "negative") {
                p_M_H <- 418/477
                p_M_Other <- 59/477
            }

            # Function to calculate the probability for each combination
            adduct_form_calculate_probability <- function(column) {
                p_adduct <- sapply(seq_along(column), function(i) {
                    if (adduct[i] == "M+H" | adduct[i] == "M-H") {
                        ifelse(column[i] == 1, p_M_H, 1 - p_M_H)
                    } else {
                        ifelse(column[i] == 1, p_M_Other, 1 - p_M_Other)
                    }
                })
                return(prod(p_adduct))
            }


            ################################# adduct correlation probability calculation #################################

            # Get the adduct correlation cluster for the current result
            adduct_correlation_cluster <- data$adduct_corr_cluster

            # Define the probabilities for each adduct correlation cluster being present (1)
            p_within_adduct_cor <- 63/81
            p_outside_adduct_cor <- 40/160

            # Function to calculate the probability for each combination
            adduct_correlation_calculate_probability <- function(column) {
                p_adduct_correlation <- sapply(seq_along(column), function(i) {
                    if (!is.na(adduct_correlation_cluster[i])) {
                        ifelse(column[i] == 1, p_within_adduct_cor, 1 - p_within_adduct_cor)
                    } else {
                        ifelse(column[i] == 1, p_outside_adduct_cor, 1 - p_outside_adduct_cor)
                    }
                })
                return(prod(p_adduct_correlation))
            }

            ################################# isotope correlation probability calculation #################################

            # Get the isotope correlation for the current result
            isotope_correlation <- data$mz_isotope

            # Define the probabilities for each isotope correlation being present (1)
            p_isotope_yes <- 30/30
            p_isotope_no <- 14/114

            # Function to calculate the probability for each combination
            isotope_correlation_calculate_probability <- function(column) {
                p_isotope_correlation <- sapply(seq_along(column), function(i) {
                    if (!is.na(isotope_correlation[i])) {
                        ifelse(column[i] == 1, p_isotope_yes, 1 - p_isotope_yes)
                    } else {
                        ifelse(column[i] == 1, p_isotope_no, 1 - p_isotope_no)
                    }
                })
                return(prod(p_isotope_correlation))
            }

            # create a final combination probability data frame to store the results
            combinations_probability_final = combinations_probability

            # Apply the adduct_form_calculate_probability function to each combination to calculate the adduct_formation_probability of adduct formation for each metabolite adduct
            combinations_probability_final$adduct_formation_probability <- apply(combinations_probability, 1, adduct_form_calculate_probability)
            # Apply the adduct_correlation_calculate_probability function to each combination to calculate the adduct_correlation_probability
            combinations_probability_final$adduct_correlation_probability <- apply(combinations_probability, 1, adduct_correlation_calculate_probability)
            # Apply the isotope_correlation_calculate_probability function to each combination to calculate the probability
            combinations_probability_final$isotope_correlation_probability <- apply(combinations_probability, 1, isotope_correlation_calculate_probability)

            # normalize the adduct_formation_probability to sum to 1
            combinations_probability_final$adduct_formation_probability <- combinations_probability_final$adduct_formation_probability / sum(combinations_probability_final$adduct_formation_probability)
            # normalize the adduct_correlation_probability to sum to 1
            combinations_probability_final$adduct_correlation_probability <- combinations_probability_final$adduct_correlation_probability / sum(combinations_probability_final$adduct_correlation_probability)
            # normalize the isotope_correlation_probability to sum to 1
            combinations_probability_final$isotope_correlation_probability <- combinations_probability_final$isotope_correlation_probability / sum(combinations_probability_final$isotope_correlation_probability)

            # calculate the final probability for each combination by multiplying them together
            combinations_probability_final$mean_raw_posterior <- combinations_probability_final$adduct_formation_probability * combinations_probability_final$adduct_correlation_probability * combinations_probability_final$isotope_correlation_probability

            # normalize the mean_raw_posterior to sum to 1
            combinations_probability_final$mean_normalized_posterior <- combinations_probability_final$mean_raw_posterior / sum(combinations_probability_final$mean_raw_posterior)

            # arrange the combinations_probability_final by the mean_normalized_posterior in descending order
            combinations_probability_final <- combinations_probability_final[order(-combinations_probability_final$mean_normalized_posterior),]

            # if detail is TRUE, then save the combinations_probability_final to a csv file
            if (detail) {
                combinations_probability_final = as.data.frame(combinations_probability_final)
                write_csv(combinations_probability_final, paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids, "/", current_kegg_ids, "_bayesian_probability.csv"))
            }

            max_mean_normalized_posterior = max(combinations_probability_final$mean_normalized_posterior)

            # select the rows with the highest mean_normalized_posterior
            highest_posterior_probability_adduct <- combinations_probability_final %>%
                dplyr::filter(mean_normalized_posterior == max_mean_normalized_posterior) %>%
                ## remove the adduct_formation_probability, adduct_correlation_probability, isotope_correlation_probability, enzyme_correlation_probability, mean_raw_posterior, and mean_normalized_posterior columns
                dplyr::select(-adduct_formation_probability, -adduct_correlation_probability, -isotope_correlation_probability, -mean_raw_posterior, -mean_normalized_posterior)

            # when there are multiple rows with the same highest mean_normalized_posterior, then merge the rows into a single row following this logic:
            ## take the average of all columns and round up to the nearest integer. This essentially resolves the conflicts when one row has 1 for one metabolite adduct and the other has 0 for the same metabolite adduct
            highest_posterior_probability_adduct_final <- highest_posterior_probability_adduct %>%
                summarise(across(everything(), ~ceiling(mean(.)))) %>%
                ## transform the data frame into a long format with the col_names_annotated as the column name and the value as the value
                pivot_longer(cols = everything(), names_to = "col_names_annotated", values_to = "MSMICA_identification")

            # left join the highest_posterior_probability_adduct with the original data to get the final result using the col_names_annotated as the key
            final_result <- left_join(data, highest_posterior_probability_adduct_final, by = "col_names_annotated", relationship = "many-to-many")

            return(final_result)
        }

        # apply the bayesian_probability_calculation_cluster_child function to the data_adduct_corr_3 data frame
        final_result_cluster <- bayesian_probability_calculation_cluster_child(data_adduct_corr_3)

        return(final_result_cluster)
    }

    print("Calculating Bayesian probability for each metabolite using adduct correlation, isotope correlation, and adduct formation:")

    # Apply the bayesian_probability_calculation_cluster function to each group
    results_cluster <- lapply(data_split_cluster, function(x) bayesian_probability_calculation_cluster(x, MSMICA_cor_input))

    # Combine all data frames in the list into one final data frame
    final_results_cluster <- bind_rows(results_cluster)
    final_results_cluster = tibble(final_results_cluster)

    # select only the rows with MSMICA_identification == 1
    final_results_cluster_identified = final_results_cluster %>%
        filter(MSMICA_identification == 1)

    # select only the KEGGID and Name columns
    KEGG_database_simple = KEGG_database %>%
        dplyr::select(KEGGID, Name) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1))

    # left join the final_results_cluster_identified with the KEGG_database_simple to get the final result using KEGGID
    final_results_cluster_identified_2 = left_join(final_results_cluster_identified, KEGG_database_simple, by = "KEGGID", relationship = "many-to-many")

    # select the following columns to form the identified metabolite list: ion_mode, Name, KEGGID, Adduct, mz_annotated, time_annotated
    final_results_cluster_identified_3 = final_results_cluster_identified_2 %>%
        dplyr::select(ion_mode, Name, KEGGID, Adduct_annotated, mz_annotated, time_annotated) %>%
        ## rename the columns: Name as metabolite, Adduct as adduct, mz_annotated as mz, time_annotated as time
        rename(metabolite = Name, adduct = Adduct_annotated, mz = mz_annotated, time = time_annotated)

    # round the mz to 4 decimal places
    final_results_cluster_identified_3$mz = round(final_results_cluster_identified_3$mz, 4)
    # round the time to 1 decimal place
    final_results_cluster_identified_3$time = round(final_results_cluster_identified_3$time, 1)

    # round the feature table mz to 4 decimal places
    met_raw_wide$mz = round(met_raw_wide$mz, 4)
    # round the feature table time to 1 decimal place
    met_raw_wide$time = round(met_raw_wide$time, 1)

    # left join the final_results_test_identified_3 with the met_raw_wide to get the intensity values for the identified metabolites using mz and time as the keys
    final_results_cluster_identified_4 = left_join(final_results_cluster_identified_3, met_raw_wide, by = c("mz", "time"), relationship = "many-to-many")
    final_results_cluster_identified_4 = tibble(final_results_cluster_identified_4)

    # create a column, MSMICA_identification, with value 1
    final_results_cluster_identified_4$MSMICA_identification = 1

    # move the MSMICA_identification column to the 7th column, right after the time column
    final_results_cluster_identified_4 = final_results_cluster_identified_4 %>%
        dplyr::relocate(MSMICA_identification, .after = time)

    print("Here are the identified metabolites by using the evidence of adduct correlation, isotope correlation, and adduct formation:")
    print(final_results_cluster_identified_4)

    # check if the metabolite_identified is provided by the user and is not NULL
    if (!is.null(metabolite_identified)) {
        # if the metabolite_identified is provided by the user, then use the metabolite_identified as the starting point
        print("The metabolite_identified is provided by the user")
        print(head(metabolite_identified))

        metabolite_identified_2 = metabolite_identified

        # create a column, MSMICA_identification, with value 0
        metabolite_identified_2$MSMICA_identification = 0
    
        # remove the metabolites present in the user-provided metabolite_identified from the final_results_cluster_identified_4
        final_results_cluster_identified_4 = final_results_cluster_identified_4 %>%
            filter(!(KEGGID %in% metabolite_identified$KEGGID))

        # move the MSMICA_identification column to the 7th column, right after the time column
        metabolite_identified_2 = metabolite_identified_2 %>%
            dplyr::relocate(MSMICA_identification, .after = time)

        # remove duplicates
        metabolite_identified_2 = distinct(metabolite_identified_2)

        ################################## select identified metabolites as the starting point to test the MSMICA algorithm ##################################
        # if the identified_target is not NULL, then select the identified_target as the starting point
        if (!is.null(identified_target)) {
            print(paste0("The identified metabolites provided by the user are selected as the starting point for the MSMICA algorithm"))
            print(identified_target)

            metabolite_identified_2 = metabolite_identified_2 %>%
                filter(KEGGID %in% identified_target)
        } 
        # if the identified_target is NULL, then elect all identified metabolites as the starting point
        else {
            print("All identified metabolites provided by the user are selected as the starting point for the MSMICA algorithm")
        }

        ## add column_mode_KEGGID_adduct_mz_time_sample column by combining ion_mode, KEGGID, Adduct, mz_sample, time_sample, MSMICA_identification, identification_method, separated by "_"
        metabolite_identified_2 = metabolite_identified_2 %>%
            mutate(
                column_mode_KEGGID_adduct_mz_time_sample = paste0("identified_", ion_mode, "_", metabolite, "_", KEGGID, "_", adduct, "_", round(mz, 4), "_", round(time, 1), "_", MSMICA_identification, "_", "user")
            )
    } 
    ## if the identified metabolite list is not provided by user, then initiate an empty data frame
    else {
        print("The user does not provide any identified metabolites as the starting point for the MSMICA algorithm")
        metabolite_identified_2 = data.frame()
    }

    ## add column_mode_KEGGID_adduct_mz_time_sample column by combining ion_mode, KEGGID, Adduct, mz_sample, time_sample, MSMICA_identification, identification_method, separated by "_"
    final_results_cluster_identified_4 = final_results_cluster_identified_4 %>%
        mutate(
            column_mode_KEGGID_adduct_mz_time_sample = paste0("identified_", ion_mode, "_", metabolite, "_", KEGGID, "_", adduct, "_", round(mz, 4), "_", round(time, 1), "_", MSMICA_identification, "_", "clustering of adducts and isotopes AND precursor-product correlation")
        )

    # combine the final_results_cluster_identified_4 with the metabolite_identified_2 by rows
    metabolite_identified_FINAL = bind_rows(metabolite_identified_2, final_results_cluster_identified_4)

    # convert to tibble
    metabolite_identified_FINAL_1 = tibble(metabolite_identified_FINAL)

    # dplyr::select only the column_mode_KEGGID_adduct_mz_time_sample and columns representing intensity (now the column_mode_KEGGID_adduct_mz_time_sample is in the last column so it is included)
    metabolite_identified_FINAL_2 = metabolite_identified_FINAL_1[, c(8:ncol(metabolite_identified_FINAL_1))]

    # reorder the columns by moving column_mode_KEGGID_adduct_mz_time_sample to the first column
    metabolite_identified_FINAL_2 = metabolite_identified_FINAL_2 %>%
        dplyr::select(column_mode_KEGGID_adduct_mz_time_sample, everything())

    # remove duplicates based on column_mode_KEGGID_adduct_mz_time_sample
    metabolite_identified_FINAL_2 = metabolite_identified_FINAL_2 %>%
        distinct(column_mode_KEGGID_adduct_mz_time_sample, .keep_all = TRUE)

    # transform the data from wide to long format by setting column_mode_KEGGID_adduct_mz_time_sample as the id column and all other columns as value columns
    metabolite_identified_FINAL_2_long = t(metabolite_identified_FINAL_2[,2:ncol(metabolite_identified_FINAL_2)])
    metabolite_identified_FINAL_2_long = as.data.frame(metabolite_identified_FINAL_2_long)
    colnames(metabolite_identified_FINAL_2_long) = metabolite_identified_FINAL_2$column_mode_KEGGID_adduct_mz_time_sample
    metabolite_identified_FINAL_2_long_tb = tibble(metabolite_identified_FINAL_2_long)
    # combine the data frames
    total_identified_long <- metabolite_identified_FINAL_2_long_tb

    # remove previously created data frames to save memory
    rm(metabolite_identified_FINAL_2_long_tb, metabolite_identified_FINAL_2_long, metabolite_identified_FINAL_1, metabolite_identified_FINAL_2)

    # exclude those annotated metabolites that are identified already in the metabolite_identified_FINAL
    met_raw_wide_final_1 = met_raw_wide_final_1 %>%
        ## include only those rows in the KEGGID column that are not in the metabolite_identified
        filter(!(KEGGID %in% metabolite_identified_FINAL$KEGGID))

    # dplyr::select only the column_mode_KEGGID_adduct_mz_time_sample and columns starting at 14th column
    met_raw_wide_final_2 = met_raw_wide_final_1[, c(14:ncol(met_raw_wide_final_1))]

    # reorder the columns by moving column_mode_KEGGID_adduct_mz_time_sample to the first column
    met_raw_wide_final_2 = met_raw_wide_final_2 %>%
        dplyr::select(column_mode_KEGGID_adduct_mz_time_sample, everything())
    
    # remove duplicates based on column_mode_KEGGID_adduct_mz_time_sample
    met_raw_wide_final_2 = met_raw_wide_final_2 %>%
        distinct(column_mode_KEGGID_adduct_mz_time_sample, .keep_all = TRUE)

    # transform the data from wide to long format by setting column_mode_KEGGID_adduct_mz_time_sample as the id column and all other columns as value columns
    met_raw_long = t(met_raw_wide_final_2[,2:ncol(met_raw_wide_final_2)])
    met_raw_long = as.data.frame(met_raw_long)
    colnames(met_raw_long) = met_raw_wide_final_2$column_mode_KEGGID_adduct_mz_time_sample
    met_raw_long_tb = tibble(met_raw_long)
    # remove previously created data frames to save memory
    rm(met_raw_long, met_raw_wide_final_1)

    # if detail = TRUE, then save the final_results_1 as a csv file temporarily
    if (detail == TRUE) {
        final_results_1_file_name <- paste0(new_identified_folder_name, "/", "new_identified_metabolite", "_", 0, ".csv")
        metabolite_identified_FINAL_simple = metabolite_identified_FINAL[, 1:7]
        write_csv(metabolite_identified_FINAL_simple, final_results_1_file_name)
    }

    ############################# finally, combine the previous MSMICA_cor_input with the with the identified metabolite data if needed
    MSMICA_cor_input = combine_data_frames(total_identified_long, met_raw_long_tb)
    # remove previously created data frames to save memory
    rm(total_identified_long, met_raw_long_tb)

    # create a data frame with the column names of MSMICA_cor_input
    MSMICA_cor_input_col_names = data.frame(col_names = colnames(MSMICA_cor_input))
    
    # split the column names into multiple columns: identification_type, ion_mode, KEGGID, Adduct, mz, time, MSMICA_identification, identification_method
    MSMICA_cor_input_col_names_2 = MSMICA_cor_input_col_names %>%
        separate(col = col_names, into = c("identification_type", "ion_mode", "identified_Name", "KEGGID", "Adduct", "mz", "time", "MSMICA_identification", "identification_method"), sep = "_", remove = FALSE, fill = "right")

    # remove previously created data frames to save memory
    rm(MSMICA_cor_input_col_names)

    # select only the columns with "identified" in the identification_type column
    identified_metabolite = MSMICA_cor_input_col_names_2[MSMICA_cor_input_col_names_2$identification_type == "identified", ] %>%
        arrange(KEGGID) %>%
        tibble()

    # select only the columns with "annotated" in the identification_type column
    annotated_metabolite = MSMICA_cor_input_col_names_2[MSMICA_cor_input_col_names_2$identification_type == "annotated", ] %>%
        arrange(KEGGID) %>%
        tibble()


    # remove previously created data frames to save memory
    rm(MSMICA_cor_input_col_names_2)

    ## create a data frame with the mean intensity values
    mean_MSMICA_1 = colMeans(MSMICA_cor_input, na.rm = TRUE)
    mean_MSMICA_2 = data.frame(
        col_names = names(mean_MSMICA_1),
        ## round the mean intensity values to integers
        mean_intensity = round(mean_MSMICA_1, 0)
        )
    # remove previously created data frames to save memory
    rm(mean_MSMICA_1)

    # Before the while loop, create a new progress bar:
    total_iterations <- nrow(annotated_metabolite)
    pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

    print("This is current identified metabolite as the starting points for the MSMICA algorithm:")
    print(identified_metabolite)

    # add correlation values as NA in the correlation column to the identified metabolites
    identified_metabolite = identified_metabolite %>%
        mutate(correlation = NA) %>%
        ## move correlation column right after time
        relocate(correlation, .after = "time")

    ########## connect previously identified or MSMICA identified metabolites to the currently annotated metabolites based on kegg enzyme reaction connection ###############
    ###### this is a loop until all the identified or MSMICA identified metabolites are connected to the annotated metabolites

    # set the iteration number to 0 to count the number of iterations
    a = 0

    # set the while loop to run until there are no more annotated metabolites or until the final_results_1 is empty (which means no more newly identified metabolites)
    while(nrow(annotated_metabolite) > 0) {
        # add 1 to the iteration number
        a = a + 1
        print(paste0("This is iteration number ", a))

        # In the loop
        setTxtProgressBar(pb, total_iterations - nrow(annotated_metabolite))
        print(paste0("There are ", nrow(annotated_metabolite), " annotated metabolite adducts left to be identified"))

        # Perform a left join of the kegg_connection database to the identified_metabolite database
        # Allowing cartesian join by creating all possible combinations of rows
        MSMICA_col_names_connection <- identified_metabolite %>%
            left_join(kegg_connection, by = c("KEGGID" = "connection_1"), relationship = "many-to-many")

        # Arrange the data by the connection_2 column in ascending order
        MSMICA_col_names_connection_simplifed_2 <- MSMICA_col_names_connection %>%
            arrange(connection_2)

        # reorganize the columns in MSMICA_col_names_connection_simplifed: col_names, identification_type, KEGGID, Adduct, mz, time, identification_method, connection_2
        MSMICA_col_names_connection_simplifed_2 = MSMICA_col_names_connection_simplifed_2[, c("col_names", "identification_type", "KEGGID", "Adduct", "mz", "time", "identification_method", "connection_2")]

        # Ensure that all columns used for joining are characters
        MSMICA_col_names_connection_simplifed_2$connection_2 <- as.character(MSMICA_col_names_connection_simplifed_2$connection_2)
        annotated_metabolite$KEGGID <- as.character(annotated_metabolite$KEGGID)

        # inner join the MSMICA_TCA_col_names_connection with the annotated data
        MSMICA_col_names_connection_annotated = inner_join(
            MSMICA_col_names_connection_simplifed_2, 
            annotated_metabolite, 
            by = c("connection_2" = "KEGGID"),
            suffix = c("_identified", "_annotated"),
            relationship = "many-to-many"
            )

        print(paste0("There are ", nrow(MSMICA_col_names_connection_annotated), " precursor-product correlations built"))

        # if there is no new metabolite identified via enzyme based reaction on Kegg, then break the loop
        if (nrow(MSMICA_col_names_connection_annotated) == 0) {
            print("There is no more new metabolite that can be identified using enzyme-based reaction correlation ranking. The MSMICA algorithm has finished.")
            break
        }


        ############ calculation correlation coefficients and p-values only for the necessary identified and annotated metabolites ############

        # write a loop to compare the identified metabolites with the annotated metabolites by using the correlation test and create a data frame with the correlation values and p-values
        # create a data frame to store the correlation values and p-values
        cor_MSMICA_data = data.frame(
            col_names_identified = character(),
            col_names_annotated = character(),
            correlation = numeric(),
            p_value = numeric()
        )

        ## loop through each row of the MSMICA_col_names_connection_annotated to perform correlation test between col_names_identified and col_names_annotated by using subset the corresponding columns from MSMICA_cor_input
        for (i in 1:nrow(MSMICA_col_names_connection_annotated)){
            current_row = MSMICA_col_names_connection_annotated[i, ]

            identified_metabolite_current = current_row$col_names_identified
            annotated_metabolite_current = current_row$col_names_annotated

            identified_matrix <- unlist(MSMICA_cor_input[, identified_metabolite_current])
            other_matrix <- unlist(MSMICA_cor_input[, annotated_metabolite_current])

            cor_MSMICA = cor(identified_matrix, other_matrix, method="spearman", use="complete.obs")
            cor_p_MSMICA <- suppressWarnings(cor.test(identified_matrix, other_matrix, method = "spearman", use = "complete.obs"))

            # create a data frame with col_names_identified, col_names_annotated, correlation, p_value
            cor_MSMICA_data_current = data.frame(
                col_names_identified = identified_metabolite_current,
                col_names_annotated = annotated_metabolite_current,
                correlation = cor_MSMICA,
                p_value = cor_p_MSMICA$p.value
            )

            # combine the cor_MSMICA_data with cor_MSMICA_data_1
            cor_MSMICA_data = rbind(cor_MSMICA_data, cor_MSMICA_data_current)
        }

        # filter out those with p_value > 0.05
        cor_MSMICA_1 = cor_MSMICA_data %>%
            filter(p_value < 0.05)

        # left join the cor_MSMICA_1 to the MSMICA_col_names_connection_annotated to add correlation values
        ## now, we are basically adding the correlation values between the identified metabolites and the annotated metabolites from the previously created correlation matrix to the MSMICA_col_names_connection_annotated
        ## here, we just know that the current row of the MSMICA_col_names_connection_annotated is connected to all the identified metabolites, but we don't know which identified metabolite is connected to which annotated metabolite based on enzyme-based reaction on KEGG
        ## this is the first step to allow us to see the correlation values between the identified metabolites and the annotated metabolites for enzyme-based reaction correlation ranking
        MSMICA_col_names_connection_annotated_2 = left_join(MSMICA_col_names_connection_annotated, cor_MSMICA_1, by = c("col_names_annotated", "col_names_identified"), relationship = "many-to-many") %>%
            # left join the mean_MSMICA_2 to the MSMICA_TCA_col_names_connection_annotated_1 to add the mean intensity values of the annotated metabolites
            left_join(mean_MSMICA_2, by = c("col_names_annotated" = "col_names"), relationship = "many-to-many")

        # if there is no new metabolite identified via enzyme based reaction on Kegg, then break the loop
        if (nrow(MSMICA_col_names_connection_annotated_2) == 0) {
            print("There is no more new metabolite that can be identified using enzyme-based reaction correlation ranking. The MSMICA algorithm has finished.")
            break
        }

        # remove rows with NA in the correlation column
        MSMICA_col_names_connection_annotated_2 = MSMICA_col_names_connection_annotated_2 %>%
            filter(!is.na(correlation))

        # if detail = TRUE, then save the MSMICA_col_names_connection_annotated_2 as a csv file temporarily
        if (detail == TRUE) {
            MSMICA_col_names_connection_annotated_2_file_name <- paste0(simplified_connection_annotation_folder_name, "/", "MSMICA_enzyme_reaction_annotated", "_", a, ".csv")
            write_csv(MSMICA_col_names_connection_annotated_2, MSMICA_col_names_connection_annotated_2_file_name)
        }

        # arrange the correlation from high to low by each KEGGID and connection_2
        MSMICA_col_names_connection_annotated_2 = MSMICA_col_names_connection_annotated_2 %>%
            arrange(KEGGID, connection_2, desc(correlation))

        # arrange the data by connection_2 and correlation from high to low
        MSMICA_col_names_connection_annotated_2 = MSMICA_col_names_connection_annotated_2 %>%
            arrange(connection_2, desc(correlation))
        
        # round the correlation to 4 decimal places
        MSMICA_col_names_connection_annotated_2$correlation = round(MSMICA_col_names_connection_annotated_2$correlation, 4)

        # Filter out those self-correlation values (features with correlation of 1 with themselves)
        ## if correlation == 1, make it NA
        MSMICA_col_names_connection_annotated_2 = MSMICA_col_names_connection_annotated_2 %>%
            mutate(correlation = ifelse(correlation == 1, NA, correlation))

        # Summarize the correlation to mean_correlation by connection_2, col_names_annotated, identification_type_annotated, ion_mode, identified_Name, Adduct_annotated, mz_annotated, time_annotated
        MSMICA_col_names_connection_annotated_2 = MSMICA_col_names_connection_annotated_2 %>%
            group_by(connection_2, col_names_annotated, identification_type_annotated, ion_mode, identified_Name, Adduct_annotated, mz_annotated, time_annotated) %>%
            summarize(correlation = mean(correlation, na.rm = TRUE), .groups = "keep") %>%
            ungroup()
        
        # filter out those correlation values that are NA
        MSMICA_col_names_connection_annotated_2 = MSMICA_col_names_connection_annotated_2 %>%
            filter(!is.na(correlation))

        # Splitting the data by connection_2
        data_split <- split(MSMICA_col_names_connection_annotated_2, MSMICA_col_names_connection_annotated_2$connection_2)


        # if detail = TRUE, then create a sub-folder to store the MSMICA_decision_input CSV files for each run
        if (detail == TRUE) {
            # Create a sub-folder to store the MSMICA_decision_input CSV files for each run
            MSMICA_decision_input_run_folder_name <- paste0(MSMICA_decision_input_folder_name, "/", "MSMICA_decision_input_run", "_", a)
            dir.create(MSMICA_decision_input_run_folder_name)
        }

        # Function of the core of MSMICA algorithm
        bayesian_probability_calculation_enzyme <- function(data) {

            # select the current annotated metabolite's KEGGID
            current_kegg_ids = unique(data$connection_2)

            # Make the correlation as the absolute value of the correlation - this is to ensure that the correlation ranking is correct for both positive and negative correlation
            data$correlation = abs(data$correlation)

            # Convert mz_annotated from character to numeric
            data$mz_annotated = as.numeric(data$mz_annotated)

            # Convert time_annotated from character to numeric
            data$time_annotated <- as.numeric(data$time_annotated)

            # return NULL if no rows in the data frame
            if (nrow(data) == 0) {
                return(NULL)
            }

            # if detail = TRUE, then create a folder for the current KEGGID
            if (detail) {
                dir.create(paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids), showWarnings = FALSE)
            }

            # if there is only 1 row in the data frame, then return data because there is no need to calculate the likelihood of valid metabolite identification using enzyme correlation ranking if it is the only annotation using m/z 5 ppm tolerance and passed the correlation test
            if (nrow(data) == 1) {
                # add precursor-product correlation to identification_method and set MSMICA_identification to 1
                data = data %>%
                    mutate(
                        MSMICA_identification = 1,
                        identification_method = "precursor-product correlation"
                    )
                
                # if detail is TRUE, then save the data with enzyme_correlation_rank to a csv file
                if (detail) {
                    write_csv(data, paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids, "/", current_kegg_ids, "_enzyme_based_correlation_ranking.csv"))
                }

                return(data)
            }


            ####################### enzyme-based reaction correlation analysis #######################
            # the enzyme-based reaction correlation ranking is from high to low. so we need to first create a column to rank the correlation in descending order. Also, we need to rank the correlation by Adduct_annotated, because the correlation ranking is designed to compare different adducts for the same metabolite, not the same adduct
            data = data %>%
                arrange(desc(correlation)) %>%
                mutate(
                    enzyme_correlation_rank = row_number(desc(correlation))
                )

            # if detail is TRUE, then save the data with enzyme_correlation_rank to a csv file
            if (detail) {
                write_csv(data, paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids, "/", current_kegg_ids, "_enzyme_based_correlation_ranking.csv"))
            }

            # Get the enzyme reaction ranking for the current result
            enzyme_reaction_ranking <- data$enzyme_correlation_rank

            # initiate the data frame to store the probability combination results with 0 and 1 as the values
            ## Generate all combinations of 0 and 1 for each column
            combinations_probability <- expand.grid(rep(list(0:1), nrow(data)))
            colnames(combinations_probability) <- data$col_names_annotated                

            # remove the rows with all 0s
            combinations_probability <- combinations_probability[rowSums(combinations_probability) > 0,]

            # Define the probabilities for each enzyme reaction ranking being present (1)
            p_enzyme_first <- 71671/126239
            p_enzyme_second <- 35589/126239
            p_enzyme_other <- 18979/158635

            # Function to calculate the probability for each combination
            bayesian_probability_calculation_enzyme_child <- function(column) {
                p_enzyme_reaction_ranking <- sapply(seq_along(column), function(i) {
                ## if the enzyme_reaction_ranking is 1, then the probability is p_enzyme_first, otherwise, it is 1 - p_enzyme_first
                if (enzyme_reaction_ranking[i] == 1) {
                    ifelse(column[i] == 1, p_enzyme_first, 1 - p_enzyme_first)
                } 
                ## if the enzyme_reaction_ranking is 2, then the probability is p_enzyme_second, otherwise, it is 1 - p_enzyme_second
                else if (enzyme_reaction_ranking[i] == 2) {
                    ifelse(column[i] == 1, p_enzyme_second, 1 - p_enzyme_second)
                }
                else {
                    ifelse(column[i] == 1, p_enzyme_other, 1 - p_enzyme_other)
                }
                })
                return(prod(p_enzyme_reaction_ranking))
            }

            # create a final combination probability data frame to store the results
            combinations_probability_final = combinations_probability

            # Apply the bayesian_probability_calculation_enzyme_child function to each combination to calculate the enzyme_correlation_probability
            combinations_probability_final$enzyme_correlation_probability <- apply(combinations_probability, 1, bayesian_probability_calculation_enzyme_child)

            # normalize the enzyme_correlation_probability to sum to 1
            combinations_probability_final$enzyme_correlation_probability <- combinations_probability_final$enzyme_correlation_probability / sum(combinations_probability_final$enzyme_correlation_probability)

            # calculate the final probability for each combination by multiplying them together
            combinations_probability_final$mean_raw_posterior <- combinations_probability_final$enzyme_correlation_probability

            # normalize the mean_raw_posterior to sum to 1
            combinations_probability_final$mean_normalized_posterior <- combinations_probability_final$mean_raw_posterior / sum(combinations_probability_final$mean_raw_posterior)

            # arrange the combinations_probability_final by the mean_normalized_posterior in descending order
            combinations_probability_final <- combinations_probability_final[order(-combinations_probability_final$mean_normalized_posterior),]

            # if detail is TRUE, then save the combinations_probability_final to a csv file
            if (detail) {
                combinations_probability_final = as.data.frame(combinations_probability_final)
                write_csv(combinations_probability_final, paste0(MSMICA_decision_input_run_folder_name, "/", current_kegg_ids, "/", current_kegg_ids, "_bayesian_probability.csv"))
            }

            # select the rows with the highest mean_normalized_posterior
            highest_posterior_probability_adduct <- combinations_probability_final %>%
                filter(mean_normalized_posterior == max(mean_normalized_posterior)) %>%
                ## remove the adduct_formation_probability, adduct_correlation_probability, isotope_correlation_probability, enzyme_correlation_probability, mean_raw_posterior, and mean_normalized_posterior columns
                select(-enzyme_correlation_probability, -mean_raw_posterior, -mean_normalized_posterior) %>%
                ## transform the data frame into a long format with the col_names_annotated as the column name and the value as the value
                pivot_longer(cols = everything(), names_to = "col_names_annotated", values_to = "MSMICA_identification")

            # left join the highest_posterior_probability_adduct with the original data to get the final result using the col_names_annotated as the key
            final_result <- left_join(data, highest_posterior_probability_adduct, by = "col_names_annotated", relationship = "many-to-many")

            # add the identification_method column with value "precursor-product correlation"
            final_result$identification_method <- "precursor-product correlation"

            return(final_result)
        }

        # Apply tthe bayesian_probability_calculation_enzyme function to each group
        results <- lapply(data_split, function(x) bayesian_probability_calculation_enzyme(x))

        # Combine all data frames in the list into one final data frame
        final_results <- bind_rows(results)

        # if there is no new metabolite identified via enzyme based reaction on Kegg, then break the loop
        if (nrow(final_results) == 0) {
            print("There is no more new metabolite that can be identified using enzyme-based reaction correlation ranking. The MSMICA algorithm has finished.")
            break
        }

        # select only the rows with MSMICA_identification == 1
        final_results <- final_results %>%
            filter(MSMICA_identification == 1)

        # only if there is any new metabolite identified, then update the MSMICA_confidence column in the identified_metabolite
        if (nrow(final_results) > 0) {
            # print the number of new metabolites identified
            print(paste0("There are ", length(unique(final_results$connection_2)), " new metabolites identified"))

            # select specific columns from final_results: col_names_annotated, identification_type_annotated, ion_mode, identified_Name, Adduct_annotated, mz_annotated, time_annotated, correlation, MSMICA_identification, identification_method
            final_results_2 = final_results %>%
                dplyr::select(col_names_annotated, identification_type_annotated, ion_mode, identified_Name, connection_2, Adduct_annotated, mz_annotated, time_annotated, correlation, MSMICA_identification, identification_method) %>%
                ## rename columns: col_names_annotated as col_names, identification_type_annotated as identification_type, connection_2 as KEGGID, Adduct_annotated as Adduct, mz_annotated as mz, time_annotated as time
                rename(col_names = col_names_annotated, identification_type = identification_type_annotated, KEGGID = connection_2, Adduct = Adduct_annotated, mz = mz_annotated, time = time_annotated)

            # if detail = TRUE, then save the final_results_1 as a csv file temporarily
            if (detail == TRUE) {
                final_results_1_file_name <- paste0(new_identified_folder_name, "/", "new_identified_metabolite", "_", a, ".csv")
                write_csv(final_results_2, final_results_1_file_name)
            }

            print("This is the new identified metabolite:")
            print(final_results_2)

            # append the final_results_1 to the MSMICA_TCA_col_names_identified
            identified_metabolite = rbind(identified_metabolite, final_results_2) %>%
                arrange(desc(identification_type), KEGGID) %>%
                tibble()

            # dplyr::select only the columns with "annotated" in the identification_type column and KEGGID not in the identified_metabolite$KEGGID
            ## this removes the annotated metabolites that have been identified by MSMICA based on KEGGID
            annotated_metabolite = annotated_metabolite %>%
                filter(identification_type == "annotated" & !(KEGGID %in% identified_metabolite$KEGGID)) %>%
                tibble()
        }


        # remove the original cor_MSMICA to save memory
        rm(cor_MSMICA)

        # if all the m/z 5 ppm annotated metabolites have been identified by MSMICA, then break the loop
        if(nrow(annotated_metabolite) == 0) {
            print("All the annotated metabolites have been identified by MSMICA. The MSMICA algorithm has finished.")
            break
        }
        # if there is no new metabolite identified via enzyme based reaction on Kegg, then break the loop
        if (nrow(final_results) == 0) {
            print("There is no more new metabolite that can be identified using enzyme-based reaction. The MSMICA algorithm has finished.")
            break
        }
    }

    # After the loop
    close(pb)

    # add mean intensity values to all confirmed and identified metabolite by joining the mean_MSMICA_2
    identified_metabolite = identified_metabolite %>%
        left_join(mean_MSMICA_2, by = c("col_names" = "col_names"), relationship = "many-to-many") %>%
        ## move mean_intensity column right after rt
        relocate(mean_intensity, .after = "time")

    # select the metabolite adduct with the highest mean intensity for each metabolite, if there are multiple adducts for the same metabolite
    identified_metabolite = identified_metabolite %>%
        group_by(KEGGID) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup()

    # replace all "annotated" in the identification_type column with "identified"
    identified_metabolite = identified_metabolite %>%
        mutate(identification_type = ifelse(identification_type == "annotated", "identified", identification_type))

    # remove duplicates
    identified_metabolite = identified_metabolite %>%
        distinct(KEGGID, .keep_all = TRUE)

    ######################## to resolve some duplicated metabolites for some features, we use enzyme-based reaction correlation to filter out those metabolites identified by adduct correlation or isotope correlation but do not have any significant correlations with their precursor/product metabolites.
    # Perform a left join of the kegg_connection database to the identified_metabolite database
    # Allowing cartesian join by creating all possible combinations of rows
    MSMICA_col_names_connection <- identified_metabolite %>%
        left_join(kegg_connection, by = c("KEGGID" = "connection_1"), relationship = "many-to-many")

    # Arrange the data by the connection_2 column in ascending order
    MSMICA_col_names_connection_simplifed_2 <- MSMICA_col_names_connection %>%
        arrange(connection_2)

    # reorganize the columns in MSMICA_col_names_connection_simplifed: col_names, identification_type, KEGGID, Adduct, mz, time, MSMICA_identification, identification_method, connection_2
    MSMICA_col_names_connection_simplifed_2 = MSMICA_col_names_connection_simplifed_2[, c("col_names", "identification_type", "KEGGID", "Adduct", "mz", "time", "MSMICA_identification","identification_method", "connection_2")]

    # inner join the MSMICA_TCA_col_names_connection with the annotated data
    MSMICA_col_names_connection_identified_final = inner_join(
        MSMICA_col_names_connection_simplifed_2, 
        identified_metabolite, 
        by = c("connection_2" = "KEGGID"),
        suffix = c("_identified", "_identified_final"),
        relationship = "many-to-many"
        )
    
    ############ calculation correlation coefficients and p-values only for the necessary identified and annotated metabolites ############
    # write a loop to compare the identified metabolites with the annotated metabolites by using the correlation test and create a data frame with the correlation values and p-values
    # create a data frame to store the correlation values and p-values
    cor_MSMICA_data = data.frame(
        col_names_identified = character(),
        col_names_identified_final = character(),
        correlation = numeric(),
        p_value = numeric()
    )

    ## loop through each row of the MSMICA_col_names_connection_identified_final to perform correlation test between col_names_identified and col_names_identified_final by using subset the corresponding columns from MSMICA_cor_input
    for (i in 1:nrow(MSMICA_col_names_connection_identified_final)){
        current_row = MSMICA_col_names_connection_identified_final[i, ]

        identified_metabolite_current = current_row$col_names_identified
        annotated_metabolite_current = current_row$col_names_identified_final

        # check if identified_metabolite_current or annotated_metabolite_current are NA. If so, pass
        if (is.na(identified_metabolite_current) | is.na(annotated_metabolite_current)){
            next
        }

        identified_matrix <- unlist(MSMICA_cor_input[, identified_metabolite_current])
        other_matrix <- unlist(MSMICA_cor_input[, annotated_metabolite_current])

        cor_MSMICA = cor(identified_matrix, other_matrix, method="spearman", use="complete.obs")
        cor_p_MSMICA <- suppressWarnings(cor.test(identified_matrix, other_matrix, method = "spearman", use = "complete.obs"))

        # create a data frame with col_names_identified, col_names_identified_final, correlation, p_value
        cor_MSMICA_data_current = data.frame(
            col_names_identified = identified_metabolite_current,
            col_names_identified_final = annotated_metabolite_current,
            correlation = cor_MSMICA,
            p_value = cor_p_MSMICA$p.value
        )

        # combine the cor_MSMICA_data with cor_MSMICA_data_1
        cor_MSMICA_data = rbind(cor_MSMICA_data, cor_MSMICA_data_current)
    }

    # make correlation values absolute
    cor_MSMICA_data$correlation = abs(cor_MSMICA_data$correlation)

    # filter out those with p_value > 0.05
    cor_MSMICA_1 = cor_MSMICA_data %>%
        filter(p_value < 0.05)

    # remove the correlation column from MSMICA_col_names_connection_identified_final
    MSMICA_col_names_connection_identified_final = MSMICA_col_names_connection_identified_final %>%
        select(-correlation)

    # left join the cor_MSMICA_1 to the MSMICA_col_names_connection_identified_final to add correlation values
    ## now, we are basically adding the correlation values between the identified metabolites and the annotated metabolites from the previously created correlation matrix to the MSMICA_col_names_connection_identified_final
    ## here, we just know that the current row of the MSMICA_col_names_connection_identified_final is connected to all the identified metabolites, but we don't know which identified metabolite is connected to which annotated metabolite based on enzyme-based reaction on KEGG
    ## this is the first step to allow us to see the correlation values between the identified metabolites and the annotated metabolites for enzyme-based reaction correlation ranking
    MSMICA_col_names_connection_identified_final_2 = left_join(MSMICA_col_names_connection_identified_final, cor_MSMICA_1, by = c("col_names_identified_final", "col_names_identified"), relationship = "many-to-many")

    # arrange the correlation from high to low by each KEGGID and connection_2
    MSMICA_col_names_connection_identified_final_2 = MSMICA_col_names_connection_identified_final_2 %>%
        arrange(KEGGID, connection_2, desc(correlation))

    # arrange the data by connection_2 and correlation from high to low
    MSMICA_col_names_connection_identified_final_2 = MSMICA_col_names_connection_identified_final_2 %>%
        arrange(connection_2, desc(correlation))
    
    # round the correlation to 4 decimal places
    MSMICA_col_names_connection_identified_final_2$correlation = round(MSMICA_col_names_connection_identified_final_2$correlation, 4)

    # Filter out those self-correlation values (features with correlation of 1 with themselves)
    ## if correlation == 1, make it NA
    MSMICA_col_names_connection_identified_final_2 = MSMICA_col_names_connection_identified_final_2 %>%
        mutate(correlation = ifelse(correlation == 1, NA, correlation))

    # Summarize the correlation to mean_correlation by connection_2, col_names_identified_final, identification_type_identified_final, ion_mode, identified_Name, Adduct_identified_final, mz_identified_final, time_identified_final
    MSMICA_col_names_connection_identified_final_3 = MSMICA_col_names_connection_identified_final_2 %>%
        group_by(connection_2, col_names_identified_final, identification_type_identified_final, ion_mode, identified_Name, Adduct_identified_final, mz_identified_final, time_identified_final, MSMICA_identification_identified_final, identification_method_identified_final) %>%
        summarize(correlation = mean(correlation, na.rm = TRUE), .groups = "keep") %>%
        ungroup()
    
    # round the correlation to 4 decimal places
    MSMICA_col_names_connection_identified_final_3$correlation = round(MSMICA_col_names_connection_identified_final_3$correlation, 4)

    # filter out those with all NA correlation values
    MSMICA_col_names_connection_identified_final_3 = MSMICA_col_names_connection_identified_final_3 %>%
        filter(!is.na(correlation))
        
    # select specific columns from final_results: col_names_annotated, identification_type_annotated, ion_mode, identified_Name, Adduct_identified_final, mz_identified_final, time_identified_final, correlation, MSMICA_identification_identified_final, identification_method_identified_final
    MSMICA_col_names_connection_identified_final_4 = MSMICA_col_names_connection_identified_final_3 %>%
        dplyr::select(col_names_identified_final, ion_mode, identification_type_identified_final, identified_Name, connection_2, Adduct_identified_final, mz_identified_final, time_identified_final, correlation, MSMICA_identification_identified_final, identification_method_identified_final) %>%
        ## rename columns: col_names_annotated as col_names, identification_type_annotated as identification_type, connection_2 as KEGGID, Adduct_identified_final as Adduct, mz_identified_final as mz, time_identified_final as time
        rename(col_names = col_names_identified_final, identification_type = identification_type_identified_final, KEGGID = connection_2, Adduct = Adduct_identified_final, mz = mz_identified_final, time = time_identified_final, MSMICA_identification = MSMICA_identification_identified_final, identification_method = identification_method_identified_final) %>%
        ## remove the col_names column
        dplyr::select(-col_names)

    # select only KEGGID, Name, Exact_mass, and isotopic_mass columns from KEGG_database
    KEGG_database_simple = KEGG_database %>%
        dplyr::select(KEGGID, Name, Exact_mass, isotopic_mass) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1)) %>%
        ## rename Name as Identified_Name
        rename(Identified_Name = Name)
    
    # left join the KEGG_database_simple to the identified_metabolite to add the Identified_Name, Exact_mass, and isotopic_mass columns using KEGGID
    MSMICA_col_names_connection_identified_final_4 = left_join(MSMICA_col_names_connection_identified_final_4, KEGG_database_simple, by = c("KEGGID" = "KEGGID"), relationship = "many-to-many") %>%
        ## move Identified_Name column right after KEGGID
        relocate(Identified_Name, .after = "KEGGID") %>%
        ## move Exact_mass column right after Identified_Name
        relocate(Exact_mass, .after = "Identified_Name") %>%
        ## move isotopic_mass column right after Exact_mass
        relocate(isotopic_mass, .after = "Exact_mass") %>%
        ## move ion_mode before identification_type
        relocate(ion_mode, .before = "identification_type") %>%
        ## remove the identified_Name column
        dplyr::select(-identified_Name)
    
    MSMICA_col_names_connection_identified_final_4$mz = as.numeric(MSMICA_col_names_connection_identified_final_4$mz)
    MSMICA_col_names_connection_identified_final_4$time = as.numeric(MSMICA_col_names_connection_identified_final_4$time)

    if (prefix != ""){
        # add the prefix to saved file names
        write_csv(MSMICA_col_names_connection_identified_final_4, paste0(prefix, "_MSMICA_identified_metabolites.csv"))
        write_csv(annotated_metabolite, paste0(prefix, "_MSMICA_not_identified_metabolites.csv"))
    } else {
        # save the data without prefix
        write_csv(MSMICA_col_names_connection_identified_final_4, "identified_MSMICA_identified_metabolites.csv")
        write_csv(annotated_metabolite, "MSMICA_not_identified_metabolites.csv")
    }

    # print where the MSMICA_TCA_result_identified.csv and MSMICA_TCA_result_not_identified.csv are saved
    print(paste0("The identified_MSMICA_identified_metabolites.csv was saved in ", getwd()))
}