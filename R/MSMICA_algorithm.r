#' MSMICA algorithm
#' 
#' This function is used to perform the MSMICA algorithm for metabolite identification using the metabolomics feature table and the KEGG database.
#' @rawNamespace import(dplyr, except = c(first, last, between))
#' @param met_raw_wide a metabolomics feature table in wide format with mz as the first column, time as the second column, and intensity values as the remaining columns.
#' @param class_file a class file in wide format with the first column name as metabolomics raw file name, the second column name as subject ID, and the third column name as class label (study sample or reference standard sample).
#' @param LC a character value indicating which liquid chromatography (LC) column to be used to predict the retention time of the metabolites. Default is "HILIC" (hydrophilic interaction liquid chromatography). Other options is "RP" or "C18" (reversed phase liquid chromatography, also called C18).
#' @param LC_run_time a numeric value indicating the run time of the liquid chromatography. Default is 5 minutes.
#' @param mz_threshold the m/z threshold for the metabolite identification. Default is 10 ppm.
#' @param biospecimen a character value indicating the biospecimen of the study samples. Default is "Blood". The other options include "Urine", "Feces", "Cerebrospinal Fluid", "Saliva", "Breast Milk", "Sweat", "Cellular Cytoplasm", "Amniotic Fluid", "Aqueous Humour", "Ascites Fluid", "Lymph", "Tears",  "Bile", "Semen", "Pericardial Effusion"
#' @param hmdb_detection_preference a logical value indicating whether to use the HMDB detection preference for the adduct identification. Default is TRUE. If TRUE, then only the metabolites noted as "detected" in the HMDB database will be used for MSMICA algorithm. If FALSE, then all the metabolites specified in the metabolite_database will be used for MSMICA algorithm.
#' @param All_Adduct the adduct forms of the metabolites. Default is c("M+H","M+Na","M+2Na-H","M+H-H2O","M+H-NH3","M+ACN+H","M+ACN+2H","2M+H","M+2H","M+H-2H2O") for the positive mode. This includes primary and secondary adducts.
#' @param metabolite_database a character value indicating the metabolite database to be used. Default is "KEGG_HMDB".
#' @param reaction_database a character vector specifying the reaction database to be used. Default is c("mammalia"). Other option is c("general").
#' @param backpropagation_correlation_direction the direction of the backpropagation precursor-product/transportercorrelation coefficient. Default is "positive". Other options are "both". If positive, then only the positive backpropagation correlation coefficient is used. If both, then both the positive and negative backpropagation correlation coefficients are used.
#' @param adduct_correlation_r_threshold the correlation threshold for adduct correlation analysis. Default is 0.39 (spearman correlation).
#' @param adduct_correlation_time_threshold the retention time threshold for adduct correlation analysis. Default is 6 (seconds).
#' @param isotopic_correlation_r_threshold the correlation threshold for isotopic correlation analysis. Default is 0.71 (spearman correlation).
#' @param isotopic_correlation_time_threshold the retention time threshold for isotopic correlation analysis. Default is 4 (seconds).
#' @param imputation_method the method to be used for missing value imputation. Default is "half_min". Other options are "QRILC" (QRILC is better for triplicate samples) and NA. If NA, then no imputation is performed.
#' @param prefix a prefix to be added to the output files. Default is "".
#' @param ion_mode the ionization mode of the metabolomics data. Default is "positive". Other options are "negative".
#' @param detail a logical value indicating whether to save the intermediate results as csv files. Default is FALSE. WARNING, this can create thousands of files with a lot of space. Use with caution.
#' @param save_unidentified a logical value indicating whether the unidentified features should be saved. Default is FALSE.
#' @param progress_log a logical value indicating whether to save the log of all printings and messages to a text file. Default is TRUE.
#' @export MSMICA_algorithm

MSMICA_algorithm = function(met_raw_wide, class_file = NULL, LC = "HILIC", LC_run_time, mz_threshold = 10, biospecimen = "Blood", hmdb_detection_preference = TRUE,  All_Adduct = c("M+H","M+Na","M+2Na-H","M+H-H2O","M+H-NH3","M+ACN+H","M+ACN+2H","2M+H","M+2H","M+H-2H2O"), metabolite_database = "KEGG_HMDB", reaction_database = c("mammalia"), backpropagation_correlation_direction = "positive", adduct_correlation_r_threshold = 0.39, adduct_correlation_time_threshold = 6, isotopic_correlation_r_threshold = 0.71, isotopic_correlation_time_threshold = 4, imputation_method = "half_min", prefix = "", ion_mode = "positive", detail = FALSE, save_unidentified = FALSE, progress_log = FALSE) {
    # Load all required packages only if necessary
    library(dplyr)
    library(readr)
    library(tidyr)
    library(data.table)
    library(mgcv)
    library(pracma)
    library(MetaboCoreUtilsAdduct)
    # imputeLCMD is not needed to be loaded here because we only use a few functions from it
    # data.table is not needed to be loaded here because find.Overlapping.mzs will load it later
    # preprocessCore is not needed to be loaded here because we only use a few functions from it

    # if progress_log = TRUE, then establish a log file (txt format) to store all the printing and messages
    if (progress_log) {
        log_file = file("MSMICA_algorithm_log.txt", open = "wt")

        sink(log_file, type = "output")
        sink(log_file, type = "message")

        on.exit({
            sink(type = "message")
            sink(type = "output")
            close(log_file)
        }, add = TRUE)
    }

    # log the start time of the MSMICA algorithm
    start_time = Sys.time()
    message("The MSMICA algorithm starts at: ", start_time)

    # make sure the mz and time are numeric in the met_raw_wide
    met_raw_wide$mz = as.numeric(met_raw_wide$mz)
    met_raw_wide$time = as.numeric(met_raw_wide$time)

    # remove the rows with NA in the mz or time column in the met_raw_wide
    met_raw_wide = met_raw_wide %>%
        dplyr::filter(!is.na(mz) & !is.na(time))

    # make sure LC is within the possible options
    if (!LC %in% c("HILIC", "RP", "C18")) {
        stop("The LC is not within the possible options. The possible options are: HILIC, RP, C18. Please check the LC parameter.")
    }

    ## if LC_run_time is higher than 100, put a warning message because it is not realistic
    if (LC_run_time > 100) {
        stop("The LC run time is higher than 100 minutes. This is not realistic. The LC_run_time parameter is based on minutes, not seconds. Please check the LC_run_time parameter.")
    }

    # if the class_file is not NULL, then rename the class_file column name as file_name, subject_id, and class_label
    if (!is.null(class_file)) {

        # select only the first 3 columns of the class_file
        class_file = class_file %>%
            select(1:3)

        # rename the column names of the class_file
        colnames(class_file) = c("file_name", "subject_id", "class_label")

        # find the sample with the class_label as "study"
        study_sample = class_file %>%
            filter(class_label == "study") %>%
            select(file_name) %>%
            pull()

        # clean up the column names of the met_raw_wide: if it has ".mzXML", ".mzML", ".raw", ".cdf" remove it
        colnames(met_raw_wide) = gsub(".mzXML", "", colnames(met_raw_wide))
        colnames(met_raw_wide) = gsub(".mzML", "", colnames(met_raw_wide))
        colnames(met_raw_wide) = gsub(".raw", "", colnames(met_raw_wide))
        colnames(met_raw_wide) = gsub(".cdf", "", colnames(met_raw_wide))
        
        # remove study sample that is not in the met_raw_wide's column names
        study_sample = study_sample[study_sample %in% colnames(met_raw_wide)]

        # keep only the study samples in the met_raw_wide's column names
        met_raw_wide = met_raw_wide %>%
            select(mz, time, all_of(study_sample))
        
        print("Only the study samples are kept in the feature table.")
    }

    ##################### prepare the hmdb metabolite concentration database #####################
    # import the hmdb_metabolites_concentrations_average data frame
    data(hmdb_metabolites_concentrations_average)
    message("Biological concentration prior knowledge from HMDB is used.")
    # check if the biospecimen is within the possible options
    if (!biospecimen %in% hmdb_metabolites_concentrations_average$Biospecimen) {
        stop("The biospecimen is not within the possible options. The possible options are: ", paste(unique(hmdb_metabolites_concentrations_average$Biospecimen), collapse = ", "))
    }
    # select the indicated Biospecimen type 
    hmdb_metabolites_concentrations_average = hmdb_metabolites_concentrations_average %>%
        filter(Biospecimen == biospecimen)
    message("The selected biospecimen is: ", biospecimen)
    hmdb_metabolites_concentrations_average = hmdb_metabolites_concentrations_average %>%
        select(HMDBID, KEGG_ID, InChIKey, Name, Monisotopic_molecular_weight, Concentration_average, Concentration_sd, Concentration_units)
    # keep only the first 14 characters of the InChIKey in the hmdb_metabolites_concentrations_average
    hmdb_metabolites_concentrations_average$InChIKey = substr(hmdb_metabolites_concentrations_average$InChIKey, 1, 14)
    # group by InChIKey and keep the rows with the highest Concentration_average
    hmdb_metabolites_concentrations_average = hmdb_metabolites_concentrations_average %>%
        group_by(InChIKey) %>%
        filter(Concentration_average == max(Concentration_average)) %>%
        ungroup() %>%
        distinct(InChIKey, .keep_all = TRUE)
    # simplify the hmdb_metabolites_concentrations_average data frame
    hmdb_metabolites_concentrations_average_simple = hmdb_metabolites_concentrations_average %>%
        dplyr::select(InChIKey, Concentration_average, Concentration_units)

    ##################### prepare the metabolite database #####################
    # import the specified metabolite database metabolite_database
    ## if metabolite_database is "KEGG_HMDB", then import the KEGG_HMDB_database_mainchiral
    if (metabolite_database == "KEGG_HMDB") {
        data(KEGG_HMDB_database_mainchiral)
        metabolite_database = KEGG_HMDB_database_mainchiral
        message("KEGG_HMDB compound database with main chiral preference is used.")
    } 
    # keep the original InChIKey in the metabolite_database by selecting Name, PubChem_compound_id, HMDB_ID, KEGG_ID, and InChIKey, and using the concat of Name, PubChem_compound_id, HMDB_ID, KEGG_ID as the unique identifier to keep the original InChIKey
    metabolite_database_original_InChIKey = metabolite_database %>% 
        select(Name, PubChem_compound_id, HMDB_ID, KEGG_ID, InChIKey) %>%
        mutate(unique_identifier = paste0(Name, "_", PubChem_compound_id, "_", HMDB_ID, "_", KEGG_ID)) %>%
        select(unique_identifier, InChIKey)
    # keep only the first 14 characters of the InChIKey in the metabolite_database
    metabolite_database$InChIKey = substr(metabolite_database$InChIKey, 1, 14)
    # group by InChIKey and remove the duplicates based on InChIKey (the data is sorted by HMDB_ID first, then KEGG_ID. The first record will be kept and the rest will be removed)
    metabolite_database = metabolite_database %>%
        group_by(InChIKey) %>%
        filter(row_number() == 1) %>%
        ungroup()

    ##################### prepare the hmdb detection preference database #####################
    # if hmdb_detection_preference is TRUE, then load the hmdb_endogenous_metabolites_detected_and_quantified database for prioritization
    if (hmdb_detection_preference == TRUE) {
        data(hmdb_endogenous_metabolites_detected_and_quantified)

        # extract non-NA values from the HMDBID and InChIKey columns in the hmdb_endogenous_metabolites_detected_and_quantified dataframe
        hmdb_endogenous_metabolites_detected_and_quantified_HMDBID = hmdb_endogenous_metabolites_detected_and_quantified %>%
            filter(!is.na(HMDBID)) %>%
            select(HMDBID) %>%
            unique() %>%
            pull()

        hmdb_endogenous_metabolites_detected_and_quantified_InChIKey = hmdb_endogenous_metabolites_detected_and_quantified %>%
            filter(!is.na(InChIKey)) %>%
            select(InChIKey) %>%
            unique() %>%
            pull()
    }
    # keep only the first 14 characters of the InChIKey in the hmdb_endogenous_metabolites_detected_and_quantified
    hmdb_endogenous_metabolites_detected_and_quantified_InChIKey = substr(hmdb_endogenous_metabolites_detected_and_quantified_InChIKey, 1, 14)
    # remove duplicates
    hmdb_endogenous_metabolites_detected_and_quantified_InChIKey = unique(hmdb_endogenous_metabolites_detected_and_quantified_InChIKey)
    # append the records from hmdb_metabolites_concentrations_average
    hmdb_endogenous_metabolites_detected_and_quantified_HMDBID = unique(c(hmdb_endogenous_metabolites_detected_and_quantified_HMDBID, hmdb_metabolites_concentrations_average$HMDBID))
    hmdb_endogenous_metabolites_detected_and_quantified_InChIKey = unique(c(hmdb_endogenous_metabolites_detected_and_quantified_InChIKey, hmdb_metabolites_concentrations_average$InChIKey))

    ##################### prepare the precursor-product/transporter relationship database #####################
    # import the precursor-product reaction database according to the reaction_database input for the initial propagation of precursor-product correlation
    if (reaction_database == "general") {
        message("Generic biochemical reaction from Recon3D, BKMS-react, rhea, and KEGG database is used.")
        data(Recon3D_BKMS_react_rhea_KEGG_connection_global)
        reaction_connection = Recon3D_BKMS_react_rhea_KEGG_connection_global %>%
            select(connection_1_InChIKey, connection_2_InChIKey, react_id, enzyme, source)
    } else if (reaction_database == "mammalia") {
        message("Mammalia biochemical reaction from Recon3D, BKMS-react, rhea, and KEGG database is used.")
        data(Recon3D_BKMS_react_rhea_KEGG_connection_mammalia)
        reaction_connection = Recon3D_BKMS_react_rhea_KEGG_connection_mammalia %>%
            select(connection_1_InChIKey, connection_2_InChIKey, react_id, enzyme, source)
    }

    # load the custom biochemical reaction dataset
    custom_biochemical_reaction = custom_biochemical_reaction_loading()
    # select the connection_1_InChIKey, connection_2_InChIKey, react_id, enzyme, source columns
    custom_biochemical_reaction = custom_biochemical_reaction %>%
        select(connection_1_InChIKey, connection_2_InChIKey, react_id, enzyme, source)
    message("Custom biochemical reaction database is used.")
    # copy the custom_biochemical_reaction and flip the connection_1_InChIKey and connection_2_InChIKey
    custom_biochemical_reaction_flipped = custom_biochemical_reaction %>%
        mutate(connection_1_InChIKey = connection_2_InChIKey, connection_2_InChIKey = connection_1_InChIKey)
    # combine the custom_biochemical_reaction and custom_biochemical_reaction_flipped
    custom_biochemical_reaction = bind_rows(custom_biochemical_reaction, custom_biochemical_reaction_flipped)

    # import the human transporter database
    data(Recon3D_unitprot_Deo_human_transporter)
    message("Human transporter database is used.")
    # group by Enzyme_Transporter and keep only the rows with more than 1 Enzyme_Transporter
    Recon3D_unitprot_Deo_human_transporter = Recon3D_unitprot_Deo_human_transporter %>%
        group_by(Enzyme_Transporter) %>%
        filter(n() > 1) %>%
        ungroup()

    # expand Recon3D_unitprot_Deo_human_transporter into a long format: group by Enzyme_Transporter, and then create a dataframe that with different combination of Metabolite_inchikey (one as Connection_1_InChIKey and the other as Connection_2_InChIKey)
    Recon3D_unitprot_Deo_human_transporter_long = Recon3D_unitprot_Deo_human_transporter %>%
        # Select only relevant columns to avoid duplicated metadata
        select(source, Enzyme_Transporter, Metabolite_inchikey) %>%
        # Join the table to itself by the transporter
        inner_join(., ., by = "Enzyme_Transporter", relationship = "many-to-many") %>%
        # Rename for clarity
        rename(connection_1_InChIKey = Metabolite_inchikey.x, 
                connection_2_InChIKey = Metabolite_inchikey.y, source = source.x) %>%
        # Optional: Remove rows where a metabolite is paired with itself (e.g., AA, BB)
        filter(connection_1_InChIKey != connection_2_InChIKey) %>%
        mutate(react_id = NA) %>%
        select(connection_1_InChIKey, connection_2_InChIKey, react_id,Enzyme_Transporter, source)
    
    # unify the fourth column name as enzyme_transporter
    colnames(Recon3D_unitprot_Deo_human_transporter_long)[4] = "enzyme_transporter"
    colnames(reaction_connection)[4] = "enzyme_transporter"
    colnames(custom_biochemical_reaction)[4] = "enzyme_transporter"

    # combine the reaction_connection and custom_biochemical_reaction
    reaction_connection = bind_rows(reaction_connection, custom_biochemical_reaction, Recon3D_unitprot_Deo_human_transporter_long) 
    
    # keep only the first 14 characters of the InChIKey in the reaction_connection
    reaction_connection$connection_1_InChIKey = substr(reaction_connection$connection_1_InChIKey, 1, 14)
    reaction_connection$connection_2_InChIKey = substr(reaction_connection$connection_2_InChIKey, 1, 14)

    # remove duplicates
    reaction_connection = reaction_connection %>%
        distinct(connection_1_InChIKey, connection_2_InChIKey, .keep_all = TRUE)
    
    # if detail = TRUE, then create the following folders
    if (detail == TRUE){
        # Create a master folder with the study name to store all the detailed results
        master_folder_name = paste0(prefix, "_", mz_threshold, "ppm")
        dir.create(master_folder_name, showWarnings = FALSE)

        # Create a folder to store the result: m/z matching for all adducts with or without isotopes
        mz_matching_folder_name = paste0(master_folder_name, "/step1_mz_matching", "_", mz_threshold, "ppm")
        dir.create(mz_matching_folder_name, showWarnings = FALSE)

        # Create a folder to store the results: rt mapping anchor metabolites between current study and retention time reference library
        rt_mapping_anchor_metabolites_folder_name = paste0(master_folder_name, "/step2_algorithm_parameter_estimation", "_", mz_threshold, "ppm")
        dir.create(rt_mapping_anchor_metabolites_folder_name, showWarnings = FALSE)

        # Create a folder to store the result: Clustering of adducts and isotopes
        clustering_folder_name = paste0(master_folder_name, "/step3_clustering", "_", mz_threshold, "ppm")
        dir.create(clustering_folder_name, showWarnings = FALSE)

        # Create a folder to store result: local optimization per feature
        local_optimization_per_feature_folder_name = paste0(master_folder_name, "/step4_local_optimization_per_feature", "_", mz_threshold, "ppm")
        dir.create(local_optimization_per_feature_folder_name, showWarnings = FALSE)
    }

    ##################### prepare the metabolite database for m/z matching #####################
    # Cross join All_Adduct with metabolite_database to form a adduct list for each unique KEGG_ID (primary and secondary adducts)
    metabolite_database_2 = tidyr::crossing(metabolite_database, All_Adduct)

    # rename All_Adduct as Adduct
    metabolite_database_2 = metabolite_database_2 %>%
        dplyr::rename(Adduct = All_Adduct)

    # remove the following adducts if the following conditions are met. This helps to reduce the number of adducts to be considered for the m/z matching given the metabolite structure.
    # M+H-NH3 and has_NH3 is FALSE
    # M-H2O+H or M+H-H2O or M-H2O-H or M-H-H2O and n_H2O >= 1
    # M-2H2O+H or M+H-2H2O and n_H2O >= 2
    metabolite_database_2 = metabolite_database_2 %>%
        # Remove rows based on adduct-specific conditions and include the adducts that are not affected by the conditions
        filter(
            (Adduct == "M+H-NH3" & has_NH3 == TRUE) |
            (Adduct %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O") & n_H2O >= 1) |
            (Adduct %in% c("M-2H2O+H", "M+H-2H2O") & n_H2O >= 2) |
            (!Adduct %in% c("M+H-NH3", "M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O", "M-2H2O+H", "M+H-2H2O"))
        ) %>%
        dplyr::select(-c(has_NH3, has_H2O, n_NH3, n_H2O))

    # calculate mz using Mono_mass and Adduct
    metabolite_database_2_mz = MetaboCoreUtilsAdduct::mass2mz_df(mass=metabolite_database_2$Mono_mass, adduct=metabolite_database_2$Adduct)
    # add mz to the metabolite_database_2
    metabolite_database_2$mz = metabolite_database_2_mz$mz

    # calculate the heavy isotope mz using Most_abundant_isotopologue_mass and and Adduct
    metabolite_database_2_mz_heavy = MetaboCoreUtilsAdduct::mass2mz_df(mass=metabolite_database_2$Most_abundant_isotopologue_mass, adduct=metabolite_database_2$Adduct)
    # add heavy isotope mz to the metabolite_database_2
    metabolite_database_2$mz_isotope = metabolite_database_2_mz_heavy$mz

    # remove mz less than 0
    metabolite_mz_library = metabolite_database_2 %>%
        dplyr::filter(mz > 0) 
    rm(metabolite_database_2, metabolite_database_2_mz, metabolite_database_2_mz_heavy)
    gc()

    # for met_raw_wide, replace all 0 values with NA
    met_raw_wide[met_raw_wide == 0] = NA

    # if there are duplicates in the met_raw_wide based on mz and time, keep the feature with the highest mean intensity
    met_raw_wide = met_raw_wide %>%
        # calculate the mean intensity values for each feature across all samples
        mutate(mean_intensity = rowMeans(select(., -mz, -time), na.rm = TRUE)) %>%
        # group by mz and time and select the feature with highest mean_intensity if there are duplicates. These are the duplicates after the mz and time rounding
        group_by(mz, time) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup()

    # extract the mz, time, and mean_intensity columns for met_raw_wide_original_mean_intensity
    met_raw_wide_original_mean_intensity = met_raw_wide %>%
        dplyr::select(mz, time, mean_intensity)

    # remove the mean_intensity column
    met_raw_wide = met_raw_wide %>%
        dplyr::select(-mean_intensity)

    # save the met_raw_wide as met_raw_wide_original for later use
    met_raw_wide_original = met_raw_wide

    ##################### perform metabolomicsdata transformation ##################
    # perform log2 transformation
    met_raw_wide = data_transformation(met_raw_wide = met_raw_wide)

    ##################### perform metabolomics data imputation #####################
    # use the function to imputate missing data in the metabolomics feature table (those 0 intensity data)
    met_raw_wide = data_imputation(met_raw_wide = met_raw_wide, imputation_method = imputation_method)
    met_raw_wide = tibble(met_raw_wide)

    ##################### perform metabolomics data normalization ##################
    # perform quantile normalization
    met_raw_wide = quantile_normalization(met_raw_wide = met_raw_wide)

    ##################### perform mz matching ######################################
    # use the mz_matching function to perform mz matching for all primary and secondary adducts
    met_raw_wide_final = mz_matching(met_raw_wide = met_raw_wide, metabolite_mz_library = metabolite_mz_library, mz_threshold = mz_threshold, detail = detail,
        mean_intensity_ref = met_raw_wide_original_mean_intensity, ion_mode = ion_mode,
        output_folder = if (detail) mz_matching_folder_name else NULL)

    # print the number of rows in the met_raw_wide_final
    message("Here is the m/z matching results for all primary and secondary adducts:")
    print(paste0("The m/z matching threshold is: ", mz_threshold, " ppm"))
    print(met_raw_wide_final)

    ##################### perform mz isotope matching ###############################
    # use the mz_isotope_matching function to perform mz matching for all adducts with heavy isotopes
    met_raw_wide_final_isotope = mz_isotope_matching(met_raw_wide = met_raw_wide, metabolite_mz_library = metabolite_mz_library, mz_threshold = mz_threshold, detail = detail,
        mean_intensity_ref = met_raw_wide_original_mean_intensity, ion_mode = ion_mode,
        output_folder = if (detail) mz_matching_folder_name else NULL)

    # print the number of rows in the met_raw_wide_final_isotope
    message("Here is the m/z matching results for all primary and secondary adducts' isotopologues by considering isotopic distribution. Only the most abundant isotopologue is used:")
    print(paste0("The m/z matching threshold is: ", mz_threshold, " ppm"))
    print(met_raw_wide_final_isotope)

    # round the mz to 4 decimal places
    met_raw_wide$mz = round(met_raw_wide$mz, 4)
    # round the time to 0 decimal places
    met_raw_wide$time = round(met_raw_wide$time, 0)

    # convert the feature table to the long format so that later functions can use it directly ######################
    ## add mz_time column by combining mz and time
    met_raw_long = met_raw_wide %>%
        mutate(mz_time = paste0(mz, "_", time)) %>%
        relocate(mz_time, .after = time)

    # Calculate the mean intensity values for each feature across all samples
    met_raw_long = met_raw_long %>%
        mutate(mean_intensity = rowMeans(select(., -mz, -time, -mz_time), na.rm = TRUE)) %>%
        # group by mz and time and select the feature with highest mean_intensity if there are duplicates. These are the duplicates after the mz and time rounding
        group_by(mz_time) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup() %>%
        select(-mean_intensity)

    ## transform the data from wide to long format
    MSMICA_cor_input = t(met_raw_long[,4:ncol(met_raw_long)])
    MSMICA_cor_input = as.data.frame(MSMICA_cor_input)
    ## add the mz_time column to the data frame
    colnames(MSMICA_cor_input) = met_raw_long$mz_time
    ## convert the data frame to a tibble
    MSMICA_cor_input = tibble(MSMICA_cor_input)

    # remove duplicates after mz and time rounding for typical adducts
    ## add InChIKey_mz_time column by combining InChIKey, mz (rounded to 4 decimals) and time (rounded to 1 decimal), separated by "_"
    ## group by InChIKey_mz_time and select the feature with highest mean_intensity if there are duplicates. 
    met_raw_wide_final_1 = met_raw_wide_final %>%
        mutate(
            InChIKey_mz_time = paste0(InChIKey, "_", mz_sample, "_", time_sample)
        )  %>%
        group_by(InChIKey_mz_time) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup() %>%
        ## remove the InChIKey_mz_time column
        dplyr::select(-InChIKey_mz_time)

    # remove duplicates after mz and time rounding for isotopic adducts
    ## add InChIKey_mz_time column by combining InChIKey, mz (rounded to 4 decimals) and time (rounded to 1 decimal), separated by "_"
    ## group by InChIKey_mz_time and select the feature with highest mean_intensity if there are duplicates. 
    met_raw_wide_final_isotope_1 = met_raw_wide_final_isotope %>%
        mutate(
            InChIKey_mz_time = paste0(InChIKey, "_", mz_sample_isotope, "_", time_sample_isotope)
        )  %>%
        group_by(InChIKey_mz_time) %>%
        filter(mean_intensity == max(mean_intensity)) %>%
        ungroup() %>%
        ## remove the InChIKey_mz_time column
        dplyr::select(-InChIKey_mz_time)

    # select the necessary columns for bayesian probability calculation for metabolite adduct clusters for typical adducts
    met_raw_wide_final_monomass = met_raw_wide_final_1 %>%
        # create a new column, identification_type_annotated, with value "identified"
        mutate(identification_type_annotated = "identified") %>%
        # rename Name to Confirmed_Name, Adduct to Adduct_annotated, mz_sample to mz_annotated, time_sample to time_annotated
        rename(Confirmed_Name = Name) %>%
        # only select the following columns: Mono_mass, InChIKey, identification_type_annotated, ion_mode, Confirmed_Name, Adduct, mz_sample, time_sample, mz_time_sample
        dplyr::select(Mono_mass, InChIKey, identification_type_annotated, ion_mode, Confirmed_Name, Adduct, mz_sample, time_sample, mz_time_sample) %>%
        # rename Adduct to Adduct_annotated, mz_sample to mz_annotated, time_sample to time_annotated
        rename(Adduct_annotated = Adduct, mz_annotated = mz_sample, time_annotated = time_sample, mz_time_annotated = mz_time_sample) %>%
        # round the Mono_mass to 4 decimals
        mutate(Mono_mass = round(Mono_mass, 4))

    # select the necessary columns for bayesian probability calculation for metabolite adduct clusters for isotopic adducts
    met_raw_wide_final_monomass_isotope = met_raw_wide_final_isotope_1 %>%
        # create a new column, identification_type_annotated, with value "identified"
        mutate(identification_type_annotated = "identified") %>%
        # rename Name to Confirmed_Name, Adduct to Adduct_annotated, mz_sample to mz_annotated, time_sample to time_annotated
        rename(Confirmed_Name = Name) %>%
        # only select the following columns: Mono_mass, InChIKey, identification_type_annotated, ion_mode, Confirmed_Name, Adduct, mz_sample_isotope, time_sample_isotope, mz_time_sample_isotope
        dplyr::select(Mono_mass, InChIKey, identification_type_annotated, ion_mode, Confirmed_Name, Adduct, mz_sample_isotope, time_sample_isotope, mz_time_sample_isotope) %>%
        # rename Adduct to Adduct_annotated, mz_sample to mz_annotated, time_sample to time_annotated
        rename(Adduct_annotated = Adduct, mz_annotated = mz_sample_isotope, time_annotated = time_sample_isotope, mz_time_annotated = mz_time_sample_isotope) %>%
        # round the Mono_mass to 4 decimals
        mutate(Mono_mass = round(Mono_mass, 4))

    print("Results are kept if either adduct correlation or isotopic correlation is possible")

    # simplify the data frame for typical adducts
    ## arrange by Adduct_annotated first and then Mono_mass
    met_raw_wide_final_monomass = met_raw_wide_final_monomass %>%
        arrange(Adduct_annotated) %>%
        arrange(Mono_mass)
    ## craete another data frame and remove InChIKey and Confirmed_Name columns. Then, remove duplicates
    ## here, we basically use the monoisotopic mass to present all the metabolites with the same monoisotopic mass. This greatly enhances the efficiency of the clustering algorithm because the clustering algorithm is based on the monoisotopic mass and there is no need to repeat the clustering algorithm with all the metabolites with the same monoisotopic mass.
    met_raw_wide_final_monomass_simple = met_raw_wide_final_monomass %>%
        dplyr::select(-InChIKey, -Confirmed_Name) %>%
        distinct()

    # simplify the data frame for isotopic adducts
    ## arrange by Adduct_annotated first and then Mono_mass
    met_raw_wide_final_monomass_isotope = met_raw_wide_final_monomass_isotope %>%
        arrange(Adduct_annotated) %>%
        arrange(Mono_mass)
    ## craete another data frame and remove InChIKey and Confirmed_Name columns. Then, remove duplicates
    ## here, we basically use the monoisotopic mass to present all the metabolites with the same monoisotopic mass.
    met_raw_wide_final_monomass_isotope_simple = met_raw_wide_final_monomass_isotope %>%
        dplyr::select(-InChIKey, -Confirmed_Name) %>%
        distinct()

    # select only the Mono_mass, Name, KEGG_ID, HMDB_ID, InChIKey, and Subclass columns
    metabolite_database_simple = metabolite_database %>%
        dplyr::select(Mono_mass, Name, KEGG_ID, HMDB_ID, InChIKey, Subclass) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1)) %>%
        ## round the Mono_mass to 4 decimals
        mutate(Mono_mass = round(Mono_mass, 4))
    
    message("Using metabolites with clean chromatography and no isomers (only 1 feature in the feature table with the same m/z) for predicting the retention time of the metabolites in the current study by mapping the current study retention timewith the retention time reference library.")

    # create a data frame to store the subclass information for the metabolites
    metabolite_database_simple_subclass = metabolite_database %>%
        dplyr::select(InChIKey, Subclass) %>%
        distinct()
    
    # within the m/z matching results, find where there is only 1 row for each Mono_mass. These are the metabolites with only 1 feature in the feature table with similar m/z (very clean chromatography and no isomers).
    met_raw_wide_mz_single_match = met_raw_wide_final_monomass %>%
        ## inner join with hmdb_metabolites_concentrations_average_simple based on InChIKey
        inner_join(hmdb_metabolites_concentrations_average_simple, by = "InChIKey") %>%
        ## inner join with the metabolite_database_simple_subclass based on InChIKey
        inner_join(metabolite_database_simple_subclass, by = "InChIKey") %>%
        ## filter the rows where:
        ### the Adduct_annotated is in c("M-H","M+H") and Subclass is not "Carbohydrates and carbohydrate conjugates"
        ### the Adduct_annotated is in c("M+Na","M+K","M+Cl") and Subclass is "Carbohydrates and carbohydrate conjugates"
        ### the Adduct_annotated is in c("M+H-H2O", "M+H-2H2O") and Subclass is "Bile acids, alcohols and derivatives"
        ### the Adduct_annotated is in c("M+H-H2O") and Subclass is "Retinoids"
        ### the Adduct_annotated is in c("M+H-H2O", "M+H-2H2O") and Subclass is "Pregnane steroids"
        dplyr::filter(
            Adduct_annotated %in% c("M-H","M+H")| 
            (Adduct_annotated %in% c("M+Na","M+K","M+Cl") & Subclass == "Carbohydrates and carbohydrate conjugates") | 
            (Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & Subclass == "Bile acids, alcohols and derivatives") | 
            (Adduct_annotated %in% c("M+H-H2O") & (Subclass == "Retinoids")) | 
            (Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & Subclass == "Pregnane steroids")
            ) %>%
        ## group by Mono_mass, Confirmed_Name, Adduct_annotated, and Subclass
        group_by(Mono_mass, Confirmed_Name, Adduct_annotated) %>%
        # remove the rows where there are multiple features with the same Mono_mass
        filter(n() == 1) %>%
        ungroup() %>%
        dplyr::select(-Subclass)
        
    # add KEGG_ID and HMDB_ID to the met_raw_wide_mz_single_match by left joining with the metabolite_database_simple based on Mono_mass, Confirmed_Name=Name, and InChIKey
    met_raw_wide_mz_single_match = left_join(met_raw_wide_mz_single_match, metabolite_database_simple, by = c("Mono_mass", "Confirmed_Name"="Name", "InChIKey"), relationship = "many-to-many")
    
    # arrange by mz_annotated and time_annotated
    met_raw_wide_mz_single_match = met_raw_wide_mz_single_match %>%
        arrange(mz_annotated, time_annotated)
    
    # Add identification_method as "mz matching" to the met_raw_wide_mz_single_match
    met_raw_wide_mz_single_match$identification_method = "mz matching"

    # Add match_category = "single". This indicates that the metabolite is matched with only 1 feature in the feature table.
    met_raw_wide_mz_single_match$match_category = "single"

    # rename the columns: mz_annotated to mz, time_annotated to time, mz_time_annotated to mz_time, Confirmed_Name to identified_Name
    met_raw_wide_mz_single_match = met_raw_wide_mz_single_match %>%
        rename(mz = mz_annotated, time = time_annotated, mz_time = mz_time_annotated, identified_Name = Confirmed_Name)

    if (detail == TRUE) {
        # save the met_raw_wide_mz_single_match dataframe to a csv file in the rt_mapping_anchor_metabolites_folder_name
        write_csv(met_raw_wide_mz_single_match, paste0(rt_mapping_anchor_metabolites_folder_name, "/RT_mapping_anchor_metabolites.csv"))
    }

    # load the hmdb_metabolites_reference_retention_time.rda dataframe
    data(hmdb_metabolites_reference_retention_time)

    # keep only the first 14 characters of the InChIKey in the hmdb_metabolites_reference_retention_time
    hmdb_metabolites_reference_retention_time$INCHIKEY = substr(hmdb_metabolites_reference_retention_time$INCHIKEY, 1, 14)

    # select the appropriate liquid chromatography setting for RT mapping based on the "LC" parameter
    if (LC == "HILIC"){
        print("Using HILIC reference library for retention time mapping")
        # keep only INCHIKEY, HILIC_RT, and HILIC_LC_run_time columns
        hmdb_metabolites_reference_retention_time = hmdb_metabolites_reference_retention_time %>%
            dplyr::select(INCHIKEY, HILIC_RT, HILIC_LC_run_time, HILIC_RT_source) %>%
            dplyr::rename(REF_RT = HILIC_RT, REF_LC_run_time = HILIC_LC_run_time, REF_RT_source = HILIC_RT_source)
    } else if (LC == "RP" | LC == "C18"){
        print("Using C18 reference library for retention time mapping")
        # keep only INCHIKEY, C18_RT, and C18_LC_run_time columns
        hmdb_metabolites_reference_retention_time = hmdb_metabolites_reference_retention_time %>%
            dplyr::select(INCHIKEY, C18_RT, C18_LC_run_time, C18_RT_source) %>%
            dplyr::rename(REF_RT = C18_RT, REF_LC_run_time = C18_LC_run_time, REF_RT_source = C18_RT_source)
    }

    # arrange by REF_RT_source and remove duplicates: this is to prioritize the experimental retention time over the predicted retention time when there are multiple retention time predictions for the same InChIKey due to metabolite inchikey duplicates
    hmdb_metabolites_reference_retention_time = hmdb_metabolites_reference_retention_time %>%
        arrange(desc(REF_RT_source)) %>%
        distinct(INCHIKEY, .keep_all = TRUE)

    # create a new column, REF_RT_normalized, by dividing the REF_RT by the REF_LC_run_time
    hmdb_metabolites_reference_retention_time$REF_RT_normalized = hmdb_metabolites_reference_retention_time$REF_RT / hmdb_metabolites_reference_retention_time$REF_LC_run_time

    # create a new column, time_normalized, by dividing the time by specified LC_run_time in the met_raw_wide_mz_single_match
    # update the LC_run_time from minutes to seconds
    LC_run_time = LC_run_time * 60
    print(paste0("The LC run time is: ", LC_run_time, " seconds"))

    # now we have very confident metabolites with clean chromatography and no isomers. This is the training dataset for the retention time prediction model.
    met_raw_wide_mz_single_match$time_normalized = met_raw_wide_mz_single_match$time / LC_run_time
    # add retention time from the reference library to the current study's metabolites by left join the InChIKey==INCHIKEY column
    met_raw_wide_mz_single_match_RT_mapping = inner_join(met_raw_wide_mz_single_match, hmdb_metabolites_reference_retention_time, by = c("InChIKey"="INCHIKEY"))

    if (detail == TRUE) {
        # save the met_raw_wide_mz_single_match dataframe to a csv file in the rt_mapping_anchor_metabolites_folder_name
        write_csv(met_raw_wide_mz_single_match_RT_mapping, paste0(rt_mapping_anchor_metabolites_folder_name, "/RT_prediction_parameter_estimation.csv"))
    }

    ############## Monotonically Constrained GAM + Robust Weighting retention time mapping algorithm (adapted from PredRet) ##############
    # Fit Retention Time Prediction Model 1: High Confidence Matches (Training Data)
    met_raw_wide_mz_single_match_RT_mapping_predret = met_raw_wide_mz_single_match_RT_mapping

    # Fit the PredRet GAM model
    RT_mapping_fit_predret = fit_predret_gam(
        anchors = met_raw_wide_mz_single_match_RT_mapping_predret,
        rt_col  = "time_normalized",
        rt_ref_col = "REF_RT_normalized"
    )

    # now we have the sigma value for those very confident metabolites with clean chromatography and no isomers. This is the sigma value of the retention time prediction among mostly correct metabolite identities (sigma_correct).
    # Statistics Calculation for Retention Time Prediction Model (Training Data)
    if(!is.null(RT_mapping_fit_predret)) {

        # Get R-squared (Approximate from the initial unconstrained fit used for weights)
        rt_mapping_r_squared_predret = RT_mapping_fit_predret$r.sq

        # Calculate Sigma (Standard Deviation of Residuals)
        # In PredRet, this includes outliers but they were weighted down during fit.
        rt_mapping_sigma_predret = sd(RT_mapping_fit_predret$residuals, na.rm = TRUE)
        
        # Convert sigma to seconds
        rt_mapping_sigma_predret = rt_mapping_sigma_predret * LC_run_time

        # set a minimum sigma value of 25 seconds, in case it becomes too optimistic for the retention time prediction model (e.g. 4 seconds due to coincidences)
        if (rt_mapping_sigma_predret < 25) {
            rt_mapping_sigma_predret = 25
        }
        
        print(paste0("The sigma value for the RT mapping fit among metabolite adducts with clean chromatography (training data): ", round(rt_mapping_sigma_predret, 2), " seconds"))
        
        # Predict retention times
        met_raw_wide_mz_single_match_RT_mapping_predret$time_predicted = predict_predret_gam(
            met_raw_wide_mz_single_match_RT_mapping_predret$REF_RT_normalized, 
            RT_mapping_fit_predret
        )
        
        # Denormalize predictions
        met_raw_wide_mz_single_match_RT_mapping_predret$time_predicted = met_raw_wide_mz_single_match_RT_mapping_predret$time_predicted * LC_run_time
        
        # Calculate MAE
        mean_difference_predret = mean(abs(met_raw_wide_mz_single_match_RT_mapping_predret$time - met_raw_wide_mz_single_match_RT_mapping_predret$time_predicted), na.rm=TRUE)
        print(paste0("The mean absolute error for retention time prediction in the 'training' data is: ", round(mean_difference_predret, 1), " seconds"))
    
    } else {
        print("Failed to fit training model (insufficient points).")
    }

    ############# establish parameters for the precursor-product/transporter correlation analysis #############
    ## average and standard deviation of the correlation coefficients
    # prepare the input data for the precursor-product correlation analysis among those background data: met_raw_wide_mz_single_match_RT_mapping_predret
    met_raw_wide_mz_single_match_RT_mapping_predret_2 = met_raw_wide_mz_single_match_RT_mapping_predret %>%
        rename(identification_type = identification_type_annotated, Adduct = Adduct_annotated)

    # prepare the input data for the precursor-product correlation analysis among those training data
    MSMICA_col_names_connection_training = precursor_product_correlation_preparation(precursor_product_correlation_preparation_input_data = met_raw_wide_mz_single_match_RT_mapping_predret_2,
        rxn_connection = reaction_connection)

    ############ calculation correlation coefficients and p-values only between the proposed precursor-product relationships ############
    # apply the precursor_product_correlation function to the MSMICA_col_names_connection_training
    precursor_product_correlation_result_training_unique = precursor_product_correlation(precursor_product_correlation_input_data = MSMICA_col_names_connection_training, precursor_col = "mz_time_identified", product_col = "mz_time_identified_final",
        cor_input = MSMICA_cor_input)

    # keep only significant precursor-product correlations
    precursor_product_correlation_result_training_unique = precursor_product_correlation_result_training_unique %>%
        filter(p_value < 0.05)

    # if detail is TRUE, save the precursor_product_correlation_result_training_unique to a file
    if (detail == TRUE) {
        write_csv(precursor_product_correlation_result_training_unique, paste0(rt_mapping_anchor_metabolites_folder_name, "/precursor_product_transporter_correlation_parameter_estimation.csv"))
    }

    # Fisher Z-transformation: this stretches the ends of the correlation scale to infinity. This makes the sampling distribution of correlations Normal (Gaussian), allowing you to validly use dnorm() for your likelihoods
    pp_z_obs = atanh(abs(precursor_product_correlation_result_training_unique$correlation))

    # extract mu and sigma from the z_obs
    pp_mu = mean(pp_z_obs)
    pp_sigma = sd(pp_z_obs)

    print(paste0("The mean of the Fisher Z-transformed precursor-product correlations for the training data is: ", round(pp_mu, 2)))
    print(paste0("The standard deviation of the Fisher Z-transformed precursor-product correlations for the training data is: ", round(pp_sigma, 2)))

    # Splitting the data by Mono_mass
    data_split_cluster = split(met_raw_wide_final_monomass_simple, met_raw_wide_final_monomass_simple$Mono_mass)

    # Before the lapply function, create a new progress bar for bayesian probability calculation for clustering patterns:
    total_iterations_cluster = length(data_split_cluster)
    pb_cluster = txtProgressBar(min = 0, max = total_iterations_cluster, style = 3)

    # performing the clustering of adducts and isotopes algorithm by calculating the Bayesian probability for each metabolite using adduct correlation, isotope correlation, and adduct formation
    message("*******************************************************************************")
    message("Performing the clustering of adduct and isotope algorithm using adduct correlation, isotope correlation, and adduct formation:")
    message("All the features annotated as the same monoisotopic mass will be clustered together.")

    # Apply the bayesian_probability_calculation_cluster function to each metabolite cluster to find features protable to be the metabolites with the same monoisotopic mass
    results_cluster = lapply(seq_along(data_split_cluster), function(j) {
        res = bayesian_probability_calculation_cluster(
            data = data_split_cluster[[j]],
            MSMICA_cor_input = MSMICA_cor_input,
            adduct_corr_time_thresh = adduct_correlation_time_threshold,
            adduct_corr_r_thresh = adduct_correlation_r_threshold,
            isotopic_corr_time_thresh = isotopic_correlation_time_threshold,
            isotopic_corr_r_thresh = isotopic_correlation_r_threshold,
            primary_df_simple = met_raw_wide_final_monomass_simple,
            isotope_df_simple = met_raw_wide_final_monomass_isotope_simple
        )
        setTxtProgressBar(pb_cluster, j)
        res
    })

    # Combine all data frames in the list into one final data frame
    final_results_cluster = bind_rows(results_cluster)

    # clear memory
    rm(results_cluster)

    # convert the final_results_cluster to a tibble
    final_results_cluster = tibble(final_results_cluster)

    # rename correlation as isotopic_correlation
    final_results_cluster = final_results_cluster %>%
        rename(isotopic_correlation = correlation)
    
    # add mean_intensity to met_raw_wide_mz_match_4
    met_raw_wide_original_mean_intensity_2 = met_raw_wide_original_mean_intensity %>% 
        mutate(mz_time = paste0(round(mz, 4), "_", round(time, 0)), mean_intensity = round(mean_intensity, 0)) %>%
        select(mz_time, mean_intensity)

    # create mz_time_isotope_annotated column by combining mz_isotope and time_isotope
    final_results_cluster$mz_time_isotope_annotated = paste0(round(final_results_cluster$mz_isotope, 4), "_", round(final_results_cluster$time_isotope, 0))
    # relocate the mz_time_isotope_annotated column right after the time_isotope column
    final_results_cluster = final_results_cluster %>%
        relocate(mz_time_isotope_annotated, .after = time_isotope)

    # add the mean_intensity to the final_results_cluster by left joining with the met_raw_wide_original_mean_intensity_2 based on mz_time_annotated == mz_time (mean intensity of the features without considering isotopologues)
    final_results_cluster = final_results_cluster %>%
        left_join(met_raw_wide_original_mean_intensity_2, by = c("mz_time_annotated" = "mz_time"))

    # add the mean_intensity to the final_results_cluster by left joining with the met_raw_wide_original_mean_intensity_2 based on mz_time_isotope_annotated == mz_time (mean intensity of the features considering isotopologues)
    final_results_cluster = final_results_cluster %>%
        left_join(met_raw_wide_original_mean_intensity_2, by = c("mz_time_isotope_annotated" = "mz_time"))

    # rename the mean_intensity.x column to mean_intensity and mean_intensity.y column to mean_intensity_isotope
    final_results_cluster = final_results_cluster %>%
        rename(mean_intensity = mean_intensity.x, mean_intensity_isotope = mean_intensity.y) %>%
        ## relocate the mean_intensity and mean_intensity_isotope columns right after the time_isotope column
        relocate(mean_intensity, .after = mz_time_annotated) %>%
        relocate(mean_intensity_isotope, .after = mz_time_isotope_annotated)

    # if detail == TRUE, save the final_results_cluster in the clustering_folder_name folder
    if (detail == TRUE){
        write_csv(final_results_cluster, paste0(clustering_folder_name, "/", "clustering_of_adduct_isotope_all_result.csv"))
    }

    # extract the adduct correlation among the features with the same monoisotopic mass and close retention time
    adduct_correlation_long_df = extract_adduct_correlation_from_cluster_results(
        final_results_cluster = final_results_cluster,
        MSMICA_cor_input = MSMICA_cor_input,
        time_threshold = adduct_correlation_time_threshold
    )

    # if detail == TRUE, save the adduct_correlation_long_df to a file
    if (detail == TRUE) {
        write_csv(adduct_correlation_long_df, paste0(clustering_folder_name, "/adduct_correlation_result.csv"))
    }

    # select only protonated, deprotonated, sodium, potassium, and chloride adducts
    ## because protonated and deprotonated adducts are the most common adducts for many metabolites and sodium, potassium and chloride adducts are common for carbohydrate metabolites, we prioritize them to mitigate the mz confidence issue, which means the secondary adducts of abundant metabolites may be wrongly assigned to the protonated or deprotonated adduct of less abundant metabolites with the same m/z.
    ## the adducts include positive and negative adducts, which should not matter here, because the adduct choice are preselected based on the All_Adduct parameter in the MSMICA_algorithm function,  as the oppositely charged adducts will not be considered at the front end.
    final_results_cluster_identified = final_results_cluster %>%
        dplyr::filter(
            (Adduct_annotated %in% c("M-H","M+H") & MSMICA_identification == 1) | Adduct_annotated %in% c("M+Na","M+K","M+Cl") | Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") | Adduct_annotated %in% c("M+2Na-H","M+2K-H")
            )

    # select only the Mono_mass, Name, KEGG_ID, HMDB_ID, InChIKey, and Subclass columns
    metabolite_database_simple = metabolite_database %>%
        dplyr::select(Mono_mass, Name, KEGG_ID, HMDB_ID, InChIKey, Subclass, has_NH3, has_H2O, n_NH3, n_H2O) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1)) %>%
        ## round the Mono_mass to 4 decimals
        mutate(Mono_mass = round(Mono_mass, 4))
    
    # select only InChIKey and Subclass columns
    metabolite_database_simple_subclass = metabolite_database_simple %>%
        dplyr::select(InChIKey, Subclass) %>%
        distinct()

    # inner join the final_results_cluster_identified with the metabolite_database_simple to get the final result using Mono_mass
    final_results_cluster_identified_2 = final_results_cluster_identified %>% 
        inner_join(metabolite_database_simple, by = "Mono_mass", relationship = "many-to-many") %>%
        ## left join with the hmdb_metabolites_concentrations_average_simple based on InChIKey
        left_join(hmdb_metabolites_concentrations_average_simple, by = "InChIKey", relationship = "many-to-many")

    # remove the following adducts if the following conditions are met. This helps to reduce the number of adducts to be considered for the m/z matching given the metabolite structure.
    final_results_cluster_identified_2 = final_results_cluster_identified_2 %>%
        # Remove rows based on adduct-specific conditions and include the adducts that are not affected by the conditions
        filter(
            (Adduct_annotated == "M+H-NH3" & has_NH3 == TRUE) |
            (Adduct_annotated %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O") & n_H2O >= 1) |
            (Adduct_annotated %in% c("M-2H2O+H", "M+H-2H2O") & n_H2O >= 2) |
            (!Adduct_annotated %in% c("M+H-NH3", "M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O", "M-2H2O+H", "M+H-2H2O"))
        ) %>%
        dplyr::select(-c(has_NH3, has_H2O, n_NH3, n_H2O))

    # select only protonated and deprotonated adducts unless the subclass is "Carbohydrates and carbohydrate conjugates" and "Bile acids, alcohols and derivatives" and "Retinoids"
    ## this is because the carbohydrates usually form M+Na (in plasma), M+K (in tissue), and M+Cl (in plasma) adducts, while the other subclasses usually form M-H and M+H adducts
    ## bile acids, alcohols and derivatives and retinoids are likely to form M+H-H2O adducts
    ## Pregnane steroids are likely to form M+H-2H2O and M+H-H2O adducts
    ## gluconic-acid-like compounds are especially prone to strong sodium adduction, and M+2Na−H is chemically very plausible for it in positive-mode MS. Classic chemistry on metal–gluconate complexes shows that gluconate is a strong metal-binding ligand
    final_results_cluster_identified_2 = final_results_cluster_identified_2 %>%
        dplyr::filter(
            (Adduct_annotated %in% c("M-H","M+H")) | 
            (Adduct_annotated %in% c("M+Na","M+K","M+Cl") & Subclass == "Carbohydrates and carbohydrate conjugates") | 
            (Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & Subclass == "Bile acids, alcohols and derivatives") | 
            (Adduct_annotated %in% c("M+H-H2O") & (Subclass == "Retinoids")) | 
            (Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & Subclass == "Pregnane steroids") |
            (Adduct_annotated %in% c("M+Na", "M+K", "M+2Na-H", "M+2K-H") & InChIKey %in% c("RGHNJXZEOKUKBD", "BIRSGZKFKXLSJQ", "NBFWIISVIFCMDK", "PALQXFMLVVWXFD", "IWHLYPDWHHPVAA", "SBCIXDBITAKZCS", "YGMNHEPVTNXLLS", "KWMLJOLKUYYJFJ"))
            )

    # rename "identification_type_annotated" as "identification_type", "Adduct_annotated" as "Adduct", "mz_annotated" as "mz", "time_annotated" as "time"
    final_results_cluster_identified_2 = final_results_cluster_identified_2 %>%
        rename(identification_type = identification_type_annotated, Adduct = Adduct_annotated, mz = mz_annotated, time = time_annotated, identified_Name = Name) %>%
        # create the identification_method column with value "mz matching"
        mutate(identification_method = paste0("mz matching", "; clustering of adducts and isotopes"))

    # select specific columns from final_results
    final_results_cluster_identified_3 = final_results_cluster_identified_2 %>%
        mutate(match_category = "multiple") %>%
        dplyr::select(Mono_mass, identification_type, ion_mode, identified_Name, KEGG_ID, HMDB_ID, InChIKey, Adduct, mz, time, MSMICA_identification, identification_method, Concentration_average, Probability, match_category) %>%
        filter(!is.na(InChIKey))

    # round mz to 4 decimal places and time to 0 decimal places
    final_results_cluster_identified_3$mz = round(final_results_cluster_identified_3$mz, 4)
    final_results_cluster_identified_3$time = round(final_results_cluster_identified_3$time, 0)
    # create mz_time column in the final_results_cluster_identified_3
    final_results_cluster_identified_3$mz_time = paste0(final_results_cluster_identified_3$mz, "_", final_results_cluster_identified_3$time)
    # relocate the mz_time column right after the time column
    final_results_cluster_identified_3 = final_results_cluster_identified_3 %>%
        dplyr::relocate(mz_time, .after = time)

    # if detail == TRUE, save the final_results_cluster_identified_3 in the clustering_folder_name folder
    if (detail == TRUE){
        write_csv(final_results_cluster_identified_3, paste0(clustering_folder_name, "/", "clustering_of_adduct_isotope_main_result.csv"))
    }

    # extract the clustered features so that they can be excluded from later analyses to reduce the risk of m/z coincidence ########################################
    print("features selected by the clustering algorithm and single-match m/z and concentration are excluded for the later analysis")

    # select only the InChIKey and Name columns
    metabolite_database_simple = metabolite_database %>%
        dplyr::select(InChIKey, HMDB_ID, KEGG_ID, Name, Formula, Mono_mass, has_NH3, has_H2O, n_NH3, n_H2O) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1)) %>%
        # round the Mono_mass to 4 decimal places
        mutate(Mono_mass = round(Mono_mass, 4))

    # find the Mono_mass of the mz_time_clustered_feature in the final_results_cluster
    ## here, the mz and time are those metabolites selected by the clustering algorithm (typically M+H or M-H adducts)
    Mono_mass_clustered_feature = final_results_cluster_identified_3 %>%
        select(InChIKey, Adduct, mz, time) %>%
        left_join(metabolite_database_simple, by = c("InChIKey" = "InChIKey"))

    # remove the following adducts if the following conditions are met. This helps to reduce the number of adducts to be considered for the m/z matching given the metabolite structure.
    Mono_mass_clustered_feature = Mono_mass_clustered_feature %>%
        filter(
            (Adduct == "M+H-NH3" & has_NH3 == TRUE) |
            (Adduct %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O") & n_H2O >= 1) |
            (Adduct %in% c("M-2H2O+H", "M+H-2H2O") & n_H2O >= 2) |
            (!Adduct %in% c("M+H-NH3", "M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O", "M-2H2O+H", "M+H-2H2O"))
        ) %>%
        dplyr::select(-c(has_NH3, has_H2O, n_NH3, n_H2O))

    # create mz_time_adduct column by combining mz and Adduct
    Mono_mass_clustered_feature_mz_time_adduct = paste0(Mono_mass_clustered_feature$mz, "_", Mono_mass_clustered_feature$time, "_", Mono_mass_clustered_feature$Adduct)

    # find those features with secondary adducts that are within the time threshold of the primary adduct
    final_results_cluster_multiple_match = Mono_mass_clustered_feature  %>%
        inner_join(final_results_cluster, by = c("Mono_mass" = "Mono_mass"), relationship = "many-to-many") %>%
        # remove those rows with abs(time-time_annotated) > adduct_correlation_time_threshold or abs(time_isotope-time_annotated) > isotopic_correlation_time_threshold
        filter(abs(time-time_annotated) <= adduct_correlation_time_threshold | abs(time_isotope-time) <= isotopic_correlation_time_threshold) %>%
        # create mz_time_adduct_annotated column by combining mz_annotated and time_annotated and Adduct_annotated
        mutate(mz_time_adduct_annotated = paste0(mz_annotated, "_", time_annotated, "_", Adduct_annotated)) %>%
        # remove where mz_time_adduct_annotated is in Mono_mass_clustered_feature_mz_time_adduct
        filter(!(mz_time_adduct_annotated %in% Mono_mass_clustered_feature_mz_time_adduct)) %>%
        # remove the mz_time_adduct_annotated column
        dplyr::select(-mz_time_adduct_annotated) %>%
        # remove Adduct, mz, time
        dplyr::select(-Adduct, -mz, -time) %>%
        # remove these adducts,  M+H-2H2O  M+H-H2O  M+H-NH3, as they may be the protonated forms of the metabolites (e.g. 3-Hydroxydecanoyl carnitine (M-H2O+H) and Acylcarnitine C10:1 (M+H) have the same m/z)
        filter(!(Adduct_annotated %in% c("M+H-2H2O", "M+H-H2O", "M+H-NH3")))

    # round the mz_annotated to 4 decimal places and time_annotated to 0 decimal places
    final_results_cluster_multiple_match$mz_annotated = round(final_results_cluster_multiple_match$mz_annotated, 4)
    final_results_cluster_multiple_match$time_annotated = round(final_results_cluster_multiple_match$time_annotated, 0)
    final_results_cluster_multiple_match$mz_isotope = round(final_results_cluster_multiple_match$mz_isotope, 4)
    final_results_cluster_multiple_match$time_isotope = round(final_results_cluster_multiple_match$time_isotope, 0)

    # get the mz_time_adduct data by combining mz_annotated and time_annotated
    mz_time_adduct = paste0(final_results_cluster_multiple_match$mz_annotated, "_", final_results_cluster_multiple_match$time_annotated)
    # get the mz_time_isotope data by combining mz_isotope and time_isotope
    mz_time_isotope = paste0(final_results_cluster_multiple_match$mz_isotope, "_", final_results_cluster_multiple_match$time_isotope)
    # combine the mz_time_adduct and mz_time_isotope to get mz_time_adduct_isotope
    mz_time_adduct_isotope = c(mz_time_adduct, mz_time_isotope)
    # remove "NA_NA"
    mz_time_adduct_isotope = mz_time_adduct_isotope[mz_time_adduct_isotope != "NA_NA"]

    # combine the mz_time_clustered_feature and mz_time_adduct_isotope
    mz_time_clustered_feature = mz_time_adduct_isotope
    # remove duplicates
    mz_time_clustered_feature = unique(mz_time_clustered_feature)

    ##### use m/z matching, retention time prediction, and biospecimen-specific concentration to annotate metabolites for metabolomics features ############################# 
    # exclude clustered features by the clustering of adduct and isotope algorithm
    met_raw_wide_mz_match = met_raw_wide_final_monomass %>%
        ## left join with the hmdb_metabolites_concentrations_average_simple based on InChIKey
        left_join(hmdb_metabolites_concentrations_average_simple, by = "InChIKey", relationship = "many-to-many") %>%
        ## inner join with the metabolite_database_simple_subclass based on InChIKey
        inner_join(metabolite_database_simple_subclass, by = "InChIKey", relationship = "many-to-many") %>%
        dplyr::filter(
            Adduct_annotated %in% c("M-H","M+H")| 
            (Adduct_annotated %in% c("M+Na","M+K","M+Cl") & Subclass == "Carbohydrates and carbohydrate conjugates") | 
            (Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & Subclass == "Bile acids, alcohols and derivatives") | 
            (Adduct_annotated %in% c("M+H-H2O") & (Subclass == "Retinoids")) | 
            (Adduct_annotated %in% c("M+H-H2O", "M+H-2H2O") & Subclass == "Pregnane steroids") |
            (Adduct_annotated %in% c("M+Na", "M+K", "M+2Na-H", "M+2K-H") & InChIKey %in% c("RGHNJXZEOKUKBD", "BIRSGZKFKXLSJQ", "NBFWIISVIFCMDK", "PALQXFMLVVWXFD", "IWHLYPDWHHPVAA", "SBCIXDBITAKZCS", "YGMNHEPVTNXLLS", "KWMLJOLKUYYJFJ"))
            ) %>%
        dplyr::select(-Subclass) %>%
        # check if the feature (mz_time_annotated) has any quantified metabolites annotated
        group_by(mz_time_annotated) %>%
        mutate(has_any_concentration = any(!is.na(Concentration_average))) %>%
        ungroup() %>%
        ## for those Concentration_average as NA, exclude mz_time_sample that are in the mz_time_clustered_feature
        filter(!(mz_time_annotated %in% mz_time_clustered_feature & !has_any_concentration)) %>%
        dplyr::select(-has_any_concentration)


    # add KEGG_ID and HMDB_ID to the met_raw_wide_mz_match by left joining with the metabolite_database_simple based on Mono_mass, Confirmed_Name=Name, and InChIKey
    met_raw_wide_mz_match = left_join(met_raw_wide_mz_match, metabolite_database_simple, by = c("Mono_mass", "Confirmed_Name"="Name", "InChIKey"), relationship = "many-to-many")
    
    # arrange by mz_annotated and time_annotated
    met_raw_wide_mz_match = met_raw_wide_mz_match %>%
        arrange(mz_annotated, time_annotated)
    
    # Add identification_method as "mz matching" to the met_raw_wide_mz_match
    met_raw_wide_mz_match$identification_method = "mz matching"

    # within the met_raw_wide_mz_match, group by Mono_mass, Confirmed_Name, Adduct_annotated, and then create match_category = "single" if there is only 1 row for each Mono_mass, Confirmed_Name, Adduct_annotated
    met_raw_wide_mz_match = met_raw_wide_mz_match %>%
        group_by(Mono_mass, Confirmed_Name, Adduct_annotated) %>%
        mutate(match_category = if_else(n() == 1, "single", "multiple")) %>%
        ungroup() %>%
        rename(mz = mz_annotated, time = time_annotated, mz_time = mz_time_annotated, identified_Name = Confirmed_Name)

    # rename met_raw_wide_mz_match
    met_raw_wide_mz_match = met_raw_wide_mz_match %>%
        rename(identification_type = identification_type_annotated, Adduct = Adduct_annotated)
    
    # remove columns from final_results_cluster_identified_3: MSMICA_identification, Probability
    final_results_cluster_identified_4 = final_results_cluster_identified_3 %>% dplyr::select(-MSMICA_identification, -Probability)

    # order the columns of met_raw_wide_mz_match like final_results_cluster_identified_4: Mono_mass, identification_type, ion_mode, identified_Name, KEGG_ID, HMDB_ID, InChIKey, Adduct, mz, time, mz_time, identification_method, Concentration_average, match_category
    met_raw_wide_mz_match_2 = met_raw_wide_mz_match %>% select(Mono_mass, identification_type, ion_mode, identified_Name, KEGG_ID, HMDB_ID, InChIKey, Adduct, mz, time, mz_time, identification_method, Concentration_average, match_category)

    # combine final_results_cluster_identified_4 and met_raw_wide_mz_match_2
    met_raw_wide_mz_match_3 = bind_rows(final_results_cluster_identified_4, met_raw_wide_mz_match_2) %>%
        arrange(mz, time) %>% 
        distinct(Mono_mass, KEGG_ID, HMDB_ID, InChIKey, Adduct, mz, time, .keep_all = TRUE)

    # add predicted retention time to met_raw_wide_mz_match from the retention time reference library based on InChIKey
    met_raw_wide_mz_match_3 = inner_join(met_raw_wide_mz_match_3, hmdb_metabolites_reference_retention_time, by = c("InChIKey"="INCHIKEY"))

    # predict the retention time using the predict_predret_gam function and RT_mapping_fit_predret
    met_raw_wide_mz_match_3$time_predicted = predict_predret_gam(met_raw_wide_mz_match_3$REF_RT_normalized, RT_mapping_fit_predret)
    met_raw_wide_mz_match_3$time_predicted = met_raw_wide_mz_match_3$time_predicted * LC_run_time

    # calculate the absolute difference between time and time_predicted
    met_raw_wide_mz_match_3$time_difference = abs(met_raw_wide_mz_match_3$time - met_raw_wide_mz_match_3$time_predicted)

    # among all results with the same Mono_mass, if all results have time_difference >= 200 seconds, keep them all; Otherwise, exclude those results with time_difference >= 200 seconds and keep the rest
    met_raw_wide_mz_match_4 = met_raw_wide_mz_match_3 %>%
        group_by(Mono_mass) %>%
        mutate(all_time_difference_200 = all(time_difference >= 200)) %>%
        ungroup() %>%
        filter(all_time_difference_200 | time_difference < 200) %>%
        dplyr::select(-all_time_difference_200)

    # if Concentration_average is NA, set it as 1e-6
    met_raw_wide_mz_match_4$Concentration_average[is.na(met_raw_wide_mz_match_4$Concentration_average)] = 1e-6

    # left join met_raw_wide_mz_match_4 with met_raw_wide_original_mean_intensity_2 based on mz_time
    met_raw_wide_mz_match_5 = left_join(met_raw_wide_mz_match_4, met_raw_wide_original_mean_intensity_2, by = "mz_time")

    # create mz_time column by pasting the mz and time columns together
    met_raw_wide_mz_match_5$mz_time = paste0(met_raw_wide_mz_match_5$mz, "_", met_raw_wide_mz_match_5$time)

    # apply the precursor_product_correlation_processing function to the identified_metabolite_final
    ## precursor-product correlations
    cor_MSMICA_precursor_product = precursor_product_correlation_processing(met_raw_wide_mz_match_5,
        rxn_connection = reaction_connection,
        cor_input = MSMICA_cor_input)

    # apply the backpropagation_correlation_analysis function to the MSMICA_col_names_connection_identified_final and cor_MSMICA_precursor_product
    ## this is for regular backpropagation correlation analysis among precursors and products or metabolites sharing the same transporters. Therefore, the backpropagation_correlation_direction is set to "positive" and FDR_correction is set to FALSE.
    MSMICA_col_names_connection_identified_final_precursor_product = backpropagation_correlation_analysis(MSMICA_col_names_connection_12 = cor_MSMICA_precursor_product$MSMICA_col_names_connection_12, cor_MSMICA = cor_MSMICA_precursor_product$cor_MSMICA, backpropagation_correlation_direction = backpropagation_correlation_direction, FDR_correction = FALSE, duplicate_removal = TRUE)

    # select only the $MSMICA_col_names_connection_identified_final_4 within the MSMICA_col_names_connection_identified_final_precursor_product
    MSMICA_col_names_connection_identified_final_precursor_product = MSMICA_col_names_connection_identified_final_precursor_product$MSMICA_col_names_connection_identified_final_4
    # arrange by mz and time
    MSMICA_col_names_connection_identified_final_precursor_product = MSMICA_col_names_connection_identified_final_precursor_product %>%
        arrange(mz, time)

    # extract the identified_Name, mz_time, Adduct, and correlation columns from MSMICA_col_names_connection_identified_final_precursor_product
    precursor_product_transporter_correlation_extracted = MSMICA_col_names_connection_identified_final_precursor_product %>%
        select(identified_Name, mz_time, Adduct, correlation) %>%
        ## exclude correlation that is NA
        filter(!is.na(correlation))

    # left join met_raw_wide_mz_match_5 and precursor_product_transporter_correlation_extracted based on identified_Name, mz_time, and Adduct
    met_raw_wide_mz_match_6 = left_join(met_raw_wide_mz_match_5, precursor_product_transporter_correlation_extracted, by = c("identified_Name", "mz_time", "Adduct"))

    # group by mz_time and priortize those metabolites that are detected and quantified in HMDB
    met_raw_wide_mz_match_6_clean = met_raw_wide_mz_match_6 %>%
        group_by(mz_time) %>%
        # 1. Group by each unique feature, and then mark rows that are in the detected+quantified HMDB/InChIKey sets
        mutate(
            in_detected = HMDB_ID %in% hmdb_endogenous_metabolites_detected_and_quantified_HMDBID |
                        InChIKey %in% hmdb_endogenous_metabolites_detected_and_quantified_InChIKey,
            any_in_detected = any(in_detected)
        ) %>%
        # 2. If any_in_detected is TRUE for this mz_time, keep only in_detected rows;
        #    otherwise keep everything.
        filter(if_else(any_in_detected, in_detected, TRUE)) %>%
        # 3. Clean up helper columns
        select(-in_detected, -any_in_detected) %>%
        ungroup()

    # for all results in each mz_time, prioritize the metabolites that quantified (Concentration_average > 1e-6)
    met_raw_wide_mz_match_6_clean = met_raw_wide_mz_match_6_clean %>%
        group_by(mz_time) %>%
        # 1. Group by each unique feature, and then mark rows that are in the detected+quantified HMDB/InChIKey sets
        mutate(
            quantified = Concentration_average > 1e-6,
            any_quantified = any(quantified)
        ) %>%
        # 2. If any_quantified is TRUE for this mz_time, keep only quantified rows;
        #    otherwise keep everything.
        filter(if_else(any_quantified, quantified, TRUE)) %>%
        # 3. Clean up helper columns
        select(-quantified, -any_quantified) %>%
        ungroup()

    # calculate the log_concentration_prior and log_mean_intensity
    met_raw_wide_mz_match_6_clean = met_raw_wide_mz_match_6_clean %>%
        mutate(
            log_concentration_prior = log(case_when(
                !is.na(Concentration_average) & Concentration_average > 1e-6 ~ Concentration_average,
                TRUE ~ 1e-6
            )),
            log_mean_intensity = log(case_when(
                !is.na(mean_intensity) & mean_intensity > 0 ~ mean_intensity,
                TRUE ~ 1
            ))
        )

    # group by Mono_mass and if there are M+H or M-H, keep only these rows. If not, keep all the rows.
    met_raw_wide_mz_match_6_clean = met_raw_wide_mz_match_6_clean %>%
        group_by(Mono_mass) %>%
        filter(
            if (any(Adduct %in% c("M+H", "M-H"), na.rm = TRUE)) {
                Adduct %in% c("M+H", "M-H")
            } else {
                TRUE
            }
        ) %>%
        ungroup()

    # performing the clustering of adducts and isotopes algorithm by calculating the Bayesian probability for each metabolite using adduct correlation, isotope correlation, and adduct formation
    message("*******************************************************************************")
    message("Performing the local optimization per feature using common metabolomics adducts:")
    # Features annotated as metabolites with the same monoisotopic mass will be evaluated using the Bayesian statistics with the evidence of retention time prediction, precursor-product/transporter correlation, and biospecimen-specific concentration.
    # Most likely metabolite identities will be selected for each metabolomics feature

    # iterative optimization for each feature
    local_opt_output = run_local_optimization_iteratively(
        input_data = met_raw_wide_mz_match_6_clean,
        rt_sigma = rt_mapping_sigma_predret,
        pp_mu = pp_mu,
        pp_sigma = pp_sigma,
        rxn_connection = reaction_connection,
        cor_input = MSMICA_cor_input,
        adduct_corr_time_thresh = adduct_correlation_time_threshold,
        detail = detail,
        output_folder = local_optimization_per_feature_folder_name,
        max_iter = 99
    )

    # extract the result from the local_opt_output
    MSMICA_local_optimization_per_feature_result_combined = local_opt_output$result

    # remove correlation, log_concentration_prior, and log_mean_intensity columns
    MSMICA_local_optimization_per_feature_result_combined = MSMICA_local_optimization_per_feature_result_combined %>%
        dplyr::select(-correlation, -log_concentration_prior, -log_mean_intensity)

    # save the MSMICA_local_optimization_per_feature_result_combined in the local_optimization_per_feature_folder_name folder
    if (detail == TRUE){
        write_csv(MSMICA_local_optimization_per_feature_result_combined, paste0(local_optimization_per_feature_folder_name, "/", "MSMICA_local_optimization_per_feature_main_adduct_result_combined.csv"))
    }

    # select the following columns to update the identified_metabolites
    MSMICA_results_filtered_final = MSMICA_local_optimization_per_feature_result_combined %>%
        dplyr::select(identification_type, ion_mode, identified_Name, KEGG_ID, HMDB_ID, InChIKey, Adduct, mz, time, time_predicted, time_difference, identification_method, Probability, match_category)

    # round Probability column to 1 decimal places
    MSMICA_results_filtered_final$Probability = round(MSMICA_results_filtered_final$Probability, 1)

    # if Probability is > 95, set it to 95: this is because we want to be conservative in the identification of metabolites in case it may be a unknown metabolite that is not discovered yet.
    MSMICA_results_filtered_final$Probability = ifelse(MSMICA_results_filtered_final$Probability > 95, 95, MSMICA_results_filtered_final$Probability)

    # if mz_time is not available, create it by pasting the mz and time columns together
    if (!("mz_time" %in% colnames(met_raw_wide_original_mean_intensity))) {
        met_raw_wide_original_mean_intensity$mz = round(met_raw_wide_original_mean_intensity$mz, 4)
        met_raw_wide_original_mean_intensity$time = round(met_raw_wide_original_mean_intensity$time, 0)
        met_raw_wide_original_mean_intensity$mz_time = paste0(met_raw_wide_original_mean_intensity$mz, "_", met_raw_wide_original_mean_intensity$time)
        # remove the mz and time columns
        met_raw_wide_original_mean_intensity = met_raw_wide_original_mean_intensity %>%
            dplyr::select(-mz, -time)
    }

    # create the mz_time column by pasting the mz and time columns together
    MSMICA_results_filtered_final$mz_time = paste0(MSMICA_results_filtered_final$mz, "_", MSMICA_results_filtered_final$time)

    # add mean intensity values to all confirmed and identified metabolite by joining the met_raw_wide_mean_intensity based on mz and time
    MSMICA_results_filtered_final_2 = MSMICA_results_filtered_final %>%
        left_join(met_raw_wide_original_mean_intensity, by = "mz_time") %>%
        ## round mean_intensity to integers
        mutate(mean_intensity = round(mean_intensity, 0)) %>%
        ## move mean_intensity column right after mz_time
        relocate(mean_intensity, .after = "mz_time") %>%
        ## remove identified_Name column
        dplyr::select(-identified_Name)

    # select only necessary columns from metabolite_database
    metabolite_database_simple = metabolite_database %>%
        dplyr::select(PubChem_compound_id, LIPID_MAPS_ID, ChEBI_ID, BioCyc_ID, DrugBank_ID, Name, Mono_mass, Formula, Charge, SMILES, InChIKey, Kingdom, Superclass, Class, Subclass, Direct_parent) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1)) %>%
        ## round the Mono_mass to 4 decimal places
        mutate(Mono_mass = round(Mono_mass, 4))
    
    # left join the metabolite_database_simple to the MSMICA_results_filtered_final_2 using InChIKey
    MSMICA_results_filtered_final_3 = left_join(MSMICA_results_filtered_final_2, metabolite_database_simple, by = c("InChIKey" = "InChIKey"), relationship = "many-to-many") %>%
        # prioritize the metabolite_database_simple columns
        dplyr::select(PubChem_compound_id, HMDB_ID, KEGG_ID, LIPID_MAPS_ID, ChEBI_ID, BioCyc_ID, DrugBank_ID, Name, Mono_mass, Formula, Charge, SMILES, InChIKey, Kingdom, Superclass, Class, Subclass, Direct_parent, everything())
    
    # rename columns: feature_metabolite_category to metabolite_within_feature_number
    MSMICA_results_filtered_final_3 = MSMICA_results_filtered_final_3 %>%
        ## remove the match_category column
        dplyr::select(-match_category) %>%
        group_by(mz_time) %>%
        ## summarize the number of metabolites within each feature
        mutate(metabolite_within_feature_number = n()) %>%
        ungroup()
    
    # relocate Probability right after metabolite_within_feature_number
    MSMICA_results_filtered_final_3 = MSMICA_results_filtered_final_3 %>%
        dplyr::relocate(Probability, .after = "metabolite_within_feature_number")
    
    # round the mz to 4 decimal places and time to integers in met_raw_wide_original
    met_raw_wide_original$mz = round(met_raw_wide_original$mz, 4)
    met_raw_wide_original$time = round(met_raw_wide_original$time, 0)

    # left join the MSMICA_results_filtered_final_3 with the met_raw_wide_original to get all intensity values
    MSMICA_results_filtered_final_4 = MSMICA_results_filtered_final_3 %>%
        left_join(met_raw_wide_original, by = c("mz" = "mz", "time" = "time"))
    
    # arrange the MSMICA_results_filtered_final_4 first by HMDB_ID, then by KEGG_ID
    MSMICA_results_filtered_final_4 = MSMICA_results_filtered_final_4 %>%
        arrange(KEGG_ID, HMDB_ID)

    # remove duplicates
    MSMICA_results_filtered_final_4 = MSMICA_results_filtered_final_4 %>%
        distinct()

    # round time_predicted and time_difference to 1 decimal place
    MSMICA_results_filtered_final_4$time_predicted = round(MSMICA_results_filtered_final_4$time_predicted, 1)
    MSMICA_results_filtered_final_4$time_difference = round(MSMICA_results_filtered_final_4$time_difference, 1)

    # reorder the identification_method column based on the canonical order
    canonical_order = c(
        "mz matching",
        "retention time prediction",
        "biospecimen-specific concentration",
        "precursor-product/transporter correlation",
        "clustering of adducts and isotopes"
    )

    reorder_methods = function(x, sep = ";\\s*") {
        parts = unlist(strsplit(x, sep))
        paste(canonical_order[canonical_order %in% parts], collapse = "; ")
    }

    MSMICA_results_filtered_final_4$identification_method = vapply(MSMICA_results_filtered_final_4$identification_method, reorder_methods, character(1))

    # update the identification_type column to indicate which Schymanski level should be used for each annotation based on identification_method column
    ## if identification_method equal to "mz matching; retention time prediction" or "mz matching; retention time prediction; clustering of adducts and isotopes", then it should be "Schymanski level 3b" because it only has mass spectromery evidence (m/z matching) and chromatographic evidence from RT prediction and/or clustering of adducts and isotopes
    ## if identification_method equal to other cases where more pieces of evidence are used, then it should be "Schymanski level 3a" because it has mass spectromery evidence, chromatographic evidence, and biological evidence (biospecimen-specific concentration and/or precursor-product/transporter correlation)
    MSMICA_results_filtered_final_5 = MSMICA_results_filtered_final_4 %>%
        mutate(
            identification_type = case_when(
                identification_method == "mz matching; retention time prediction" | identification_method == "mz matching; retention time prediction; clustering of adducts and isotopes" ~ "Schymanski level 3b",
                TRUE ~ "Schymanski level 3a"
            )
        ) %>%
        # relocate identification_type right after identification_method
        relocate(identification_type, .after = "identification_method")

    # see which metabolites have secondary adducts and/or isotopic adducts in the final_results_cluster with the same Mono_mass ##########################
    metabolite_database_simple = metabolite_database %>%
        dplyr::select(Mono_mass, Name, InChIKey, has_NH3, has_H2O, n_NH3, n_H2O) %>%
        ## select only the values in the Name column before the first ;
        mutate(Name = strsplit(Name, ";") %>% sapply(`[[`, 1)) %>%
        ## round the Mono_mass to 4 decimals
        mutate(Mono_mass = round(Mono_mass, 4))

    # inner join the final_results_cluster with the metabolite_database_simple to get the final result using Mono_mass
    final_results_cluster_simple = inner_join(final_results_cluster, metabolite_database_simple, by = "Mono_mass", relationship = "many-to-many")

    # remove the following adducts if the following conditions are met. This helps to reduce the number of adducts to be considered for the m/z matching given the metabolite structure.
    final_results_cluster_simple = final_results_cluster_simple %>%
        filter(
            (Adduct_annotated == "M+H-NH3" & has_NH3 == TRUE) |
            (Adduct_annotated %in% c("M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O") & n_H2O >= 1) |
            (Adduct_annotated %in% c("M-2H2O+H", "M+H-2H2O") & n_H2O >= 2) |
            (!Adduct_annotated %in% c("M+H-NH3", "M-H2O+H", "M+H-H2O", "M-H2O-H", "M-H-H2O", "M-2H2O+H", "M+H-2H2O"))
        ) %>%
        dplyr::select(-c(has_NH3, has_H2O, n_NH3, n_H2O))

    # select only the necessary columns related to isotope
    final_results_cluster_simple_isotope = final_results_cluster_simple %>%
        dplyr::select(InChIKey, Adduct_annotated, mz_isotope, time_isotope) %>%
        ## keep only where mz_isotope is not NA
        filter(!is.na(mz_isotope))

    # select only the necessary columns
    final_results_cluster_simple_2 = final_results_cluster_simple %>%
        dplyr::select(InChIKey, Adduct_annotated, mz_annotated, time_annotated, adduct_corr_cluster) %>%
        ## keep only the rows where adduct_corr_cluster is TRUE
        filter(adduct_corr_cluster == TRUE)

    # remove the intensity values to the MSMICA_results_filtered_final_4 (those columns starting at 31th column)
    MSMICA_results_filtered_final_4_simple = MSMICA_results_filtered_final_4[, 1:31]
    
    # see which metabolites have secondary adducts and/or isotopic adducts in the final_results_cluster with the same InChIKey
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = final_results_cluster_simple_2 %>% 
        inner_join(MSMICA_results_filtered_final_4_simple, by = c("InChIKey" = "InChIKey"))

    # left join the MSMICA_results_filtered_final_4_secondary_adduct_isotopic with the final_results_cluster_simple_isotope by InChIKey and Adduct_annotated
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        left_join(final_results_cluster_simple_isotope, by = c("InChIKey" = "InChIKey", "Adduct_annotated" = "Adduct_annotated"))

    # create a new column called adduct_type, which is "primary" for all rows in MSMICA_results_filtered_final_4_secondary_adduct_isotopic where Adduct_annotated == Adduct, mz_annotated == mz, time_annotated == time
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        mutate(adduct_type = ifelse(Adduct_annotated == Adduct & mz_annotated == mz & time_annotated == time, "primary", "secondary"))

    # keep only those rows with abs(time-time_annotated) <= adduct_correlation_time_threshold
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        filter(abs(time-time_annotated) <= adduct_correlation_time_threshold)

    # remove the columns: Adduct, mz, time, Then, rename Adduct_annotated as Adduct, mz_annotated as mz, and time_annotated as time
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        dplyr::select(-Adduct, -mz, -time) %>%
        rename(Adduct = Adduct_annotated, mz = mz_annotated, time = time_annotated)

    # move Adduct, mz, time, mz_isotope, and time_isotope columns right after identification_type column
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        relocate(Adduct, .after = "identification_type") %>%
        relocate(mz, .after = "Adduct") %>%
        relocate(time, .after = "mz") %>%
        relocate(mz_isotope, .after = "time") %>%
        relocate(time_isotope, .after = "mz_isotope")
    
    # round the mz to 4 decimal places and time to intergers
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic$mz = round(MSMICA_results_filtered_final_4_secondary_adduct_isotopic$mz, 4)
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time = round(MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time, 0)

    # update the mz_time column to be the mz and time columns concatenated by "_"
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        mutate(mz_time = paste0(mz, "_", time)) 
    
    # update the mean_intensity column
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        select(-mean_intensity) %>%
        left_join(met_raw_wide_original_mean_intensity, by = c("mz_time")) %>%
        ## round mean_intensity to integers
        mutate(mean_intensity = round(mean_intensity, 0)) %>%
        ## move mean_intensity column right after mz_time
        relocate(mean_intensity, .after = "mz_time")
    
    # round mz_isotope and time_isotope to 4 decimal places and integers
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic$mz_isotope = round(MSMICA_results_filtered_final_4_secondary_adduct_isotopic$mz_isotope, 4)
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time_isotope = round(MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time_isotope, 0)

    # left join the MSMICA_results_filtered_final_4_secondary_adduct_isotopic with the met_raw_wide_original by mz and time to append the intensity columns
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = left_join(MSMICA_results_filtered_final_4_secondary_adduct_isotopic, met_raw_wide_original, by = c("mz" = "mz", "time" = "time"))

    # remove duplicates
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        distinct()
    
    # arrange the MSMICA_results_filtered_final_4_secondary_adduct_isotopic first by HMDB_ID, then by KEGG_ID
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic = MSMICA_results_filtered_final_4_secondary_adduct_isotopic %>%
        arrange(HMDB_ID, KEGG_ID)

    # round time_predicted and time_difference to 1 decimal place
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time_predicted = round(MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time_predicted, 1)
    MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time_difference = round(MSMICA_results_filtered_final_4_secondary_adduct_isotopic$time_difference, 1)

    # create mz_isotope and time_isotope columns in MSMICA_results_filtered_final_4
    MSMICA_results_filtered_final_4_isotopic = MSMICA_results_filtered_final_4 %>%
        mutate(mz_isotope = NA_real_, time_isotope = NA_real_) %>%
        relocate(mz_isotope, .after = "time") %>%
        relocate(time_isotope, .after = "mz_isotope")

    # combine MSMICA_results_filtered_final_4_secondary_adduct_isotopic and MSMICA_results_filtered_final_4_isotopic and remove duplicates
    MSMICA_results_filtered_adduct_isotopic_final = bind_rows(MSMICA_results_filtered_final_4_secondary_adduct_isotopic, MSMICA_results_filtered_final_4_isotopic)
    
    # arrange the MSMICA_results_filtered_adduct_isotopic_final first by HMDB_ID, then by KEGG_ID
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        arrange(KEGG_ID, HMDB_ID)

    # remove duplicates based on InChIKey, Name, Adduct, mz, time
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        distinct(InChIKey, Name, Adduct, mz, time, .keep_all = TRUE)

    # if adduct_type is NA", set it to "primary"
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        mutate(adduct_type = ifelse(is.na(adduct_type), "primary", adduct_type))

    # update the identification_type column to indicate which Schymanski level should be used for each annotation based on identification_method column
    ## if identification_method equal to "mz matching; retention time prediction" or "mz matching; retention time prediction; clustering of adducts and isotopes", then it should be "Schymanski level 3b" because it only has mass spectromery evidence (m/z matching) and chromatographic evidence from RT prediction and/or clustering of adducts and isotopes
    ## if identification_method equal to other cases where more pieces of evidence are used, then it should be "Schymanski level 3a" because it has mass spectromery evidence, chromatographic evidence, and biological evidence (biospecimen-specific concentration and/or precursor-product/transporter correlation)
    ## special case: if adduct_type is "secondary", then it should be "Schymanski level 4" because it is secondary adduct from the MS
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        mutate(
            identification_type = case_when(
                identification_method == "mz matching; retention time prediction" | identification_method == "mz matching; retention time prediction; clustering of adducts and isotopes" ~ "Schymanski level 3b",
                TRUE ~ "Schymanski level 3a"
            )
        ) %>%
        mutate(
            identification_type = case_when(
                adduct_type == "secondary" ~ "Schymanski level 4",
                TRUE ~ identification_type
            )
        )

    # recalculate the time_difference using the time_predicted and time columns
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        mutate(time_difference = abs(time_predicted - time))

    # recalculate metabolite_within_feature_number
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        group_by(mz_time) %>%
        ## summarize the number of metabolites within each feature
        mutate(metabolite_within_feature_number = n()) %>%
        ungroup()

    # add a new column, isomer_exist, which is TRUE when metabolite_within_feature_number is greater than 1
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        mutate(isomer_exist = ifelse(metabolite_within_feature_number > 1, TRUE, FALSE))

    # reoreder the column names for better readability
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        dplyr::select(mz, time, mz_time, mz_isotope, time_isotope, Adduct, adduct_type, Name, identification_type, identification_method, Probability, metabolite_within_feature_number, isomer_exist, time_predicted, time_difference, mean_intensity, everything())

    # arrange the MSMICA_results_filtered_adduct_isotopic_final first by HMDB_ID, then by KEGG_ID
    MSMICA_results_filtered_adduct_isotopic_final = MSMICA_results_filtered_adduct_isotopic_final %>%
        arrange(KEGG_ID, HMDB_ID)

    # if save_unidentified is not FALSE
    if (save_unidentified != FALSE){
        # select those features not identified by MSMICA by using the met_raw_wide_original data frame with the anti_join function
        features_not_identified = met_raw_wide_original %>% 
            ## based on mz and time columns
            anti_join(MSMICA_results_filtered_adduct_isotopic_final, by = c("mz" = "mz", "time" = "time")) %>%
            ## based on mz_isotope and time_isotope columns
            anti_join(MSMICA_results_filtered_adduct_isotopic_final, by = c("mz" = "mz_isotope", "time" = "time_isotope"))
        
        # within features_not_identified, remove those features in the mz_time_adduct_isotope (features that ever identified in the clustering of adducts and isotopes step)
        features_not_identified = features_not_identified %>%
            mutate(mz_time_sample = paste0(mz, "_", time)) %>%
            filter(!(mz_time_sample %in% mz_time_adduct_isotope)) %>%
            dplyr::select(-mz_time_sample)

        # save the features_not_identified data frame
        if (prefix != ""){
            # add the prefix to saved file names
            write_csv(features_not_identified, paste0(prefix, "_features_not_identified.csv"), progress = FALSE)
        } else {
            # save the data without prefix
            write_csv(features_not_identified, "features_not_identified.csv", progress = FALSE)
        }

        # clear memory
        rm(features_not_identified)
        gc()
    }

    # save the MSMICA_results_filtered_adduct_isotopic_final data frame
    if (prefix != ""){
        # add the prefix to saved file names
        write_csv(MSMICA_results_filtered_adduct_isotopic_final, paste0(prefix, "_MSMICA_identified_metabolites_filtered_secondary_adduct_isotopic.csv"), progress = FALSE)
    } else {
        # save the data without prefix
        write_csv(MSMICA_results_filtered_adduct_isotopic_final, "identified_MSMICA_identified_metabolites_filtered_secondary_adduct_isotopic.csv", progress = FALSE)
    }

    # add a new column, isomer_exist, which is TRUE when metabolite_within_feature_number is greater than 1
    MSMICA_results_filtered_final_5 = MSMICA_results_filtered_final_5 %>%
        mutate(isomer_exist = ifelse(metabolite_within_feature_number > 1, TRUE, FALSE))

    # reoreder the column names for better readability
    MSMICA_results_filtered_final_5 = MSMICA_results_filtered_final_5 %>%
        dplyr::select(mz, time, mz_time, Adduct, Name, identification_type, identification_method, Probability, metabolite_within_feature_number, isomer_exist, time_predicted, time_difference, mean_intensity, everything())

    if (prefix != ""){
        # add the prefix to saved file names
        write_csv(MSMICA_results_filtered_final_5, paste0(prefix, "_MSMICA_identified_metabolites_filtered.csv"), progress = FALSE)
    } else {
        # save the data without prefix
        write_csv(MSMICA_results_filtered_final_5, "identified_MSMICA_identified_metabolites_filtered.csv", progress = FALSE)
    }

    # print where the MSMICA_TCA_result_identified.csv and MSMICA_TCA_result_not_identified.csv are saved
    print("The identified_MSMICA_identified_metabolites.csv was saved in:")
    print(getwd())

    # log the end time of the MSMICA algorithm
    end_time = Sys.time()
    message("The MSMICA algorithm ends at: ", end_time)
    message("The total running time of the MSMICA algorithm is: ", end_time - start_time)

    # return MSMICA main result data frame in case it is directly needed for further analysis
    return(MSMICA_results_filtered_adduct_isotopic_final)
}
