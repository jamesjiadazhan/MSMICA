#' Bayesian identification probability for adduct/isotope clusters
#'
#' For a group of candidate feature annotations sharing the same
#' monoisotopic mass, this function computes the posterior probability
#' that each feature is the true representative of its metabolite by
#' combining three sources of experimental evidence:
#' \enumerate{
#'   \item \strong{Adduct formation} -- whether the observed adduct is
#'     the dominant one (\code{M+H} or \code{M-H}) or a secondary form.
#'   \item \strong{Adduct correlation} -- whether this feature's
#'     intensity co-varies with another feature of the same monoisotopic
#'     mass at a nearby retention time (Spearman r ≥
#'     \code{adduct_corr_r_thresh}).
#'   \item \strong{Isotopologue correlation} -- whether a
#'     \code{[M+H]+[+1]} (or equivalent) peak co-elutes, has the
#'     expected relative abundance (≤ 100 % of the monoisotopic peak),
#'     and shows Spearman r ≥ \code{isotopic_corr_r_thresh}.
#' }
#' The likelihood ratios for all three evidence sources are derived from
#' experimentally validated MSMICA training data (hard-coded contingency
#' counts inside the child function). A uniform prior over the \emph{k}
#' candidate features in the cluster is updated with these log-likelihood
#' ratios to yield a posterior probability for each feature.
#'
#' Groups with no adduct correlation and no isotopologue match return an
#' empty data frame immediately, since such groups offer no clustering
#' evidence and are handled by the main-adduct local optimization instead.
#'
#' Requires that the following objects be defined in the calling scope
#' (set by \code{MSMICA_algorithm()}): \code{adduct_corr_time_thresh},
#' \code{adduct_corr_r_thresh}, \code{isotopic_corr_time_thresh},
#' \code{isotopic_corr_r_thresh},
#' \code{primary_df_simple}, and
#' \code{isotope_df_simple}.
#'
#' @param data Data frame of candidate adduct annotations for a single
#'   monoisotopic mass group, with columns including
#'   \code{Mono_mass}, \code{Adduct_annotated}, \code{mz_annotated},
#'   \code{time_annotated}, and \code{mz_time_annotated}.
#' @param MSMICA_cor_input Wide-format intensity matrix used for
#'   correlation calculations, with sample rows and feature (\code{mz_time})
#'   columns.
#' @return A data frame (same schema as \code{data}) augmented with
#'   \code{MSMICA_identification} (1 = MAP candidate, 0 = other) and
#'   \code{Probability} (0-100 posterior percentage), sorted by
#'   descending \code{Probability}. Returns an empty data frame if no
#'   adduct or isotopologue clustering evidence is present.
#' @keywords internal
#' @noRd
bayesian_probability_calculation_cluster = function(data, MSMICA_cor_input,
                                                    adduct_corr_time_thresh,
                                                    adduct_corr_r_thresh,
                                                    isotopic_corr_time_thresh,
                                                    isotopic_corr_r_thresh,
                                                    primary_df_simple,
                                                    isotope_df_simple) {
    # All rows in data share one monoisotopic mass
    current_monomass = unique(data$Mono_mass)

    ####################### adduct correlation analysis (correlation between different adducts with the same monoisotopic mass) #######################
    # this version is more efficient than the previous version because it groups the features by retention time and then calculates the correlation matrix for each group
    # Compute the adduct correlation matrix only if there are more than 1 rows in the data data frame
    if(nrow(data) > 1) {

        # Order data by retention time to enable efficient grouping
        data_ordered = data %>%
            arrange(time_annotated)
        
        # Group features with close retention times and assign them to the time_groups list by the mz_time_annotated column
        time_groups = list()

        # Only group the features if there are more than 1 row in the data_ordered data frame
        if (nrow(data_ordered) > 0) {
            current_group = data_ordered$mz_time_annotated[1]
            last_time = data_ordered$time_annotated[1]
            
            for (i in 2:nrow(data_ordered)) {
                if ((data_ordered$time_annotated[i] - last_time) <= adduct_corr_time_thresh) {
                    current_group = c(current_group, data_ordered$mz_time_annotated[i])
                } else {
                    if (length(current_group) > 1) {
                        time_groups[[length(time_groups) + 1]] = current_group
                    }
                    current_group = data_ordered$mz_time_annotated[i]
                }
                last_time = data_ordered$time_annotated[i]
            }
            if (length(current_group) > 1) {
                time_groups[[length(time_groups) + 1]] = current_group
            }
        }

        # Function to find high correlation pairs
        find_high_correlation_pairs = function(cor_matrix, threshold) {
            # Filter the correlation matrix by setting the lower triangle and diagonal to NA
            cor_matrix[lower.tri(cor_matrix)] = NA
            diag(cor_matrix) = NA

            # Find the indices of high correlation pairs using the given threshold
            high_cor_indices = which(cor_matrix >= threshold, arr.ind = TRUE)

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

            ## for example, here the correlation between M-H2O+H_157.109_208.2 and C00062_M+H_175.1195_210 is 0.51, so they are a high correlation pair. Thus, when we extract the column name 2 and 3, this is equal to extract the row name 2 and column name 3 in the original cor_matrix, as the row and column names are the same in the cor_matrix.

            high_cor_pairs = apply(high_cor_indices, 1, function(idx) {
                colnames(cor_matrix)[idx]
            })

            # Convert to unique pairs
            unique_pairs = unique(apply(high_cor_pairs, 2, function(pair) {
                sort(pair)
            }, simplify = FALSE))

            # Convert the list of unique pairs to a list of character vectors
            unique_pairs_list = lapply(unique_pairs, function(pair) {
                return(pair)
            })

            return(unique_pairs_list)
        }

        # use lapply to calculate the correlation matrix for each time group
        filtered_original_groups = lapply(time_groups, function(group) {
            # extract the intensity values for the current group of features with close retention times
            MSMICA_cor_input_annotated = MSMICA_cor_input[, group]
            # calculate the correlation matrix for the current group of features with close retention times
            full_cor_matrix_annotated = suppressWarnings(cor(MSMICA_cor_input_annotated, MSMICA_cor_input_annotated, method="spearman", use="pairwise.complete.obs"))
            # find the high correlation pairs for the current group of features with close retention times
            intra_correlation_all = find_high_correlation_pairs(full_cor_matrix_annotated, adduct_corr_r_thresh)
            # return the correlation pairs that are greater than the adduct_corr_r_thresh
            return(intra_correlation_all)
        })

        # assign the filtered_original_groups to the intra_correlation_all_clean
        intra_correlation_all_clean = filtered_original_groups 

        # Flatten the nested list
        flat_list = do.call(c, intra_correlation_all_clean)

        # initialize the intra_correlation_all_clean_df as NULL
        intra_correlation_all_clean_df = NULL 

        # check if any intra connection > adduct_corr_r_thresh found. If yes, then perform the following. If no, then return NULL
        if (length(flat_list) > 0) {
            # Transform to data frame
            intra_correlation_all_clean_df = data.frame(cluster = rep(1:length(flat_list), sapply(flat_list, length)),
                variable = unlist(flat_list))

            # rename the cluster column as adduct_corr_cluster
            intra_correlation_all_clean_df = intra_correlation_all_clean_df %>%
                rename(adduct_corr_cluster = cluster)

            # Merge the data with the intra_correlation_all_clean_df to add the adduct_corr_cluster column and assign values for those within the adduct correlation cluster
            data_adduct_corr = data %>%
                full_join(intra_correlation_all_clean_df, by = c("mz_time_annotated" = "variable"), relationship = "many-to-many")

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

    # remove duplicate rows in the data_adduct_corr data frame based on Mono_mass, Adduct, mz, and time
    data_adduct_corr_2 = data_adduct_corr %>%
        distinct(Mono_mass, Adduct_annotated, mz_annotated, time_annotated, .keep_all = TRUE) %>%
        # replace the numeric values in the adduct_corr_cluster with TRUE
        mutate(adduct_corr_cluster = ifelse(!is.na(adduct_corr_cluster), TRUE, NA))

    ####################### isotopic adduct correlation analysis (correlation between regular adducts and their most abundant isotologues) #######################
    # Initialize an empty data frame to store the ratios
    abundance_ratios = data.frame()

    # Filter both datasets for the current Mono_mass
    primary_df = primary_df_simple %>%
        dplyr::filter(Mono_mass %in% current_monomass)

    isotope_df = isotope_df_simple %>%
        dplyr::filter(Mono_mass %in% current_monomass)

    # if there is no isotopic adduct, then add NA to the mz_isotope and time_isotope columns
    if (nrow(isotope_df) == 0) {
        data_adduct_corr_3 = data_adduct_corr_2 %>%
            mutate(
                mz_isotope = NA,
                time_isotope = NA
                )
    } else if (any(primary_df$Adduct_annotated %in% isotope_df$Adduct_annotated) == TRUE) {
        # if the isotopic adduct annotation at least has 1 same Adduct_annotated in the primary adduct, then perform the following
        ## loop over each row in primary_df to calculate isotopic abundance ratio and correlation if there is any isotopic adduct
        for (i in 1:nrow(primary_df)) {

            # get the current primary adduct and time
            current_primary_df = primary_df[i,]
            adduct_primary = current_primary_df$Adduct_annotated
            time_primary = current_primary_df$time_annotated
            mz_time_primary = current_primary_df$mz_time_annotated

            # find the isotopic adduct that has the same adduct in the primary adduct and the time difference is within the isotopic_corr_time_thresh
            isotope_df_filtered = isotope_df %>%
                filter(Adduct_annotated == adduct_primary) %>%
                filter(time_annotated >= time_primary - isotopic_corr_time_thresh & time_annotated <= time_primary + isotopic_corr_time_thresh)

            # if the current_primary_df does not have any isotopic adduct, skip the current iteration
            if (nrow(isotope_df_filtered) == 0) {
                next
            }

            ## loop over each row in isotope_df_filtered
            for (j in 1:nrow(isotope_df_filtered)) {

                current_isotope_df = isotope_df_filtered[j,]
                adduct_isotope = current_isotope_df$Adduct_annotated
                time_isotope = current_isotope_df$time_annotated
                mz_time_isotope = current_isotope_df$mz_time_annotated

                # extract the sample intensities for the current primary and isotopic adducts from the MSMICA_cor_input data frame
                primary_adduct_intensity = MSMICA_cor_input[, mz_time_primary]
                isotope_adduct_intensity = MSMICA_cor_input[, mz_time_isotope]

                ## calculate the ratio of the isotopic adduct to the primary adduct (because the intensity is in log2 transformed, so the ratio is 2^(isotope_adduct_intensity - primary_adduct_intensity))
                ratios = 2^(isotope_adduct_intensity - primary_adduct_intensity)

                ## create a dataframe to store the results
                result_df = data.frame(current_monomass = current_monomass, mz_primary = current_primary_df$mz_annotated, time_primary = current_primary_df$time_annotated, mz_isotope = current_isotope_df$mz_annotated, time_isotope = current_isotope_df$time_annotated, Adduct = adduct_primary, primary = primary_adduct_intensity, isotope = isotope_adduct_intensity, ratios = ratios)
                
                # rename the column names:  current_monomass mz_primary time_primary mz_isotope time_isotope Adduct primary isotope ratios
                colnames(result_df) = c("current_monomass", "mz_primary", "time_primary", "mz_isotope", "time_isotope", "Adduct", "primary", "isotope", "ratios")

                ## add the result to the abundance_ratios dataframe
                abundance_ratios = rbind(abundance_ratios, result_df)
            }
        }

        # if there is no isotopic adduct, then add NA to the mz_isotope and time_isotope columns
        if (nrow(abundance_ratios) == 0) {
            data_adduct_corr_3 = data_adduct_corr_2 %>%
                mutate(
                    mz_isotope = NA,
                    time_isotope = NA
                    )
        } else {
            # Convert it back to a tibble for nicer messageing, and set row names as a new column
            abundance_ratios_2 = as_tibble(abundance_ratios)

            # Convert the ratios to numeric, and replace any non-finite values with 0
            abundance_ratios_2$ratios[!is.finite(abundance_ratios_2$ratios)] = NA

            # replace all 0 values with NA
            abundance_ratios_2$ratios[abundance_ratios_2$ratios == 0] = NA

            # multiply the ratios by 100
            abundance_ratios_2$ratios = abundance_ratios_2$ratios * 100

            # Convert all 0 to NA
            abundance_ratios_2$ratios[abundance_ratios_2$ratios == 0] = NA

            # summarize the abundance ratio using mean by current_monomass
            abundance_ratios_3 = abundance_ratios_2 %>%
                group_by(current_monomass, mz_primary, time_primary, mz_isotope, time_isotope, Adduct) %>%
                ## calculate the spearman correlation between the isotopic adduct and the primary adduct
                mutate(
                    correlation = cor(primary, isotope, method = "spearman", use = "pairwise.complete.obs")
                    ) %>%
                summarize(
                    mean_abundance_ratio = mean(ratios, na.rm = TRUE),
                    mean_correlation = mean(correlation, na.rm = TRUE), 
                    .groups = "keep"
                ) %>%
                ungroup()

            # remove the abundance_ratios_3 with mean_abundance_ratio > 100 if any (since the isotopic adduct should have lower intensity than the primary adduct overall)
            abundance_ratios_4 = abundance_ratios_3 %>%
                filter(mean_abundance_ratio <= 100) %>%
                ## remove mean_correlation < isotopic_corr_r_thresh (default is 0.71)
                filter(mean_correlation >= isotopic_corr_r_thresh)

            if (nrow(abundance_ratios_4) > 0) {
                # select only these columns: mz_primary, time_primary, current_monomass, mz_isotope, time_isotope, Adduct, mean_correlation
                abundance_ratios_4 = abundance_ratios_4[, c("mz_primary", "time_primary", "current_monomass", "mz_isotope", "time_isotope", "Adduct", "mean_correlation")]

                # rename the mean_correlation as correlation
                abundance_ratios_4 = abundance_ratios_4 %>%
                    rename(correlation = mean_correlation)

                # round the correlation to 2 decimal places
                abundance_ratios_4$correlation = round(abundance_ratios_4$correlation, 2)

                # round mz_primary to 4 decimal places
                abundance_ratios_4$mz_primary = round(abundance_ratios_4$mz_primary, 4)

                # round mz_annotated to 4 decimal places in the data_adduct_corr_2
                data_adduct_corr_2$mz_annotated = round(data_adduct_corr_2$mz_annotated, 4)

                # round time_primary to 0 decimal place
                abundance_ratios_4$time_primary = round(abundance_ratios_4$time_primary, 0)

                # round time_annotated to 1 decimal place in the data_adduct_corr_2
                data_adduct_corr_2$time_annotated = round(data_adduct_corr_2$time_annotated, 0)

                # inner join the abundance_ratios_4 with the data_adduct_corr_2 to add the abundance_ratios_4 columns to the data_adduct_corr_2 by current_monomass and Adduct_annotated
                data_adduct_corr_3 = data_adduct_corr_2 %>%
                    left_join(abundance_ratios_4, by = c("mz_annotated"="mz_primary", "time_annotated"="time_primary", "Mono_mass"="current_monomass", "Adduct_annotated"="Adduct"), relationship = "many-to-many")

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
    } else {
        # else, if there is no isotopic adduct that has the same adduct in the primary adduct, then add NA to the mz_isotope and time_isotope columns
        data_adduct_corr_3 = data_adduct_corr_2 %>%
            mutate(
                mz_isotope = NA,
                time_isotope = NA
                )
    }

    # remove duplicates based on Mono_mass, Adduct_annotated, mz_annotated, and time_annotated
    data_adduct_corr_3 = data_adduct_corr_3 %>%
        distinct(Mono_mass, Adduct_annotated, mz_annotated, time_annotated, .keep_all = TRUE)

    # if there is no adduct correlation cluster and isotope correlation cluster (all adduct_corr_cluster column's value is NA or all mz_isotope's value is NA), then just return an empty data frame without even calculating the Bayesian probability because all we care about this function is to calculate the Bayesian probability for those metabolite adducts with adduct correlation cluster and isotope correlation cluster
    if (all(is.na(data_adduct_corr_3$adduct_corr_cluster)) & all(is.na(data_adduct_corr_3$mz_isotope))) {
        return(data.frame())
    }

    # select features not with NA in the adduct_corr_cluster OR mz_isotope columns in the data_adduct_corr_3
    ## this is needed because in long LC-MS run (15-20 mins), there may be too many features with very similar m/z and different retention time. As MSMICA Bayesian probability calculation is based on 2^n, where n is the number of features with adducts assigned to a given metabolite, when n is larger than 13, that will be too large to calculation for small computers. In the end, features with adduct correlation or isotope correlation will have the highest probability. Thus, we simplify the calculation by removing features with no adduct correlation or isotope correlation.
    data_adduct_corr_3 = data_adduct_corr_3 %>%
        filter(!is.na(adduct_corr_cluster) | !is.na(mz_isotope))

    # see if more than 10 rows are present in the data_adduct_corr_3. If so, prioritize those features with both adduct correlation AND isotope correlation.
    if (nrow(data_adduct_corr_3) > 10) {
        data_adduct_corr_3_filtered = data_adduct_corr_3 %>%
            filter(!is.na(adduct_corr_cluster) & !is.na(mz_isotope))

        # if there is no row with adduct_corr_cluster and mz_isotope, then arrange the data by putting the "M+H" or "M-H" adduct at the top and then select the top 10 rows
        if (nrow(data_adduct_corr_3_filtered) == 0) {
            data_adduct_corr_3_filtered = data_adduct_corr_3 %>%
                arrange(desc(Adduct_annotated == "M+H" | Adduct_annotated == "M-H"))
            
            data_adduct_corr_3_filtered = data_adduct_corr_3_filtered[1:10,]
        }

        # update the data_adduct_corr_3 with the data_adduct_corr_3_filtered
        data_adduct_corr_3 = data_adduct_corr_3_filtered
    }

    # remove duplicates based on Mono_mass and mz_time_annotated
    data_adduct_corr_3 = data_adduct_corr_3 %>%
        distinct(Mono_mass, mz_time_annotated, .keep_all = TRUE)

    # if there is only 1 row, then return the following data as the output
    if (nrow(data_adduct_corr_3) == 1) {
        final_result_cluster = data_adduct_corr_3 %>%
            mutate(
                MSMICA_identification = 1,
                Probability = 100
            )
        return(final_result_cluster)
    }

    # calculate the Bayesian posterior probability using the following evidence:
    ## prior probability: before considering any evidence, the probability of the feature best representing the metabolite is 1/k, where k is the number of features in the clustering result
    ## evidence 1: adduct formation
    ## evidence 2: adduct correlation (correlation between different adducts with the same monoisotopic mass (e.g. M+H and M+Na))
    ## evidence 3: isotopic adduct correlation (correlation between regular adducts and their most abundant isotopologue (e.g. M+H and M+H[+1]))
    bayesian_probability_calculation_cluster_child = function(
        data_cluster,
        n_MH_right = 1259,
        n_MH_total = 1883,
        n_other_right = 1144,
        n_other_total = 6361,
        n_cor_yes_right = 1156,
        n_cor_yes_total = 2121,
        n_cor_no_right = 1235,
        n_cor_no_total = 6099,
        n_iso_yes_right = 537,
        n_iso_yes_total = 846,
        n_iso_no_right = 2431,
        n_iso_no_total = 10149,
        eps = 1e-12,
        laplace = 0
    ) {
        K = nrow(data_cluster)
        if (K == 0) return(NULL)

        # Keep probabilities away from 0 and 1 so log-ratios stay finite.
        clamp = function(p) pmin(pmax(p, eps), 1 - eps)

        # Recover counts for features that were not the correct representative.
        n_MH_wrong       = n_MH_total - n_MH_right
        n_other_wrong    = n_other_total - n_other_right
        n_cor_yes_wrong  = n_cor_yes_total - n_cor_yes_right
        n_cor_no_wrong   = n_cor_no_total - n_cor_no_right
        n_iso_yes_wrong  = n_iso_yes_total - n_iso_yes_right
        n_iso_no_wrong   = n_iso_no_total - n_iso_no_right

        add = laplace

        # Estimate P(feature state | H) and P(feature state | not H) for
        # each evidence source, where H means the row is the true parent ion.
        p_A_MH_given_H =
            clamp((n_MH_right + add) / (n_MH_right + n_other_right + 2 * add))
        p_A_other_given_H =
            clamp((n_other_right + add) / (n_MH_right + n_other_right + 2 * add))
        p_A_MH_given_notH =
            clamp((n_MH_wrong + add) / (n_MH_wrong + n_other_wrong + 2 * add))
        p_A_other_given_notH =
            clamp((n_other_wrong + add) / (n_MH_wrong + n_other_wrong + 2 * add))

        p_C_yes_given_H =
            clamp((n_cor_yes_right + add) / (n_cor_yes_right + n_cor_no_right + 2 * add))
        p_C_no_given_H =
            clamp((n_cor_no_right + add) / (n_cor_yes_right + n_cor_no_right + 2 * add))
        p_C_yes_given_notH =
            clamp((n_cor_yes_wrong + add) / (n_cor_yes_wrong + n_cor_no_wrong + 2 * add))
        p_C_no_given_notH =
            clamp((n_cor_no_wrong + add) / (n_cor_yes_wrong + n_cor_no_wrong + 2 * add))

        p_I_yes_given_H =
            clamp((n_iso_yes_right + add) / (n_iso_yes_right + n_iso_no_right + 2 * add))
        p_I_no_given_H =
            clamp((n_iso_no_right + add) / (n_iso_yes_right + n_iso_no_right + 2 * add))
        p_I_yes_given_notH =
            clamp((n_iso_yes_wrong + add) / (n_iso_yes_wrong + n_iso_no_wrong + 2 * add))
        p_I_no_given_notH =
            clamp((n_iso_no_wrong + add) / (n_iso_yes_wrong + n_iso_no_wrong + 2 * add))

        # Turn each row into three binary evidence indicators.
        A_is_MH = data_cluster$Adduct_annotated %in% c("M+H", "M-H")
        C_yes = !is.na(data_cluster$adduct_corr_cluster) &
                data_cluster$adduct_corr_cluster == TRUE
        I_yes = !is.na(data_cluster$mz_isotope)

        # Start with a uniform prior over all rows in the cluster.
        prior = rep(1 / K, K)

        # Convert each evidence source to a log-likelihood ratio so positive
        # values favor H and negative values favor not H.
        log_A = ifelse(
            A_is_MH,
            log(p_A_MH_given_H / p_A_MH_given_notH),
            log(p_A_other_given_H / p_A_other_given_notH)
        )

        log_C = ifelse(
            C_yes,
            log(p_C_yes_given_H / p_C_yes_given_notH),
            log(p_C_no_given_H / p_C_no_given_notH)
        )

        log_I = ifelse(
            I_yes,
            log(p_I_yes_given_H / p_I_yes_given_notH),
            log(p_I_no_given_H / p_I_no_given_notH)
        )

        # Combine the prior and all evidence terms in log space.
        log_w = log(prior) + log_A + log_C + log_I

        # Stabilize exponentiation, then normalize weights into probabilities
        # that sum to 100% within the cluster.
        log_w_centered = log_w - max(log_w, na.rm = TRUE)
        w = exp(log_w_centered)
        post = w / sum(w, na.rm = TRUE)

        out = data_cluster
        out$Probability = 100 * post
        # Keep the intermediate log terms for diagnostics and interpretation.
        out$log_prior = log(prior)
        out$log_A = log_A
        out$log_C = log_C
        out$log_I = log_I
        out$log_w = log_w

        # Mark every row tied for the highest posterior probability.
        max_prob = max(out$Probability, na.rm = TRUE)
        tol = 1e-8
        out$MSMICA_identification = as.integer(abs(out$Probability - max_prob) < tol)

        # Return the highest-probability candidates first.
        out = out[order(-out$Probability), ]
        rownames(out) = NULL
        out
    }

    # apply the bayesian_probability_calculation_cluster_child function to the data_adduct_corr_3 data frame
    final_result_cluster = bayesian_probability_calculation_cluster_child(data_cluster = data_adduct_corr_3)

    # remove the log_prior, log_A, log_C, log_I, and log_w columns
    final_result_cluster = final_result_cluster %>%
        dplyr::select(-log_prior, -log_A, -log_C, -log_I, -log_w)

    return(final_result_cluster)
}
