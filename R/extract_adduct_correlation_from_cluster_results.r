#' Extract pairwise adduct correlations from clustered MSMICA results
#'
#' After clustering features with the same monoisotopic mass and close
#' retention time (the \code{adduct_corr_cluster == TRUE} rows of
#' \code{final_results_cluster}), this helper reports every pairwise
#' intensity correlation between the features that belong to the same
#' RT-local cluster of the same monoisotopic mass. The cluster-within-
#' cluster structure is built explicitly by splitting first on
#' \code{Mono_mass} and then on RT gaps, so features belonging to
#' different chromatographic peaks of the same mass do not cross-
#' correlate. The output is intended for the adduct clustering report
#' in the final MSMICA results.
#'
#' @param final_results_cluster MSMICA per-feature result table carrying
#'   an \code{adduct_corr_cluster} flag.
#' @param MSMICA_cor_input Intensity matrix used for correlation
#'   (columns keyed by \code{mz_time_annotated}).
#' @param time_threshold Numeric. Maximum retention time gap (in
#'   seconds) between adjacent features in the same RT-local cluster.
#'   Defaults to the global \code{adduct_correlation_time_threshold}
#'   set by \code{MSMICA_algorithm()}.
#' @param cor_method Character. Correlation method passed to
#'   \code{cor()}. Defaults to \code{"spearman"}.
#' @return A data frame with one row per ordered-pair adduct
#'   correlation, annotated with both features' metadata.
#' @keywords internal
#' @noRd
extract_adduct_correlation_from_cluster_results = function(final_results_cluster, MSMICA_cor_input, time_threshold = adduct_correlation_time_threshold, cor_method = "spearman") {
    # Keep only rows that the clustering step marked as adduct-correlated
    final_results_cluster_corr = final_results_cluster %>%
        dplyr::filter(!is.na(adduct_corr_cluster) & adduct_corr_cluster == TRUE)

    if (nrow(final_results_cluster_corr) == 0) {
        return(data.frame())
    }

    # split by Mono_mass
    results_split = split(final_results_cluster_corr, final_results_cluster_corr$Mono_mass)

    # helper: split one mono mass into local RT groups
    split_by_rt_gap = function(df, time_threshold) {
        df = df %>%
            dplyr::arrange(time_annotated)

        if (nrow(df) == 0) {
            return(list())
        }

        groups = list()
        current_idx = 1
        current_group = 1

        if (nrow(df) == 1) {
            groups[[1]] = df
            return(groups)
        }

        for (i in 2:nrow(df)) {
            if ((df$time_annotated[i] - df$time_annotated[i - 1]) <= time_threshold) {
                current_idx = c(current_idx, i)
            } else {
                groups[[current_group]] = df[current_idx, , drop = FALSE]
                current_group = current_group + 1
                current_idx = i
            }
        }

        groups[[current_group]] = df[current_idx, , drop = FALSE]
        groups
    }

    # helper for one RT-local cluster within one Mono_mass
    extract_one_rt_cluster = function(df_one, rt_cluster_id) {
        if (nrow(df_one) < 2) {
            return(NULL)
        }

        feature_names = intersect(df_one$mz_time_annotated, colnames(MSMICA_cor_input))

        if (length(feature_names) < 2) {
            return(NULL)
        }

        df_one2 = df_one %>%
            dplyr::filter(mz_time_annotated %in% feature_names) %>%
            dplyr::distinct(mz_time_annotated, .keep_all = TRUE)

        if (nrow(df_one2) < 2) {
            return(NULL)
        }

        cor_mat = suppressWarnings(
            cor(
                MSMICA_cor_input[, df_one2$mz_time_annotated, drop = FALSE],
                method = cor_method,
                use = "pairwise.complete.obs"
            )
        )

        if (length(dim(cor_mat)) != 2 || nrow(cor_mat) < 2) {
            return(NULL)
        }

        idx = which(upper.tri(cor_mat), arr.ind = TRUE)

        if (nrow(idx) == 0) {
            return(NULL)
        }

        out = data.frame(
            Mono_mass = unique(df_one2$Mono_mass)[1],
            rt_cluster_id = rt_cluster_id,
            mz_time_annotated_1 = rownames(cor_mat)[idx[, 1]],
            mz_time_annotated_2 = colnames(cor_mat)[idx[, 2]],
            adduct_correlation = cor_mat[idx],
            stringsAsFactors = FALSE
        )

        meta1 = df_one2 %>%
            dplyr::select(
                mz_time_annotated,
                identification_type_annotated,
                ion_mode,
                Adduct_annotated,
                mz_annotated,
                time_annotated,
                mean_intensity,
                MSMICA_identification,
                Probability
            ) %>%
            dplyr::rename(
                mz_time_annotated_1 = mz_time_annotated,
                identification_type_annotated_1 = identification_type_annotated,
                ion_mode_1 = ion_mode,
                Adduct_annotated_1 = Adduct_annotated,
                mz_annotated_1 = mz_annotated,
                time_annotated_1 = time_annotated,
                mean_intensity_1 = mean_intensity,
                MSMICA_identification_1 = MSMICA_identification,
                Probability_1 = Probability
            )

        meta2 = df_one2 %>%
            dplyr::select(
                mz_time_annotated,
                identification_type_annotated,
                ion_mode,
                Adduct_annotated,
                mz_annotated,
                time_annotated,
                mean_intensity,
                MSMICA_identification,
                Probability
            ) %>%
            dplyr::rename(
                mz_time_annotated_2 = mz_time_annotated,
                identification_type_annotated_2 = identification_type_annotated,
                ion_mode_2 = ion_mode,
                Adduct_annotated_2 = Adduct_annotated,
                mz_annotated_2 = mz_annotated,
                time_annotated_2 = time_annotated,
                mean_intensity_2 = mean_intensity,
                MSMICA_identification_2 = MSMICA_identification,
                Probability_2 = Probability
            )

        out = out %>%
            dplyr::left_join(meta1, by = "mz_time_annotated_1") %>%
            dplyr::left_join(meta2, by = "mz_time_annotated_2") %>%
            dplyr::arrange(
                Mono_mass,
                rt_cluster_id,
                dplyr::desc(adduct_correlation)
            )

        rownames(out) = NULL
        out
    }

    # loop through each Mono_mass, then each RT-local cluster
    out_list = lapply(results_split, function(df_mass) {
        rt_groups = split_by_rt_gap(df_mass, time_threshold)

        rt_out = lapply(seq_along(rt_groups), function(k) {
            extract_one_rt_cluster(rt_groups[[k]], rt_cluster_id = k)
        })

        dplyr::bind_rows(rt_out)
    })

    out_all = dplyr::bind_rows(out_list)
    rownames(out_all) = NULL
    out_all
}