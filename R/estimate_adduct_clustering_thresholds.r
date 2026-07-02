#' Estimate adduct clustering thresholds from high-confidence anchors
#'
#' Uses single-match primary-adduct features as anchors, pairs them with
#' same-metabolite secondary-adduct candidate features, and chooses RT and
#' correlation thresholds that retain the requested fraction of empirical
#' anchor-secondary pairs.
#'
#' @keywords internal
#' @noRd
estimate_adduct_clustering_thresholds = function(primary_anchor_data,
                                                 annotated_adduct_data,
                                                 cor_input,
                                                 capture_fraction = 0.8,
                                                 primary_adducts = c("M+H", "M-H"),
                                                 min_pairs = 20,
                                                 max_time_difference = 15,
                                                 correlation_floor = 0.4,
                                                 cor_method = "spearman") {
    get_feature_vector = function(feature_name) {
        as.numeric(cor_input[[feature_name]])
    }

    capture_fraction = as.numeric(capture_fraction)[1]
    if (!is.finite(capture_fraction) || capture_fraction <= 0 || capture_fraction >= 1) {
        stop("capture_fraction must be a numeric value between 0 and 1.")
    }
    if (!is.null(correlation_floor)) {
        correlation_floor = as.numeric(correlation_floor)[1]
        if (!is.finite(correlation_floor) || correlation_floor < 0 || correlation_floor >= 1) {
            stop("correlation_floor must be NULL or a numeric value in [0, 1).")
        }
    }

    required_anchor_cols = c("Mono_mass", "InChIKey", "Adduct_annotated", "time", "mz_time")
    required_annotated_cols = c("Mono_mass", "InChIKey", "Adduct_annotated", "time_annotated", "mz_time_annotated")
    if (
        !all(required_anchor_cols %in% colnames(primary_anchor_data)) ||
        !all(required_annotated_cols %in% colnames(annotated_adduct_data))
    ) {
        stop("Cannot estimate empirical adduct clustering thresholds because required calibration columns are missing.")
    }

    anchors = primary_anchor_data %>%
        dplyr::filter(Adduct_annotated %in% primary_adducts) %>%
        dplyr::transmute(
            Mono_mass = round(Mono_mass, 4),
            InChIKey,
            primary_adduct = Adduct_annotated,
            primary_time = time,
            primary_mz_time = mz_time
        ) %>%
        dplyr::filter(!is.na(InChIKey), !is.na(primary_time), !is.na(primary_mz_time)) %>%
        dplyr::distinct()

    if (nrow(anchors) == 0) {
        stop("Cannot estimate empirical adduct clustering thresholds because no clean primary-adduct anchors were found.")
    }

    secondary_candidates = annotated_adduct_data %>%
        dplyr::filter(!Adduct_annotated %in% primary_adducts) %>%
        dplyr::transmute(
            Mono_mass = round(Mono_mass, 4),
            InChIKey,
            secondary_adduct = Adduct_annotated,
            secondary_time = time_annotated,
            secondary_mz_time = mz_time_annotated
        ) %>%
        dplyr::filter(!is.na(InChIKey), !is.na(secondary_time), !is.na(secondary_mz_time)) %>%
        dplyr::distinct()

    if (nrow(secondary_candidates) == 0) {
        stop("Cannot estimate empirical adduct clustering thresholds because no same-metabolite secondary-adduct candidates were found.")
    }

    candidate_pairs = anchors %>%
        dplyr::inner_join(
            secondary_candidates,
            by = c("Mono_mass", "InChIKey"),
            relationship = "many-to-many"
        ) %>%
        dplyr::filter(primary_mz_time != secondary_mz_time) %>%
        dplyr::mutate(time_difference = abs(primary_time - secondary_time)) %>%
        dplyr::filter(
            is.finite(time_difference),
            time_difference <= max_time_difference,
            primary_mz_time %in% colnames(cor_input),
            secondary_mz_time %in% colnames(cor_input)
        )

    if (nrow(candidate_pairs) == 0) {
        stop("Cannot estimate empirical adduct clustering thresholds because no anchor-secondary pairs passed the retention-time and correlation-data filters.")
    }

    pair_correlations = vapply(seq_len(nrow(candidate_pairs)), function(i) {
        suppressWarnings(cor(
            get_feature_vector(candidate_pairs$primary_mz_time[i]),
            get_feature_vector(candidate_pairs$secondary_mz_time[i]),
            method = cor_method,
            use = "pairwise.complete.obs"
        ))
    }, numeric(1))

    calibration_pairs = candidate_pairs %>%
        dplyr::mutate(adduct_correlation = pair_correlations) %>%
        dplyr::filter(is.finite(adduct_correlation), adduct_correlation > 0) %>%
        dplyr::group_by(primary_mz_time, secondary_adduct) %>%
        dplyr::arrange(dplyr::desc(adduct_correlation), time_difference, .by_group = TRUE) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

    if (nrow(calibration_pairs) < min_pairs) {
        stop(
            "Cannot estimate empirical adduct clustering thresholds because only ",
            nrow(calibration_pairs),
            " calibration pairs were found; at least ",
            min_pairs,
            " are required."
        )
    }

    threshold_pairs = calibration_pairs
    threshold_method = "adduct_empirical"
    if (!is.null(correlation_floor)) {
        high_correlation_pairs = calibration_pairs %>%
            dplyr::filter(adduct_correlation >= correlation_floor)

        threshold_method = "adduct_empirical_hard_r_floor_insufficient_pairs"
        if (nrow(high_correlation_pairs) >= min_pairs) {
            threshold_pairs = high_correlation_pairs
            threshold_method = "adduct_empirical_hard_r_floor"
        }
    }

    time_threshold = as.numeric(stats::quantile(
        threshold_pairs$time_difference,
        probs = capture_fraction,
        na.rm = TRUE,
        names = FALSE,
        type = 7
    ))
    correlation_threshold = as.numeric(stats::quantile(
        threshold_pairs$adduct_correlation,
        probs = 1 - capture_fraction,
        na.rm = TRUE,
        names = FALSE,
        type = 7
    ))

    if (!is.finite(time_threshold)) {
        stop("Cannot estimate empirical adduct retention-time threshold from calibration pairs.")
    }
    if (!is.finite(correlation_threshold)) {
        stop("Cannot estimate empirical adduct correlation threshold from calibration pairs.")
    }

    time_threshold = max(1, min(max_time_difference, time_threshold))
    correlation_threshold = max(0, min(0.99, correlation_threshold))
    if (!is.null(correlation_floor)) {
        correlation_threshold = max(correlation_floor, correlation_threshold)
    }

    summary = data.frame(
        threshold_method = threshold_method,
        capture_fraction = capture_fraction,
        n_pairs = nrow(calibration_pairs),
        n_pairs_for_threshold = nrow(threshold_pairs),
        adduct_correlation_floor = ifelse(is.null(correlation_floor), NA_real_, correlation_floor),
        adduct_correlation_time_threshold = time_threshold,
        adduct_correlation_r_threshold = correlation_threshold,
        observed_time_capture = mean(calibration_pairs$time_difference <= time_threshold, na.rm = TRUE),
        observed_correlation_capture = mean(calibration_pairs$adduct_correlation >= correlation_threshold, na.rm = TRUE),
        median_time_difference = stats::median(calibration_pairs$time_difference, na.rm = TRUE),
        median_correlation = stats::median(calibration_pairs$adduct_correlation, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    list(
        adduct_correlation_time_threshold = time_threshold,
        adduct_correlation_r_threshold = correlation_threshold,
        pairs = calibration_pairs,
        summary = summary
    )
}

#' Estimate isotope clustering thresholds from high-confidence anchors
#'
#' Uses single-match primary-adduct features as anchors, pairs them with
#' same-metabolite isotope candidate features for the same adduct, and
#' chooses RT and correlation thresholds that retain the requested
#' fraction of empirical anchor-isotope pairs.
#'
#' @keywords internal
#' @noRd
estimate_isotope_clustering_thresholds = function(primary_anchor_data,
                                                  isotope_adduct_data,
                                                  cor_input,
                                                  capture_fraction = 0.8,
                                                  primary_adducts = c("M+H", "M-H"),
                                                  min_pairs = 20,
                                                  max_time_difference = 10,
                                                  correlation_floor = 0.4,
                                                  cor_method = "spearman") {
    get_feature_vector = function(feature_name) {
        as.numeric(cor_input[[feature_name]])
    }

    capture_fraction = as.numeric(capture_fraction)[1]
    if (!is.finite(capture_fraction) || capture_fraction <= 0 || capture_fraction >= 1) {
        stop("capture_fraction must be a numeric value between 0 and 1.")
    }
    if (!is.null(correlation_floor)) {
        correlation_floor = as.numeric(correlation_floor)[1]
        if (!is.finite(correlation_floor) || correlation_floor < 0 || correlation_floor >= 1) {
            stop("correlation_floor must be NULL or a numeric value in [0, 1).")
        }
    }

    required_anchor_cols = c("Mono_mass", "InChIKey", "Adduct_annotated", "time", "mz_time")
    required_isotope_cols = c("Mono_mass", "InChIKey", "Adduct_annotated", "time_annotated", "mz_time_annotated")
    if (
        !all(required_anchor_cols %in% colnames(primary_anchor_data)) ||
        !all(required_isotope_cols %in% colnames(isotope_adduct_data))
    ) {
        stop("Cannot estimate empirical isotope clustering thresholds because required calibration columns are missing.")
    }

    anchors = primary_anchor_data %>%
        dplyr::filter(Adduct_annotated %in% primary_adducts) %>%
        dplyr::transmute(
            Mono_mass = round(Mono_mass, 4),
            InChIKey,
            Adduct_annotated,
            primary_time = time,
            primary_mz_time = mz_time
        ) %>%
        dplyr::filter(!is.na(InChIKey), !is.na(primary_time), !is.na(primary_mz_time)) %>%
        dplyr::distinct()

    isotope_candidates = isotope_adduct_data %>%
        dplyr::filter(Adduct_annotated %in% primary_adducts) %>%
        dplyr::transmute(
            Mono_mass = round(Mono_mass, 4),
            InChIKey,
            Adduct_annotated,
            isotope_time = time_annotated,
            isotope_mz_time = mz_time_annotated
        ) %>%
        dplyr::filter(!is.na(InChIKey), !is.na(isotope_time), !is.na(isotope_mz_time)) %>%
        dplyr::distinct()

    if (nrow(anchors) == 0) {
        stop("Cannot estimate empirical isotope clustering thresholds because no clean primary-adduct anchors were found.")
    }
    if (nrow(isotope_candidates) == 0) {
        stop("Cannot estimate empirical isotope clustering thresholds because no same-metabolite isotope candidates were found.")
    }

    candidate_pairs = anchors %>%
        dplyr::inner_join(
            isotope_candidates,
            by = c("Mono_mass", "InChIKey", "Adduct_annotated"),
            relationship = "many-to-many"
        ) %>%
        dplyr::filter(primary_mz_time != isotope_mz_time) %>%
        dplyr::mutate(time_difference = abs(primary_time - isotope_time)) %>%
        dplyr::filter(
            is.finite(time_difference),
            time_difference <= max_time_difference,
            primary_mz_time %in% colnames(cor_input),
            isotope_mz_time %in% colnames(cor_input)
        )

    if (nrow(candidate_pairs) == 0) {
        stop("Cannot estimate empirical isotope clustering thresholds because no anchor-isotope pairs passed the retention-time and correlation-data filters.")
    }

    pair_correlations = vapply(seq_len(nrow(candidate_pairs)), function(i) {
        suppressWarnings(cor(
            get_feature_vector(candidate_pairs$primary_mz_time[i]),
            get_feature_vector(candidate_pairs$isotope_mz_time[i]),
            method = cor_method,
            use = "pairwise.complete.obs"
        ))
    }, numeric(1))

    isotope_ratios = vapply(seq_len(nrow(candidate_pairs)), function(i) {
        primary_intensity = get_feature_vector(candidate_pairs$primary_mz_time[i])
        isotope_intensity = get_feature_vector(candidate_pairs$isotope_mz_time[i])
        ratio = 2^(isotope_intensity - primary_intensity) * 100
        ratio[!is.finite(ratio) | ratio == 0] = NA_real_
        mean(ratio, na.rm = TRUE)
    }, numeric(1))

    calibration_pairs = candidate_pairs %>%
        dplyr::mutate(
            isotopic_correlation = pair_correlations,
            mean_isotope_to_primary_percent = isotope_ratios
        ) %>%
        dplyr::filter(
            is.finite(isotopic_correlation),
            isotopic_correlation > 0,
            is.finite(mean_isotope_to_primary_percent),
            mean_isotope_to_primary_percent <= 100
        ) %>%
        dplyr::group_by(primary_mz_time, Adduct_annotated) %>%
        dplyr::arrange(dplyr::desc(isotopic_correlation), time_difference, .by_group = TRUE) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

    if (nrow(calibration_pairs) < min_pairs) {
        stop(
            "Cannot estimate empirical isotope clustering thresholds because only ",
            nrow(calibration_pairs),
            " calibration pairs were found; at least ",
            min_pairs,
            " are required."
        )
    }

    threshold_pairs = calibration_pairs
    threshold_method = "isotope_empirical"
    if (!is.null(correlation_floor)) {
        high_correlation_pairs = calibration_pairs %>%
            dplyr::filter(isotopic_correlation >= correlation_floor)

        threshold_method = "isotope_empirical_hard_r_floor_insufficient_pairs"
        if (nrow(high_correlation_pairs) >= min_pairs) {
            threshold_pairs = high_correlation_pairs
            threshold_method = "isotope_empirical_hard_r_floor"
        }
    }

    time_threshold = as.numeric(stats::quantile(
        threshold_pairs$time_difference,
        probs = capture_fraction,
        na.rm = TRUE,
        names = FALSE,
        type = 7
    ))
    correlation_threshold = as.numeric(stats::quantile(
        threshold_pairs$isotopic_correlation,
        probs = 1 - capture_fraction,
        na.rm = TRUE,
        names = FALSE,
        type = 7
    ))

    if (!is.finite(time_threshold)) {
        stop("Cannot estimate empirical isotope retention-time threshold from calibration pairs.")
    }
    if (!is.finite(correlation_threshold)) {
        stop("Cannot estimate empirical isotope correlation threshold from calibration pairs.")
    }

    time_threshold = max(1, min(max_time_difference, time_threshold))
    correlation_threshold = max(0, min(0.99, correlation_threshold))
    if (!is.null(correlation_floor)) {
        correlation_threshold = max(correlation_floor, correlation_threshold)
    }

    summary = data.frame(
        threshold_method = threshold_method,
        capture_fraction = capture_fraction,
        n_pairs = nrow(calibration_pairs),
        n_pairs_for_threshold = nrow(threshold_pairs),
        isotopic_correlation_floor = ifelse(is.null(correlation_floor), NA_real_, correlation_floor),
        isotopic_correlation_time_threshold = time_threshold,
        isotopic_correlation_r_threshold = correlation_threshold,
        observed_time_capture = mean(calibration_pairs$time_difference <= time_threshold, na.rm = TRUE),
        observed_correlation_capture = mean(calibration_pairs$isotopic_correlation >= correlation_threshold, na.rm = TRUE),
        median_time_difference = stats::median(calibration_pairs$time_difference, na.rm = TRUE),
        median_correlation = stats::median(calibration_pairs$isotopic_correlation, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    list(
        isotopic_correlation_time_threshold = time_threshold,
        isotopic_correlation_r_threshold = correlation_threshold,
        pairs = calibration_pairs,
        summary = summary
    )
}
