# Preserve the candidate columns when every mass group has no correlation
# evidence. Downstream joins still need an empty table with a stable schema.
ensure_msmica_cluster_result_schema = function(cluster_results, candidate_schema) {
    if (nrow(cluster_results) == 0L) {
        if (ncol(cluster_results) == 0L) {
            cluster_results = candidate_schema[0, , drop = FALSE]
        }
    }
    required_columns = list(
        correlation = NA_real_,
        adduct_corr_cluster = NA_integer_,
        mz_isotope = NA_real_,
        time_isotope = NA_real_,
        MSMICA_identification = NA_integer_,
        Probability = NA_real_
    )
    for (column in names(required_columns)) {
            if (!column %in% names(cluster_results)) {
                cluster_results[[column]] = rep(required_columns[[column]], nrow(cluster_results))
            }
    }
    cluster_results
}
