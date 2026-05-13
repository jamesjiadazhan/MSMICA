#' Spearman correlation between candidate precursor and product features
#'
#' Given a data frame of candidate precursor-product pairs (one pair per
#' row, each encoded as a \code{mz_time} feature identifier), this
#' helper looks up each feature's intensity vector in the globally
#' available \code{MSMICA_cor_input} matrix and computes a Spearman
#' \code{cor.test()} per pair. Pairs that appear twice in reversed
#' order (\code{(A,B)} and \code{(B,A)}) are de-duplicated by imposing
#' a canonical ordering before \code{distinct()}. Used by the Bayesian
#' scoring step to bring biochemical reaction network evidence into the
#' identification.
#'
#' @param precursor_product_correlation_input_data Data frame of
#'   candidate pairs containing the two feature-id columns named by
#'   \code{precursor_col} and \code{product_col}.
#' @param precursor_col Character. Name of the precursor feature column
#'   in the input data.
#' @param product_col Character. Name of the product feature column in
#'   the input data.
#' @return A data frame with columns \code{connection_1},
#'   \code{connection_2}, \code{correlation}, and \code{p_value}, one
#'   row per unique unordered feature pair.
#' @keywords internal
#' @noRd
precursor_product_correlation = function(precursor_product_correlation_input_data, precursor_col = NULL, product_col = NULL,
                                         cor_input) {
    # Require callers to tell us which columns encode the precursor vs. product
    if (is.null(precursor_col) || is.null(product_col)) {
        stop("precursor_col and product_col must be specified")
    }
    feature_identified = precursor_product_correlation_input_data[[precursor_col]]
    feature_annotated = precursor_product_correlation_input_data[[product_col]]
    
    # Preallocate vectors for results
    n = length(feature_identified)
    cor_values = numeric(n)
    p_values = numeric(n)
    
    # Compute one Spearman correlation per pair. MSMICA_cor_input is indexed
    # by feature id (mz_time), so we look the vectors up by column name.
    for (i in seq_len(n)) {
        connection_1_vector = cor_input[[feature_identified[i]]]
        connection_2_vector = cor_input[[feature_annotated[i]]]

        cor_test_result = cor.test(
            connection_1_vector,
            connection_2_vector,
            method = "spearman",
            alternative = "two.sided",
            exact = FALSE
        )
        cor_values[i] = cor_test_result$estimate
        p_values[i] = cor_test_result$p.value
    }
    

    # Create the final data frame
    cor_MSMICA_data = data.frame(
        connection_1 = feature_identified,
        connection_2 = feature_annotated,
        correlation = cor_values,
        p_value = p_values
    )

    # An unordered (A,B) pair and its reverse (B,A) should be counted once. We
    # build a canonical (min, max) key so distinct() collapses the duplicates.
    cor_MSMICA_data = cor_MSMICA_data %>%
        mutate(
            conn_min = pmin(connection_1, connection_2),
            conn_max = pmax(connection_1, connection_2)
        ) %>%
        distinct(conn_min, conn_max, .keep_all = TRUE) %>%
        # Drop helper columns and rename back
        select(-connection_1, -connection_2) %>%
        rename(
            connection_1 = conn_min,
            connection_2 = conn_max
        )

    return(cor_MSMICA_data)
}