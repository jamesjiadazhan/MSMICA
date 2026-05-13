#' QC_filter
#' 
#' This function is used to filter out features appear less than x percent of all samples, default is 20 percent.
#' @param x The metabolomics feature table
#' @param metabolite_start_column The column number where the intensity columns start in the feature table (typically, it is 3)
#' @param minimum_sample_appear The minimum percentage of samples that a feature should appear in to be kept in the dataset (default is 0.20), which means all features that appear in less than 20 percent of all samples will be removed.
#' @return A filtered feature table with features that appear in more than x percent of all samples
#' @examples
#' data("feature_table_exp_hilicpos")
#' QC_filter(x = feature_table_exp_hilicpos, metabolite_start_column = 3, minimum_sample_appear = 0.20)
#' @export QC_filter

QC_filter = function(x, metabolite_start_column, minimum_sample_appear = 0.20){
    sample_appear = ceiling((dim(x)[2] - (metabolite_start_column - 1)) * minimum_sample_appear)
    NumPres.All.Samples = rowSums(x != 0) - (metabolite_start_column - 1)
    x_clean = subset(x, NumPres.All.Samples > sample_appear)
    return(x_clean)
}