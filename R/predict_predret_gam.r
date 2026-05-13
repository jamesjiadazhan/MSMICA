#' Predict retention times from a fitted PredRet-style GAM
#'
#' Uses the smoothing basis (\code{sm}) and penalized coefficients
#' (\code{p}) produced by \code{fit_predret_gam()} to project new
#' reference retention times onto the observed RT scale of the current
#' LC run. Negative predictions (which cannot be physically valid) are
#' clamped to 0. If the model object is \code{NULL} (as returned when
#' too few anchor points are available) the function returns a vector of
#' \code{NA}s of the same length as the input so that downstream code
#' can skip RT prediction gracefully.
#'
#' @param new_rt_ref Numeric vector of reference retention times to
#'   project.
#' @param model_obj List returned by \code{fit_predret_gam()}, or
#'   \code{NULL}.
#' @return Numeric vector of predicted retention times on the observed
#'   RT scale, with negative values clamped to 0.
#' @keywords internal
#' @noRd
predict_predret_gam = function(new_rt_ref, model_obj) {
    if(is.null(model_obj)) return(rep(NA, length(new_rt_ref)))
    
    # Generate prediction matrix for new data
    # Note: Predict.matrix requires a dataframe with column name 'x' matching training
    fv = Predict.matrix(model_obj$sm, data.frame(x = new_rt_ref)) %*% model_obj$p
    
    y_pred = as.numeric(fv)
    y_pred[y_pred < 0] = 0 # Clamp negative RTs to 0
    return(y_pred)
}