#' Fit a robust monotonic GAM for retention time prediction (PredRet-style)
#'
#' Fits a monotonically constrained generalized additive model (GAM) that
#' maps a reference retention time (e.g. a public PredRet RT) onto the
#' observed retention time in a specific LC run. Robustness is obtained
#' in two passes: first an unconstrained thin-plate spline is fit to
#' flag and drop the worst-fitting \code{outlier_fraction} of anchors,
#' and then a penalized monotonic cubic-regression spline is fit to the
#' remaining anchors with sigmoid-weighted residuals so that leverage
#' points have bounded influence. The resulting object is consumed by
#' \code{predict_predret_gam()} to project metabolite RTs onto the study.
#'
#' @param anchors A data frame of RT anchor pairs with at least two
#'   numeric columns: the column named by \code{rt_col} (observed RT) and
#'   the column named by \code{rt_ref_col} (reference RT). Both should be
#'   on a common scale (e.g. normalized to 0-1 over the LC run).
#' @param rt_col Character. Name of the observed-RT column in
#'   \code{anchors}.
#' @param rt_ref_col Character. Name of the reference-RT column in
#'   \code{anchors}.
#' @param outlier_fraction Numeric in [0, 1]. Fraction of the anchor
#'   points with the largest residuals from the initial fit that are
#'   dropped before the final constrained fit. Defaults to 0.2.
#' @return A list carrying the fitted coefficients (\code{p}), the
#'   smoothing basis (\code{sm}), residuals of the final fit, the data
#'   used for the final fit, \code{r.sq} from the unconstrained fit, and
#'   diagnostic counts (\code{n_removed}, \code{removed_indices}).
#'   Returns \code{NULL} if fewer than 10 anchor points are provided.
#' @keywords internal
#' @noRd
fit_predret_gam = function(anchors, rt_col, rt_ref_col, outlier_fraction = 0.2) {
    # Extract vectors
    x = anchors[[rt_ref_col]] # Reference RT (0 to 1)
    y = anchors[[rt_col]]     # Observed RT (0 to 1)
    
    # Data checks
    if(length(x) < 10) {
        warning("Less than 10 points. PredRet might be unstable.")
        return(NULL)
    }
    
    # Add jitter if x is too discrete (prevents spline failures)
    if(length(unique(round(x, 2))) < 4) {
        x = jitter(x, amount = 0.01)
    }
    
    dat = data.frame(x = x, y = y)
    
    # --- Step A: Initial Unconstrained Fit to find Outliers ---
    k_val = min(length(unique(round(x, 2))), 10)
    f_init = gam(y ~ s(x, k = k_val, bs = "tp"), data = dat)
    
    # --- Step B: Remove the worst-fitting fraction of points ---
    abs_resid = abs(f_init$residuals)
    cutoff = quantile(abs_resid, probs = 1 - outlier_fraction)
    keep = abs_resid <= cutoff
    
    x_clean = x[keep]
    y_clean = y[keep]
    dat_clean = data.frame(x = x_clean, y = y_clean)
    
    # Recheck after removal
    if(length(x_clean) < 10) {
        warning("Fewer than 10 points remain after outlier removal. Using all points.")
        x_clean = x
        y_clean = y
        dat_clean = dat
        keep = rep(TRUE, length(x))
    }
    
    # --- Step C: Refit unconstrained GAM on cleaned data ---
    k_val_clean = min(length(unique(round(x_clean, 2))), 10)
    f_clean = gam(y ~ s(x, k = k_val_clean, bs = "tp"), data = dat_clean)
    
    # --- Step D: Calculate Robust Weights (Sigmoid) on cleaned data ---
    w = f_clean$residuals
    w = abs(w) / max(y_clean, na.rm = TRUE)
    w_final = pracma::sigmoid(w, a = -30, b = 0.1)
    
    # --- Step E: Monotonically Constrained Fit (PCLS) on cleaned data ---
    sm = smoothCon(s(x, k = k_val_clean, bs = "cr"), dat_clean, knots = NULL)[[1]]
    con = mono.con(sm$xp)
    
    G = list(X = sm$X, C = matrix(0,0,0), sp = f_clean$sp, p = sm$xp, y = y_clean, w = w_final)
    G$Ain = con$A; G$bin = con$b; G$S = sm$S; G$off = 0
    
    p = pcls(G)
    
    if(any(is.na(p))) {
        G$w = rep(1, length(w_final))
        p = pcls(G)
    }
    
    residuals_clean = as.numeric(y_clean - (sm$X %*% p))
    
    return(list(
        p = p, sm = sm, 
        residuals = residuals_clean,
        data = dat_clean, 
        r.sq = summary(f_clean)$r.sq,
        n_removed = sum(!keep),
        removed_indices = which(!keep)
    ))
}