#' find.Overlapping.mzs
#' 
#' This function finds overlapping m/z and/or retention time values between two data sets.
#' @import data.table
#' @param dataA A feature table. The first column is considered as mz. If the second column is included (optional), it is considered as retention time. 
#' @param dataB A feature table. The first column is considered as mz. If the second column is included (optional), it is considered as retention time.
#' @param mz.thresh The m/z threshold for matching features between two data sets. Default is 5 ppm.
#' @param time.thresh The retention time threshold for matching features between two data sets. Default is NA (retention time is not used for matching). If this is used, the recommended input is 30 (seconds).
#' @return Matrix of overlapping features with columns: index.data.A: index of overlapping m/z in dataset A mz.data.A: m/z in dataset A time.data.A: retention time in dataset A index.data.B: index of overlapping m/z in dataset B mz.data.B: m/z in dataset B time.data.B: retention time in dataset B.
#' @export find.Overlapping.mzs

find.Overlapping.mzs <- function(dataA, dataB, mz.thresh = 5, time.thresh = NA) {
    library(data.table)
    # Convert inputs to data.table
    setDT(dataA)
    setDT(dataB)

    # --- Handle column names (simplified from original) ---
    # This logic can be adapted if the complex rules for apLCMS/XCMS are strictly needed.
    # For this example, we assume the first two columns are mz and time.
    setnames(dataA, old = names(dataA)[1], new = "mz.data.A")
    setnames(dataB, old = names(dataB)[1], new = "mz.data.B")
    
    if (!is.na(time.thresh)) {
        setnames(dataA, old = names(dataA)[2], new = "time.data.A")
        setnames(dataB, old = names(dataB)[2], new = "time.data.B")
    }
    
    # Add original index columns
    dataA[, index.A := .I]
    dataB[, index.B := .I]

    # --- Perform the match ---
    # 1. Calculate m/z tolerance for each feature in dataA
    dataA[, mz_tol := mz.thresh * mz.data.A / 1000000]
    dataA[, mz_lower := mz.data.A - mz_tol]
    dataA[, mz_upper := mz.data.A + mz_tol]

    # 2. Set keys for joining
    setkey(dataA, mz_lower, mz_upper)
    setkey(dataB, mz.data.B)

    # 3. Perform a non-equi join (overlap join) to find all B rows within the m/z range of A rows
    # This is the modern, fast replacement for the old foverlaps() method
    matches <- dataA[dataB, .(index.A, mz.data.A, time.data.A = if (!is.na(time.thresh)) x.time.data.A else NA,
                               index.B, mz.data.B, time.data.B = if (!is.na(time.thresh)) i.time.data.B else NA), 
                     on = .(mz_lower <= mz.data.B, mz_upper >= mz.data.B), nomatch = 0]
    
    # If no time threshold, the job is done
    if (is.na(time.thresh)) {
        # Select and reorder columns to match original output
        final_matches <- matches[, .(index.A, mz.data.A, index.B, mz.data.B)]
        return(as.data.frame(final_matches))
    }

    # 4. If using time threshold, calculate the difference and filter
    matches[, time.difference := abs(time.data.A - time.data.B)]
    final_matches <- matches[time.difference < time.thresh]

    # Select and reorder columns to match original output
    final_matches <- final_matches[, .(index.A, mz.data.A, time.data.A, index.B, mz.data.B, time.data.B, time.difference)]
    
    return(as.data.frame(final_matches))
}