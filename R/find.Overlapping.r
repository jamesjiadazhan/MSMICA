#' find.Overlapping.mzs
#' 
#' This function finds overlapping m/z and/or retention time values between two data sets. Originally written by Dr. Karan Uppal, Emory University.
#' @param dataA A feature table, including the one created by apLCMS or XCMS. This needs to be consistent with the alignment.tool input. If alignment.tool is NA, the first column is considered as mz. If the second column is included (optional), it is considered as retention time. 
#' @param dataB A feature table, including the one created by apLCMS or XCMS. This needs to be consistent with the alignment.tool input. If alignment.tool is NA, the first column is considered as mz. If the second column is included (optional), it is considered as retention time.
#' @param mz.thresh The m/z threshold for matching features between two data sets. Default is 5 ppm.
#' @param time.thresh The retention time threshold for matching features between two data sets. Default is NA (retention time is not used for matching). If this is used, the recommended input is 30 (seconds).
#' @param alignment.tool The alignment tool used to create the feature table, inclding "apLCMS", "XCMS", or "NA". Default is NA. Use "NA" if the input matrix includes only m/z or both m/z and retnetion time values.
#' @return Matrix of overlapping features with columns: index.data.A: index of overlapping m/z in dataset A mz.data.A: m/z in dataset A time.data.A: retention time in dataset A index.data.B: index of overlapping m/z in dataset B mz.data.B: m/z in dataset B time.data.B: retention time in dataset B.
#' @export find.Overlapping.mzs

find.Overlapping.mzs <- function(dataA, dataB, mz.thresh = 5, time.thresh = NA, alignment.tool = NA) {
    # Convert input data to data frames
    data_a <- as.data.frame(dataA)
    data_b <- as.data.frame(dataB)

    # Remove unnecessary variables
    rm(dataA)
    rm(dataB)

    # Initialize variables
    com_mz_num <- 1

    # Get column names
    col.names.dataA <- colnames(data_a)
    col.names.dataB <- colnames(data_b)

    if (!is.na(alignment.tool)) {
        if (alignment.tool == "apLCMS") {
            sample.col.start <- 5
        } else if (alignment.tool == "XCMS") {
            sample.col.start <- 9
            col.names.dataA[1] <- "mz"
            col.names.dataA[2] <- "time"
            col.names.dataB[1] <- "mz"
            col.names.dataB[2] <- "time"
            colnames(data_a) <- col.names.dataA
            colnames(data_b) <- col.names.dataB
        }
    } else {
        col.names.dataA[1] <- "mz"
        col.names.dataB[1] <- "mz"
        if (!is.na(time.thresh)) {
            col.names.dataA[2] <- "time"
            col.names.dataB[2] <- "time"
            print("Using the 1st columns as \"mz\" and 2nd columns as \"retention time\"")
        } else {
            print("Using the 1st columns as \"mz\"")
        }
        colnames(data_a) <- col.names.dataA
        colnames(data_b) <- col.names.dataB
    }

    # Create header for the matrix with common features
    if (!is.na(time.thresh)) {
        mznames <- c("index.A", "mz.data.A", "time.data.A", "index.B", "mz.data.B", "time.data.B", "time.difference")
    } else {
        mznames <- c("index.A", "mz.data.A", "index.B", "mz.data.B")
    }

    # Step 1: Group features by m/z
    mz_groups <- lapply(1:dim(data_a)[1], function(j) {
        commat <- {}

        ppmb <- (mz.thresh) * (data_a$mz[j] / 1000000)
        getbind_same <- which(abs(data_b$mz - data_a$mz[j]) <= ppmb)

        # If time threshold is specified
        if (!is.na(time.thresh) && length(getbind_same) > 0) {
            all_matches <- list()

            for (comindex in 1:length(getbind_same)) {
                tempA <- cbind("index.A" = j, "mz.data.A" = data_a[j, 1], "time.data.A" = data_a[j, 2])
                tempB <- cbind("index.B" = getbind_same[comindex], "mz.data.B" = data_b[getbind_same[comindex], 1], "time.data.B" = data_b[getbind_same[comindex], 2])
                
                timediff <- abs(data_a[j, 2] - data_b[getbind_same[comindex], 2])
                temp <- cbind(tempA, tempB, "time.difference" = timediff)

                if (timediff < time.thresh) {
                    all_matches <- append(all_matches, list(as.data.frame(temp)))
                }
            }

            if (length(all_matches) > 0) {
                commat <- do.call("rbind", all_matches)
                rownames(commat) <- paste("mz", j, "_", seq(1, length(all_matches)), sep = "")
            }

        } else if (length(getbind_same) > 0) {
            for (comindex in 1:length(getbind_same)) {
                tempA <- cbind("index.A" = j, "mz.data.A" = data_a[j, 1])
                tempB <- cbind("index.B" = getbind_same[comindex], "mz.data.B" = data_b[getbind_same[comindex], 1])
                temp <- cbind(tempA, tempB)
                commat <- rbind(commat, temp)
            }
            rownames(commat) <- paste("mz", j, "_", seq(1, length(getbind_same)), sep = "")
        }

        return(as.data.frame(commat))
    })

    # Combine all mz_groups into a final data frame
    commat <- do.call("rbind", mz_groups)

    return(commat)
}