#' Assign orientation
#'
#' Function to assign orientation based on visual inspection of residual plots. Assigns orientation based on residuals for multiple depth bins,
#' simultaneously.
#'
#' @param x Data frame variables used to assign orientation.
#' @param formula Function used to cast data frame, where the depth column used for casting is the response.
#' @param resid.col Name of the column containing residual values.
#' @param resid.cut Threshold value for assigning residuals.
#' @param output.name Name of the orientation column. Default = "orientation"

#' @return Returns `x` with an new column indicating the assigned orientation.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}

assign_orientation <- function(x, formula = vessel + cruise + haul ~ cdepth,
                               resid.col,
                               resid.cut,
                               output.name = "orientation") {

  x.cast <- reshape::cast(x, formula, value = resid.col)
  x.cast$orientation <- NA

  n.iter <- length(unique(x[,which(names(x) == labels(terms(formula)))]))
  last.col <- ncol(x.cast) - n.iter

  if(length(resid.cut) > 1) {
    resid.cut <- rev(resid.cut)
  } else {
    resid.cut <- rep(resid.cut, n.iter)
  }

  iter.count <- 1

  for(i in (ncol(x.cast)-1):last.col) {
    x.cast$orientation[x.cast[,i] > resid.cut[iter.count]] <- "Good"
    x.cast$orientation[x.cast[,i] < resid.cut[iter.count]] <- "Bad"
    iter.count <- iter.count + 1
  }

  x.out <- merge(x, unique(x.cast[, c(1:(ncol(x.cast)-n.iter-1), ncol(x.cast))]))

  names(x.out)[which(names(x.out) == "orientation")] <- output.name

  return(x.out)
}
