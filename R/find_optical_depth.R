#' Depth of optical depth
#'
#' Function to estimate the depth of a user-specified optical depth.
#'
#' @param x Data frame containing optical depth and bin depth for a single cast.
#' @param od.col Name of the data frame column containing optical depth
#' @param depth.col Name of the data frame column containing depths
#' @param target.od Optical depth value of interest.
#' @param with.input Logical. Should the input data frame be returned with the optical depth? By default, false.
#' @param return.col Optional. Name of the column for output. Default = F.

#' @return Returns either a numerical vector with the optical depth for the cast, or a data frame with the input data and the optical depth data for the cast. Optical depths exceeding the maximum depth of the cast are returned as an NA.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}

find_optical_depth <- function(x,
                               od.col = "optical_depth",
                               depth.col = "cdepth",
                               target.od = 2.302585,
                               with.input = F,
                               return.col = F) {
  names(x)[which(names(x) == od.col)] <- "optical_depth"
  names(x)[which(names(x) == depth.col)] <- "cdepth"

  if(target.od <= 0) {
    stop("Invalid optical depth. Optical depth must be in the interval [0,Inf).")
  }

  if(nrow(x) > 4) {
    od_in <- x$optical_depth
    x.loess <- loess.as2(x = od_in, y = log(x$cdepth), criterion = "aicc")
    x$output<- exp(predict(x.loess, newdata = target.od))
  } else {
    warning("Insufficient data. Model fitting requires >4 observations.")
  }

  if(x$output[1] > max(x$cdepth) || is.na(x$output[1])) {
    warning("Estimated depth of water column depth greater than maximum depth bin. NA returned")
  }

  if(with.input == F) {
    return(x$output[1])
  } else {

  if(return.col != F) {
    names(x)[ncol(x)] <- return.col
  } else {
    names(x)[ncol(x)] <- "Z_optical_depth"
  }
  names(x)[which(names(x) == "optical_depth")] <- od.col
  names(x)[which(names(x) == "cdepth")] <- depth.col

  return(x)
  }
}
