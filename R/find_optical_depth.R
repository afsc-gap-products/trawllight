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
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

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

  if(nrow(x) > 1) {

    x <- x[order(x$cdepth),]

    if(any(x$optical_depth > target.od)) {
      ind <- min(which(x$optical_depth > target.od))

      aa <- (x$optical_depth[ind] - x$optical_depth[ind-1])/(x$cdepth[ind] - x$cdepth[ind-1])
      bb <- x$optical_depth[ind] - x$cdepth[ind]*aa

      zz <- (target.od-bb)/aa
      x$output <- zz

  } else {
    warning("Target OD out of cast OD range.")
    x$output <- NA
  }
  } else {
    warning("Insufficient data.")
    x$output <- NA
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
