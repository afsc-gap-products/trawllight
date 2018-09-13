#' Calculate optical depth
#'
#' \code{calculate_od} converts light measurements (photon flux density, etc.) to optical depth (Kirk 2011) for each observation, relative to the highest observed light level.
#'
#' @param x A vector of light measurements to be converted to optical depth.
#' @return Optical depth
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}
#'
#' @references Kirk, J.T.O. 2011. Light and photosynthesis in aquatic ecosystems. In 3rd edition. Cambridge University Press, New York.
#'
calculate_od <- function(x) {
  return(log(1)-log(x/max(x)))
}
