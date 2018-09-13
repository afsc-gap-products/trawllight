#' Convert archival tag light to PFD
#'
#' \code{convert_light} converts light measurements from relative units recorded by a Wildlife Computers TDR-Mk9 archival tag to photon flux density, based on the equation of Kotwicki et al. (2009).
#'
#' @param x Archival tag light measurement in relative units (0-256)
#' @return Light in units of photon flux density.
#'
#' @references Kotwicki, S., Robertis, A., von Szalay, P., and Towler, R. 2009. The effect of light intensity on the availability of walleye pollock (Theragra chalcogramma) to bottom trawl and acoustic surveys. Can. J. Fish. Aquat. Sci. 66(6): 983–994. doi:10.1139/f09-055.

convert_light <- function(x) {
  y <- 1*10^-8 * exp(0.1322*x)
  return(y)
}
