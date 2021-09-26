#' Convert archival tag light to PFD
#'
#' \code{convert_light} converts light measurements from relative units recorded by a Wildlife Computers TDR-Mk9 archival tag equipped with a Hamatsu S2387 silicone photodiode to quantum units based on the equation of Kotwicki et al. (2009), or blue light intensity reported by Wildlife Computers.
#'
#' @param x Archival tag light measurement in relative units (~25-225)
#' @param convert.method Character vector indicating which method to use. Options: kotwicki, wc. Default: kotwicki
#' @references Kotwicki, S., Robertis, A., von Szalay, P., and Towler, R. 2009. The effect of light intensity on the availability of walleye pollock (Theragra chalcogramma) to bottom trawl and acoustic surveys. Can. J. Fish. Aquat. Sci. 66(6): 983–994. doi:10.1139/f09-055.
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export


convert_light <- function(x, convert.method = "kotwicki") {
  convert.method <- tolower(convert.method)
  if(convert.method == "kotwicki") {
    y <- 1*10^-8 * exp(0.1322*x)
  } else if(convert.method == "wc") {
    y <- 10^((x - 250)/20)
  } else {
    stop("Invalid method. Valid choices: kotwicki, wc")
  }
  return(y)
}