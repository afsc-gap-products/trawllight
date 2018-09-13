#' Calculate relative light level and linear attenuation coefficient
#'
#' \code{light_proportion} calculates the proportion of light and the linear attenuation coefficient for each depth bin, relative to the highest observed light level for a cast.
#'
#' @param x A data frame which contains light and depth measurements.
#' @param light.col Name of the column containing light measurements.
#' @param depth.col Name of the column containing depth measurements.
#' @return Returns a data frame containing input values, the proportion of light relative to the highest observed light measurement as \code{light_ratio}, and the attenuation coefficient of diffuse downwelling irradiance relative to the highest observed light measurement as \code{k_linear}.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}
#'

light_proportion <- function(x, light.col = "trans_llight", depth.col = "cdepth", ...) {

  names(x)[names(x) == light.col] <- "trans_llight"
  names(x)[names(x) == depth.col] <-  "cdepth"
  # Light ratio relative to shallowest bin.
  x$light_ratio <- x$trans_llight / max(x$trans_llight, na.rm = T)

  # Linear attenuation coefficient relative to shallowest bin.
  x$k_linear <- log(x$trans_llight/max(x$trans_llight)) / (min(x$cdepth) - x$cdepth)

  # Whole column attenuation coefficient.
  x$k_column <- rep((log(min(x$trans_llight)/max(x$trans_llight))) / (min(x$cdepth) - max(x$cdepth)), nrow(x))
  names(x)[names(x) == "trans_llight"] <- light.col
  names(x)[names(x) == "cdepth"] <-  depth.col
  return(x)
}
