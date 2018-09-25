#' Indirect method for detecting tag obstruction
#'
#' \code{tag_residuals_indirect} fit a generalized additive model between surface light predictions from a global irradiance model (Frouin et al. 1989; fishmethods::astrocalc4r) and water column light, based on the position and time of the cast. Returns model residuals, estimate surface irradiance (photosynthetically active radiation), and solar position.
#'
#'
#' @param x Data frame containing light at depth, depth, latitude (decimal degrees), longitude (decimal degrees), cast start time in POSIXct format, and timezone offset relative to UTC and unique identifiers for survey vessels, years, etc.
#' @param formula an object of the class formula, passed to the generalized additive model. See documentation for \code{\link[mgcv]{gam}}.
#' @param utc.offset A numeric vector providing the timezone offset relative to UTC. See documentation for \code{\link[fishmethods]{astrocalc4r}}.
#' @param lat.col A character vector with the name of the column containing latitude.
#' @param lon.col A character vector with the name of the column containing longitude.
#' @param time.col A character vector with the name of the column containing start time.
#' @param light.col A character vector with the name of the column containing water column light measurements.
#' @param ... Additional values for gam.
#'
#' @return Returns a data frame containing input data, solar position, photosynthetically active radiation, and generalized additive model residuals.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}
#'
#' @references Frouin, R., Lingner, D.W., Gautier, C., Baker, K.S., and Smith, R.C. 1989. A simple analytical formula to compute clear sky total and photosynthetically available solar irradiance at the ocean surface. J. Geophys. Res. 94(C7): 9731. doi:10.1029/JC094iC07p09731.
#' @references  Gary A. Nelson (2017). fishmethods: Fishery Science Methods and Models in R. R package version 1.10-4.https://CRAN.R-project.org/package=fishmethods
#' @references Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society (B) 73(1):3-36

tag_residuals_indirect <- function(x, formula = log10(trans_llight) ~ s(PAR, bs = "cr"), utc.offset = -8, lat.col = "start_latitude", lon.col = "start_longitude", time.col = "start_time", light.col = "trans_llight", ...) {

  # Change column names to match processing
  names(x)[names(x) == lat.col] <- "start_latitude"
  names(x)[names(x) == lon.col] <- "start_longitude"
  names(x)[names(x) == time.col] <- "start_time"
  names(x)[names(x) == light.col] <- "trans_llight"

  x$hhour <- hour(x$start_time) + minute(x$start_time)/60

  x <- cbind(x, astrocalc4r(day = day(x$start_time), month = month(x$start_time), year = year(x$start_time), hour = hour(x$start_time) + minute(x$start_time)/60, timezone = rep(-8, nrow(x)), lat = x$start_latitude, lon = x$start_longitude, seaorland = "maritime"))

  # GAM relating Frounin et al. (1989) model output to light in a depth bin.
  LIGHT_GAM <- gam(formula = formula, data = x, ...)
  x$light_residual <- residuals(LIGHT_GAM)

  # Change column names to match input
  names(x)[names(x) == "start_latitude"] <- lat.col
  names(x)[names(x) == "start_longitude"] <-  lon.col
  names(x)[names(x) == "start_time"] <-  time.col
  names(x)[names(x) == "trans_llight"] <-  light.col

  return(x)
}
