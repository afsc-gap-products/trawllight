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
#' @return Returns a data frame containing input data, solar position, photosynthetically active radiation, and generalized additive model residuals.
#' @references Frouin, R., Lingner, D.W., Gautier, C., Baker, K.S., and Smith, R.C. 1989. A simple analytical formula to compute clear sky total and photosynthetically available solar irradiance at the ocean surface. J. Geophys. Res. 94(C7): 9731. doi:10.1029/JC094iC07p09731.
#' @references  Gary A. Nelson (2017). fishmethods: Fishery Science Methods and Models in R. R package version 1.10-4.https://CRAN.R-project.org/package=fishmethods
#' @references Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society (B) 73(1):3-36
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

tag_residuals_indirect <- function(x, formula = log10(trans_llight) ~ s(PAR, bs = "cr"),
                                   utc.offset = NULL,
                                   lat.col = "latitude",
                                   lon.col = "longitude",
                                   time.col = "start_time",
                                   light.col = "trans_llight",
                                   depth.bins = c(1, 3, 5, 7, 9), ...) {

  # Change column names to match processing
  names(x)[names(x) == lat.col] <- "latitude"
  names(x)[names(x) == lon.col] <- "longitude"
  names(x)[names(x) == time.col] <- "start_time"
  names(x)[names(x) == light.col] <- "trans_llight"
  
  if(class(x$start_time[1])[1] == "POSIXct") {
    # Use UTC times for fishmethods::astrocalc4r
    utc_start_time <- lubridate::with_tz(x$start_time, tzone = "UTC")
  }

  lout <- list()

  if(mean(c(depth.bins) %in% x$cdepth) < 1) {
    stop(paste0("tag_residuals_direct: Cannot calculate residuals. Some depth.bins not found in ", depth.col))
  }

  x <- cbind(x, fishmethods::astrocalc4r(day = lubridate::day(utc_start_time),
                                         month = lubridate::month(utc_start_time),
                                         year = lubridate::year(utc_start_time),
                                         hour = lubridate::hour(utc_start_time) + lubridate::minute(utc_start_time)/60,
                                         timezone = rep(0, nrow(x)),
                                         lat = x$latitude,
                                         lon = x$longitude,
                                         seaorland = "maritime"))

  for(i in 1:length(depth.bins)) {
    x_sub <- subset(x, cdepth == depth.bins[i])

    # GAM relating Frounin et al. (1989) model output to light in a depth bin.
    INDIRECT_GAM <- mgcv::gam(formula = formula, data = x_sub, ...)

    x_sub$indirect_residual <- residuals(INDIRECT_GAM)

    if(depth.bins[1] == depth.bins[i]) {
      output.df <- x_sub
    } else {
      output.df <- plyr::rbind.fill(output.df, x_sub)
    }

    lout[[i]] <- INDIRECT_GAM
    names(lout)[i] <- paste0("gam_dbin", depth.bins[i])

  }

  # Change column names to match input
  names(output.df)[names(output.df) == "latitude"] <- lat.col
  names(output.df)[names(output.df) == "longitude"] <-  lon.col
  names(output.df)[names(output.df) == "start_time"] <-  time.col
  names(output.df)[names(output.df) == "trans_llight"] <-  light.col

  lout$resid_df <- output.df
  return(lout)
}
