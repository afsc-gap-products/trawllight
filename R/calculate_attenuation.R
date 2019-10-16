#' Instantaneous diffuse attenuation coefficient of downwelling irradiance
#'
#' \code{calculate_attenuation} fits a loess model between depth and light using an AICc-based span selection model adapated from \code{fANCOVA::loess.as}, then estimates the first derivative of the resultant model slope to approximate the instantaneous diffuse attenuation coefficient of downwelling irradiance.
#'
#' @param x Data frame containing depth and light for a single cast.
#' @param loess.criterion Criterion for choosing the most parsimonious model. Options are bias-corrected Akaike's Information Criterion ("aicc") or generalized cross-validation ("gcv").
#' @param loess.degree Degrees for loess model. Default = 1.
#' @param kz.binsize Depth interval for estimating instantaneous diffuse attenuation coefficient of downwelling irradiance. Default = 0.2.
#' @param min.range Minimum range of depths necessary for model fitting. Default = 10.
#' @param light.predict Logical indicating whether predicted values for light should be returned.
#' @param ... Additional arguments passed to loess fitting function
#' @return Returns a list containing three data frames: \code{attenuation} contains depth and fitted attenuation values, \code{loess.fit} contains the model summary statistics, and \code{fit_residuals} contains model fit residuals.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}
#'
#' @references Hurvich, C.M., Simonoff, J.S., and Tsai, C.-L. 1998. Smoothing parameter selection in nonparametric regression using an improved Akaike information criterion. J. R. Stat. Soc. B 60(2): 271-293.
#' @references Xiao-Feng Wang (2010). fANCOVA: Nonparametric Analysis of Covariance. R package version 0.5-1. https://CRAN.R-project.org/package=fANCOVA


calculate_attenuation <- function(x,
                                  light.col = "trans_llight",
                                  depth.col = "cdepth",
                                  loess.criterion = "aicc",
                                  loess.degree = 1,
                                  kz.binsize = 0.2,
                                  min.range = 10,
                                  light.predict = F,
                                  ...) {

  names(x)[which(names(x) == light.col)] <- "trans_llight"
  names(x)[which(names(x) == depth.col)] <- "cdepth"

  # Remove profiles with only a small portion of the water column sampled

  if((max(x$cdepth) - min(x$cdepth)) > min.range) { # Do not fit if depth range is < min.range
    if(length(unique(x$cdepth)) > 4) { # Cannot fit loess to fewer than three data points

    # Fit loess model
    N_depths <- seq(min(min(x$cdepth)), max(x$cdepth), kz.binsize)
    profile_light_loess <- loess.as2(x = x$cdepth, y = log(x$trans_llight), criterion = loess.criterion, degree = loess.degree, ...)

    # Adjust for overfitting caused by omitted data
    k <- 3
    while(profile_light_loess$s == Inf) {
      profile_light_loess <- loess.as2(x = x$cdepth, y = log(x$trans_llight), criterion = loess.criterion,
                                       degree = loess.degree,
                                       min.bins = k+1)
    }

    light_fit <- predict(profile_light_loess, newdata = N_depths)


    # Output data
    output <- data.frame(depth =  N_depths[1:(length(N_depths)-1)] + kz.binsize / 2,
                         k_aicc = diff(light_fit)/kz.binsize)

    loess.fit <- data.frame(span_fit = profile_light_loess$pars$span,
                            nobs = profile_light_loess$n,
                            enp = profile_light_loess$enp,
                            rse = profile_light_loess$s,
                            smooth_trace = profile_light_loess$trace.hat,
                            fit_method = loess.criterion)

    # Output residuals
    resids <- data.frame(residual = residuals(profile_light_loess),
                         log_trans_llight = log(x$trans_llight),
                         cdepth = x$cdepth)

    if(light.predict) {
      output$predict_light <- predict(profile_light_loess, newdata = output$depth)
      resids$predict_light <- light_fit
    }

    output_dfs <- list(fit_atten = loess.fit, attenuation = output, fit_residuals = resids)
    return(output_dfs)
    } else {
      warning("calculate_attenuation: Cannot fit loess. Cast has fewer than three Unique depth bins.")
      return(NULL)
    }
  } else {
    warning("calculate_attenuation: Did not fit loess. Total depth range of cast < min.range.")
    return(NULL)
  }
}
