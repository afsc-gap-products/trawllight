### Automatic threshold selection
#'
#' \code{threshold_select} uses the first derivative of a kernel density function to estimate a threshold level at which to exclude casts with archival tag orientation errors.
#'
#' @param vals Numeric vector of values
#' @param method A character vector of length one indicating which automatic threshold detection method to use. Currently only accepts "kernel" for derivative of kernel density
#' @param bandwidth.select A vector of length one indicating which bandwidth select method to use. Currently available options are "UCV" (default), "OSCV, and "nrd0" See Details below.
#' @param mode.adjust A numeric vector of length one which adjusts the reference level for the mode for cases where the density function and derivative of the density function use different bandwidths. Default value generally should not be adjusted -0.1.
#' @param make.plots Logical vector indicating if diagnostic plots should be generated.
#' @return Returns a list which includes a threshold indicating the estimated level for incorrect archival tag orientation, bandwidth for the kernel density function (h.0), bandwidth for the derivative of the kernel density function (h.1), and the proportion of data excluded by the threshold (omit.rate)
#' @details bandwidth.select method can be one of unbiased cross-validation (UCV) implemented in the kedd package, one-sided cross validation (OSCV) implemented in the OSCV package, or the rule-of-thumb method implemented in stats::bandwidth ("nrd0").
#'
#' @references Jones, M. C., Marron, J. S. and Sheather, S. J. (1996). A brief survey of bandwidth selection for density estimation. Journal of the American Statistical Association, 91, 401–407.
#' @references Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' @references Savchuk, O.Y., Hart, J.D., and Sheather, S.J. (2010). Indirect cross-validation for density estimation. Journal of the American Statistical Association. 105, 415-423. doi: 10.1198/jasa.2010.tm08532
#' @references Scott, D.W. and George, R. T. (1987). Biased and unbiased cross-validation in density estimation. Journal of the American Statistical Association, 82, 1131–1146.
#' @references Silverman, B. W. (1986) Density Estimation. London: Chapman and Hall.
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export


threshold_select <- function(vals, method = "kernel",
                             bandwidth.select = "ICV",
                             mode.adjust = -0.1,
                             make.plots = T, silent = T) {
  if(method == "kernel") {

    if(bandwidth.select == "OSCV") {
      print("Using one-sided cross validation")
      h.0 <- OSCV::h_OSCV_dens(vals, stype = 0)
      h.1 <- h.0

    } else if(bandwidth.select == "UCV") {
      print("Using unbiased cross-validation")
      h.0 <- kedd::h.ucv(vals, deriv.order = 0)$h
      h.1 <- kedd::h.ucv(vals, deriv.order = 1)$h
      #h.1 <- h.0
    } else if(bandwidth.select == "nrd0") {
      print("Using default rule-of-thumb method")
      h.0 <- density(vals, bw = "nrd0")$bw
      h.1 <- h.0
    } else if(bandwidth.select == "ICV") {
      print("Using indirect cross-validation")
      h.0 <- ICV::h_ICV(vals)
      h.1 <- h.0
      } else {
      stop("Invalid bandwidth.select selected")
    }

    # Fit kernel density function to data
    f.0 <- kedd::dkde(vals, deriv.order = 0, h = h.0)
    f.1 <- kedd::dkde(vals, deriv.order = 1, h = h.0)

    f.0.max <- f.0$eval.points[which.max(f.0$est.fx)]

    # Find residual break point
    threshold <- mean(c(f.1$eval.points[max(which(f.1$est.fx < 0 & f.1$eval.points < (f.0.max + mode.adjust)))],
                        f.1$eval.points[max(which(f.1$est.fx < 0 & f.1$eval.points < (f.0.max + mode.adjust)))+1]))

    if(make.plots) {
      par(mfrow = c(1,2))
      print(plot(f.0, main = "Kernel density"))
      abline(h = 0, lty = 2, col = "black")
      abline(v = threshold, lty = 2, col = "red")
      print(plot(f.1, main = "First derivative of kernel density"))
      abline(h = 0, lty = 2, col = "black")
      abline(v = threshold, lty = 2, col = "red")
    }
  }
  output <- list(threshold = threshold, h.0 = h.0, h.1 = h.1, omit.rate = sum(vals[!is.na(vals)] < threshold)/length(vals[!is.na(vals)]))
  if(!silent) {
    print(output)
  }

  return(output)
}
