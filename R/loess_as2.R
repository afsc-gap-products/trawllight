#' Local polynomial regression with automatic span selection and user-defined span bounds
#'
#' \code{loess.as2} fits a local polynomial gregression with automatic smoothing paramater selection, but has been modified to facilitate user-specified dynamic span bounds based on the number of observations using the \code{min.bins} argument. The function description and code have been modified from \code{fANCOVA::loess.as} v0.5-1 (Wang 2010).
#'
#' @param x as described for \code{fANCOVA::loess.as}
#' @param y as described for \code{fANCOVA::loess.as}
#' @param degree as described for \code{fANCOVA::loess.as}
#' @param criterion as described for \code{fANCOVA::loess.as}
#' @param family as described for \code{fANCOVA::loess.as}
#' @param user.span as described for \code{fANCOVA::loess.as}
#' @param plot as described for \code{fANCOVA::loess.as}
#' @param min.bins Minimum number of depth bins over which a span can smooth. Default = 3
#' @param ... as described for \code{fANCOVA::loess.as}
#'
#' @references Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. Journal of the American Statistical Association. 74, 829–836.
#' @references Golub, G., Heath, M. and Wahba, G. (1979). Generalized cross validation as a method for choosing a good ridge parameter. Technometrics. 21, 215–224.
#' @references Hurvich, C.M., Simonoff, J.S., and Tsai, C.L. (1998), Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. Journal of the Royal Statistical Society B. 60, 271–293.
#' @references Xiao-Feng Wang (2010). fANCOVA: Nonparametric Analysis of Covariance. R package version 0.5-1. https://CRAN.R-project.org/package=fANCOVA
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

loess.as2 <- function (x, y, degree = 1, criterion = c("aicc", "gcv"), family = c("gaussian", "symmetric"), min.bins = 3, user.span = NULL, plot = FALSE, ...)
{
  span.min <- min.bins/length(x) # Minimum observations
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  x <- as.matrix(x)
  if ((ncol(x) != 1) & (ncol(x) != 2))
    stop("The predictor 'x' should be one or two dimensional!!")
  if (!is.numeric(x))
    stop("argument 'x' must be numeric!")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric!")
  if (any(is.na(x)))
    stop("'x' contains missing values!")
  if (any(is.na(y)))
    stop("'y' contains missing values!")
  if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span)))
    stop("argument 'user.span' must be a numerical number!")
  if (nrow(x) != length(y))
    stop("'x' and 'y' have different lengths!")
  if (length(y) < 3)
    stop("not enough observations!")
  data.bind <- data.frame(x = x, y = y)
  if (ncol(x) == 1) {
    names(data.bind) <- c("x", "y")
  }
  else {
    names(data.bind) <- c("x1", "x2", "y")
  }
  opt.span <- function(model, criterion = c("aicc", "gcv"),
                       span.range = c(span.min, 1)) { # span range
    as.crit <- function(x) {
      span <- x$pars$span
      traceL <- x$trace.hat
      sigma2 <- sum(x$residuals^2)/(x$n - 1)
      aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n -
                                                          traceL - 2)
      gcv <- x$n * sigma2/(x$n - traceL)^2
      result <- list(span = span, aicc = aicc, gcv = gcv)
      return(result)
    }
    criterion <- match.arg(criterion)
    fn <- function(span) {
      mod <- update(model, span = span)
      as.crit(mod)[[criterion]]
    }
    result <- optimize(fn, span.range)
    return(list(span = result$minimum, criterion = result$objective))
  }
  if (ncol(x) == 1) {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x, degree = degree, family = family,
                    data = data.bind, ...)
      span1 <- opt.span(fit0, criterion = criterion)$span
    }
    else {
      span1 <- user.span
    }
    fit <- loess(y ~ x, degree = degree, span = span1, family = family,
                 data = data.bind, ...)
  }
  else {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x1 + x2, degree = degree, family = family,
                    data.bind, ...)
      span1 <- opt.span(fit0, criterion = criterion)$span
    }
    else {
      span1 <- user.span
    }
    fit <- loess(y ~ x1 + x2, degree = degree, span = span1,
                 family = family, data = data.bind, ...)
  }
  if (plot) {
    if (ncol(x) == 1) {
      m <- 100
      x.new <- seq(min(x), max(x), length.out = m)
      fit.new <- predict(fit, data.frame(x = x.new))
      plot(x, y, col = "lightgrey", xlab = "x", ylab = "m(x)",
           ...)
      lines(x.new, fit.new, lwd = 1.5, ...)
    }
    else {
      m <- 50
      x1 <- seq(min(data.bind$x1), max(data.bind$x1), len = m)
      x2 <- seq(min(data.bind$x2), max(data.bind$x2), len = m)
      x.new <- expand.grid(x1 = x1, x2 = x2)
      fit.new <- matrix(predict(fit, x.new), m, m)
      persp(x1, x2, fit.new, theta = 40, phi = 30, ticktype = "detailed",
            xlab = "x1", ylab = "x2", zlab = "y", col = "lightblue",
            expand = 0.6)
    }
  }
  return(fit)
}
