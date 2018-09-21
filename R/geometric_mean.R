#' Geometric mean
#'
#' Calculates the geometric mean of a numeric vector.
#'
#' @param x A numeric vector
#' @param na.rm Logical. Should NA values be omitted from calculation? Default = T

geometric.mean <- function(x, na.rm = T) {

  if(any(x <= 0)) {
    stop("geometric.mean: Cannot calculate a geometric mean. Input vector contains negative values.")
  }

  x <- x[!is.na(x)]

  return(exp(mean(log(x))))
}
