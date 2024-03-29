% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tlutilities.R
\name{light_proportion}
\alias{light_proportion}
\title{Calculate relative light level and linear attenuation coefficient}
\usage{
light_proportion(x, light.col = "trans_llight", depth.col = "cdepth", ...)
}
\arguments{
\item{x}{A data frame which contains light and depth measurements.}

\item{light.col}{Name of the column containing light measurements.}

\item{depth.col}{Name of the column containing depth measurements.}
}
\value{
Returns a data frame containing input values, the proportion of light relative to the highest observed light measurement as \code{light_ratio}, and the attenuation coefficient of diffuse downwelling irradiance relative to the highest observed light measurement as \code{k_linear}.
}
\description{
\code{light_proportion} calculates the proportion of light and the linear attenuation coefficient for each depth bin, relative to the highest observed light level for a cast.
}
\keyword{internal}
