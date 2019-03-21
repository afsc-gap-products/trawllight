#' Estimate missing light measurement from reference depth bin
#'
#' Uses the linear attenuation coefficient of downwelling irradiance for the reference depth bin for cases with missing data.
#'
#' @param x Data frame containing depth, proportional light, binned light measurement, depth, and additional columns.
#' @param atten.col Name of the column containing the linear attenuation coefficient. Default = "k_linear"
#' @param ratio.col Name of the column containing proportional light level. Default = "light_ratio"
#' @param light.col Name of the column containing light measurement. Default = "trans_llight"
#' @param depth.col Name of the column containing depth. Default = "cdepth"
#' @param id.col Vector of column names which uniquely identify data for a single cast.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}

# Back-calculate light ratios to surface
estimate_reference_light <- function(x,
                             atten.col = "k_linear",
                             ratio.col = "light_ratio",
                             light.col = "trans_llight",
                             depth.col = "cdepth",
                             id.col = c("vessel", "cruise", "haul", "updown", "quality", "k_column")) {
  if(min(x$cdepth) != 1 & nrow(x) > 3) {
    names(x)[which(names(x) == ratio.col)] <- "light_ratio"
    names(x)[which(names(x) == atten.col)] <- "k_linear"
    names(x)[which(names(x) == depth.col)] <- "cdepth"
    names(x)[which(names(x) == light.col)] <- "trans_llight"

    do.not.duplicate <- c(id.col, "light_ratio", "k_linear", "cdepth")

    x.adjust <- x[which(x$cdepth == min(x$cdepth)),]

    # Mean linear attenuation coef. b/t the shallowest depth bin and second shallowest depth bin
    x.adjust$k_linear <- mean(c(x$k_linear[which(rank(x$cdepth) == 2)], x$k_linear[which(rank(x$cdepth) == 3)]))

    # Depth for the shallowest available depth bin
    rank1_depth <- x$cdepth[which(rank(x$cdepth) == 1)]

    # Determine how large the depth bin should be, calculate the depth where the middle of the reference depth bin should be.
    bin.adjust <- min(diff(x$cdepth[order(x$cdepth)]))/2

    # Calculate light ratio relative to depth for the middle of the bin
    x.adjust$light_ratio <- 1/exp(-1 * x.adjust$k_linear * (rank1_depth - bin.adjust))
    x.adjust[,which(!(names(x.adjust) %in% do.not.duplicate))] <- NA
    x.adjust$cdepth <- 1
    x.adjust$trans_llight <- x$trans_llight[which(rank(x$cdepth) == 1)] * x.adjust$light_ratio

    # Put it back
    x <- rbind(x, x.adjust)
    x$light_ratio <- x$light_ratio / max(x$light_ratio)

    names(x)[which(names(x) == "light_ratio")] <- ratio.col
    names(x)[which(names(x) == "k_linear")] <- atten.col
    names(x)[which(names(x) == "cdepth")] <- depth.col
    names(x)[which(names(x) == "trans_llight")] <- light.col

    x$estimate_ref <- T

  } else {
    x$estimate_ref <- F
  }
  return(x)

}
