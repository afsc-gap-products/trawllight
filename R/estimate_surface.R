#' Estimate missing light measurements
#'
#' Uses the linear attenuation coefficient of downwelling irradiance for the shallowest depth bins with 'good' quality light measurements to estimate light for the shallowest depth bin.
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
estimate_surface <- function(x, atten.col = "k_linear", ratio.col = "light_ratio", light.col = "trans_llight", depth.col = "cdepth", id.col = c("vessel", "cruise", "haul", "updown", "quality", "k_column")) {
  if(min(x$cdepth) != 1 & nrow(x) > 3) {
    names(x)[which(names(x) == ratio.col)] <- "light_ratio"
    names(x)[which(names(x) == atten.col)] <- "k_linear"
    names(x)[which(names(x) == depth.col)] <- "cdepth"
    names(x)[which(names(x) == light.col)] <- "trans_llight"

    do.not.duplicate <- c(id.col, "light_ratio", "k_linear", "cdepth")

    x.light_ratio_adjust <- x[which(x$cdepth == min(x$cdepth)),]
    x.light_ratio_adjust$k_linear <- mean(x$k_linear[which(rank(x$cdepth) == 2)], x$k_linear[which(rank(x$cdepth) == 3)])
    rank1_depth <- x$cdepth[which(rank(x$cdepth) == 1)]

    bin.adjust <- min(diff(x$cdepth))/2

    # Calculate light ratio relative to depth for the middle of the bin
    x.light_ratio_adjust$light_ratio <- 1/exp(-1 * x.light_ratio_adjust$k_linear * (rank1_depth - bin.adjust))
    x.light_ratio_adjust[,which(!(names(x.light_ratio_adjust) %in% do.not.duplicate))] <- NA
    x.light_ratio_adjust$cdepth <- 1
    x.light_ratio_adjust$trans_llight <- x$trans_llight[which(rank(x$cdepth) == 1)] * x.light_ratio_adjust$light_ratio
    x <- rbind(x, x.light_ratio_adjust)
    x$light_ratio <- x$light_ratio / max(x$light_ratio)

    names(x)[which(names(x) == "light_ratio")] <- ratio.col
    names(x)[which(names(x) == "k_linear")] <- atten.col
    names(x)[which(names(x) == "cdepth")] <- depth.col
    names(x)[which(names(x) == "trans_llight")] <- light.col

  }
  return(x)

}
