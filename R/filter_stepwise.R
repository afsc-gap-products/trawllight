#' Bin and filter light measurements
#'
#' \code{filter_stepwise} aggregates light measurements from a single cast into depth bins using a specificed function (e.g. 'median'), removes light measurements using a stepwise algorithm, and assigns data continuity grades based on a threshold criteria.
#'
#' @param cast.data Data frame containing light measurements and depth.
#' @param bin.size The size of the depth bin used for aggregation. Default = 2.
#' @param bin.gap The maximum size of data gap before a profile is considered to not meet continuity standards. Units are in units of depth, not the number of bins. Default = 6.
#' @param agg.fun Function used to aggregate light measurements for each depth bin.
#' @param ... Additional arguments passed to findInterval function for binning light measurements
#' @return Data frame with light by depth bin and continuity grade (1 = Good, -999 = Bad)

filter_stepwise <- function(cast.data,
                            light.col,
                            depth.col,
                            id.cols = c("vessel", "cruise", "haul", "updown"),
                            bin.size = 2,
                            bin.gap = 6,
                            agg.fun, ...) {

  names(cast.data)[which(names(cast.data) == light.col)] <- "trans_llight"
  names(cast.data)[which(names(cast.data) == depth.col)] <- "cdepth"

  max.depth <- max(ceiling(cast.data$cdepth), na.rm = T)

  # Bin by depth with bins centered
  cast.data$cdepth <- findInterval(cast.data$cdepth, seq(0, max.depth, bin.size), rightmost.closed = T, left.open = F, ...) * bin.size - bin.size/2

  # Calculate binned light level using user-specified function
  light_at_depth <- aggregate(formula = as.formula(paste("trans_llight", paste(c(id.cols, "cdepth"), collapse = "+"), sep = "~")),
                              data = cast.data,
                              FUN = agg.fun)


  light_at_depth <- light_at_depth[order(light_at_depth$cdepth),]

  # Stepwise measurement removal loop
  p2 <- 1
  while(p2 < nrow(light_at_depth) ) {
    if(nrow(light_at_depth) >= (p2 + 1)) {
      if((light_at_depth$trans_llight[p2 + 1] > light_at_depth$trans_llight[p2])) {
        light_at_depth <- light_at_depth[-p2,]
        p2 <- 0 # Index back to start
      }
    }
    p2 <- p2 + 1
  }

  # Assign data continuity codes. -999 indicates gap >= bin.gap
  if(max(diff(light_at_depth$cdepth)) <= bin.gap & min(light_at_depth$cdepth + 1) <= bin.gap) {
    light_at_depth$quality <- 1
  } else {
    light_at_depth$quality <- -999
  }

  names(light_at_depth)[which(names(light_at_depth) == "trans_llight")] <- light.col
  names(light_at_depth)[which(names(light_at_depth) == "cdepth")] <- depth.col

  return(light_at_depth)
}
