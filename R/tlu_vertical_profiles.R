#' Subsets light measurements from upcasts/downcasts
#'
#' Assigns light measurements to upcast or downcast based on downcast start time and downcast end time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @export

tlu_vertical_profiles <- function(light.data, cast.data) {

  # Remove surface data
  light.data <- subset(light.data, cdepth >= 0)

  # Create empty columns for cast direction (updown) and haul number
  light.data$updown <- rep(NA, nrow(light.data))
  light.data$haul <- rep(NA, nrow(light.data))

  # Assign upcast or downcast to times
  light.data$updown[light.data$ctime > cast.data$downcast_start & light.data$ctime < cast.data$downcast_end[1]] <- "downcast"
  light.data$updown[light.data$ctime > cast.data$upcast_start & light.data$ctime < cast.data$upcast_end[1]] <- "upcast"
  light.data$haul[light.data$ctime > cast.data$downcast_start & light.data$ctime < cast.data$upcast_end[1]] <- cast.data$haul[1]

  # Remove on-bottom and errant sampling not from casts
  light.data <- subset(light.data, is.na(updown) == F)
  return(light.data)
}
