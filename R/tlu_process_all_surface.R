#' Loop surface_light
#'
#' Calculate surface light for upcasts and downcasts.
#'
#' @param dir.structure A character vector containing filepaths for all of the directories containing light measurements from surface/deck archival tags.
#' @param adjust.time Should trawllight::tlu_time_adjustments function be used to adjust surface times to match survey time.
#' @param ... Additional arguments to be passed to surface_light or time_adjustments
#' @export

tlu_process_all_surface <- function(dir.structure, adjust.time = T, ...) {

  surface.output <- NULL

  for(t in 1:length(dir.structure)) {

    # Check for CastTimes
    if(!file.exists(paste(dir.structure[t], "/CastTimes.csv", sep = ""))) {
      stop(paste("process_all_surface: CastTimes.csv not found in" , paste(dir.structure[t])))
    }

    # Import CastTImes
    print(paste("Processing", dir.structure[t]))
    cast.times <- read.csv(paste(dir.structure[t], "/CastTimes.csv", sep = ""))

    # Find names of deck files
    deck.files <- list.files(path = dir.structure[t], pattern = "^deck.*\\.csv", full.names = T)

    # Check for CastTimes
    if(length(deck.files) < 1) {
      warning(paste("process_all_surface: Deck light measurements not found in" , paste(dir.structure[t])))
    } else {

      
    
    #Import first deck file
    deck.data <- read.csv(file = deck.files[1], header = F, skip = 3)
    deck.data$ctime <- paste(deck.data[,1], deck.data[,2], sep = " ")

    # Import additional deck files if multiple exist in one directory
    if(length(deck.files) > 1) {
      for(b in 2:length(deck.files)) {
        deck.data <- rbind(deck.data, read.csv(file = deck.files[b], header = F, skip = 3))
      }
    }

    # Convert times into POSIXct
    deck.data$ctime <- as.POSIXct(strptime(deck.data$ctime, format = "%m/%d/%Y %H:%M:%S", tz = "America/Anchorage"))

    # Convert cast times to POSIXct format, add 30 second offset to each cast time to avoid truncating cast data
    cast.times$downcast_start <- as.POSIXct(strptime(cast.times$downcast_start,
                                                     format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
    cast.times$downcast_end <- as.POSIXct(strptime(cast.times$downcast_end,
                                                   format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
    cast.times$upcast_start <- as.POSIXct(strptime(cast.times$upcast_start,
                                                   format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
    cast.times$upcast_end <- as.POSIXct(strptime(cast.times$upcast_end,
                                                 format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))

    if(adjust.time) {
    # Correct cases where there is a mismatch between survey time and tag time
    deck.data <- trawllight::tlu_time_adjustments(light.data = deck.data,
                                   cast.data = cast.times)
    }
    if(nrow(deck.data) > 0) {

      # Find surface measurements
      surface_profiles <- trawllight::tlu_surface_light(light.data = deck.data,
                                        cast.data = cast.times,
                                        ...)
      if(is.null(surface.output)) {

        surface.output <- surface_profiles

      } else {
        if(!is.null(surface_profiles)) {
          surface.output <- rbind(surface.output, surface_profiles)
        }

      }
    }
    }
  }
  return(surface.output)
}
