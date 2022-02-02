#' Wrapper function for filter_stepwise and calculate_attenuation using RACE data structure.
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. This function is designed to work with the file structure of RACE light data to subset all of the Mk9 data obtained during upscasts and downcasts. process_all runs functions trawllight::convert_light, trawllight::filter_stepwise, and trawllight::calculate_attenuation on the data.
#' 
#' @param dir.structure Vector of file paths to directories containing corr_MK9hauls.csv and CastTimes.csv files.
#' @param time.buffer Buffer around upcast_start and downcast_start
#' @param cast.dir Cast direction, either 'Downcast' or 'Upcast'
#' @param agg.fun Function to use to calculate light for a depth bin (Default = median)
#' @param binsize Bin width for depth bins
#' @param bin.gap Maximum gap in depth bins for a cast to still be considered 'good'
#' @param kz.binsize Bin size for interpolation
#' @export

tlu_process_all <- function(dir.structure,
                        time.buffer = 20,
                        cast.dir = "downcast",
                        agg.fun = median,
                        binsize = 2,
                        bin.gap = 6,
                        kz.binsize = 0.5,
                        silent = TRUE,
                        ...) {
  
  # Initialize
  loess_eval <- 1

  # Loops over directory structure ----
  for(i in 1:length(dir.structure)) {

    if(file.exists(paste(dir.structure[i], "/corr_MK9hauls.csv", sep = "")) &
       file.exists(paste(dir.structure[i], "/CastTimes.csv", sep = ""))) {

    corr_mk9hauls <- read.csv(paste(dir.structure[i], "/corr_MK9hauls.csv", sep = ""))
    casttimes <- read.csv(paste(dir.structure[i], "/CastTimes.csv", sep = ""))

    corr_mk9hauls$ctime <- as.POSIXct(strptime(corr_mk9hauls$ctime, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
    casttimes$downcast_start <- as.POSIXct(strptime(casttimes$downcast_start, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time.buffer
    casttimes$downcast_end <- as.POSIXct(strptime(casttimes$downcast_end, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time.buffer
    casttimes$upcast_start <- as.POSIXct(strptime(casttimes$upcast_start, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time.buffer
    casttimes$upcast_end <- as.POSIXct(strptime(casttimes$upcast_end, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time.buffer

    for(j in 1:nrow(casttimes)) {
      if(!silent) {
        print(paste("Cruise: ", casttimes$cruise[j], ", Vessel: ", casttimes$vessel[j], ", Haul: ", casttimes$haul[j], sep = ""))
      }
      vert <- trawllight::tlu_vertical_profiles(light.data = corr_mk9hauls,
                                cast.data = subset(casttimes, haul == casttimes$haul[j]))

      vert <- subset(vert, updown == cast.dir)

      if(nrow(vert) > 0) {
        vert$trans_llight <- trawllight::convert_light(vert$llight, ...)

        filtered <- trawllight::filter_stepwise(cast.data = vert,
                                    light.col = "trans_llight",
                                    depth.col = "cdepth",
                                    bin.size = binsize,
                                    bin.gap = bin.gap,
                                    agg.fun = agg.fun,
                                    ...)
        atten.out <- trawllight::calculate_attenuation(filtered, light.col = "trans_llight", depth.col = "cdepth", kz.binsize = kz.binsize)

        if(!is.null(atten.out)) {
          atten.out$attenuation$vessel <- vert$vessel[1]
          atten.out$attenuation$cruise <- vert$cruise[1]
          atten.out$attenuation$haul <- vert$haul[1]
          atten.out$attenuation$quality <- vert$quality[1]

          atten.out$fit_atten$vessel <- vert$vessel[1]
          atten.out$fit_atten$cruise <- vert$cruise[1]
          atten.out$fit_atten$haul <- vert$haul[1]

          atten.out$fit_residuals$vessel <- vert$vessel[1]
          atten.out$fit_residuals$cruise <- vert$cruise[1]
          atten.out$fit_residuals$haul <- vert$haul[1]

        }

        lr.out <- trawllight::light_proportion(filtered)

        if(cast.dir == "upcast") {
          lr.out$surface_time <- casttimes$upcast_end[j]
        } else if(cast.dir == "downcast") {
          lr.out$surface_time <- casttimes$downcast_start[j]
        }

        if(class(loess_eval) == "numeric") {
          loess_eval <- atten.out$fit_atten
          atten_values <- atten.out$attenuation
          resid_fit <- atten.out$fit_residuals
          light_ratios <- lr.out
        } else {
          loess_eval <- plyr::rbind.fill(loess_eval, atten.out$fit_atten)
          resid_fit <- plyr::rbind.fill(resid_fit, atten.out$fit_residuals)
          atten_values <- plyr::rbind.fill(atten_values, atten.out$attenuation)
          light_ratios <- plyr::rbind.fill(light_ratios, lr.out)
        }
      } else {
        if(silent == F) {
          print("No cast data found!")
        }
      }
    }
    } else {
      print(paste("File(s) not found in directory ", dir.structure, ". Directory skipped."))
    }
  }
  return(list(loess_eval = loess_eval,
              atten_values = atten_values,
              light_ratios = light_ratios,
              resid_fit = resid_fit))
}
