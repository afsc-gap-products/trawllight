#' Wrapper function for filter_stepwise and calculate_attenuation using RACE data structure.
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. This function is designed to work with the file structure of RACE light data to subset all of the Mk9 data obtained during upcasts and downcasts. process_all runs functions trawllight::convert_light, trawllight::filter_stepwise, and trawllight::calculate_attenuation on the data.
#' 
#' @param dir.path Vector of file paths to directories containing corr_MK9hauls.csv and CastTimes.csv files.
#' @param time.buffer Buffer around upcast_start and downcast_start
#' @param cast.dir Cast direction, either 'Downcast' or 'Upcast'
#' @param agg.fun Function to use to calculate light for a depth bin (Default = median)
#' @param binsize Bin width for depth bins
#' @param bin.gap Maximum gap in depth bins for a cast to still be considered 'good'
#' @param kz.binsize Bin size for interpolation
#' @noRd

tlu_process_all <- function(dir.path,
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
    if(file.exists(paste(dir.path, "/corr_MK9hauls.csv", sep = "")) &
       file.exists(paste(dir.path, "/CastTimes.csv", sep = ""))) {
      
      corr_mk9hauls <- read.csv(paste(dir.path, "/corr_MK9hauls.csv", sep = ""))
      casttimes <- read.csv(paste(dir.path, "/CastTimes.csv", sep = ""))
      
      ctime_1 <- corr_mk9hauls$ctime[1]
      ctime_2 <- casttimes$downcast_start[1]
      corr_mk9hauls$ctime <- as.POSIXct(strptime(corr_mk9hauls$ctime, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
      casttimes$downcast_start <- as.POSIXct(strptime(casttimes$downcast_start, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time.buffer
      casttimes$downcast_end <- as.POSIXct(strptime(casttimes$downcast_end, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time.buffer
      casttimes$upcast_start <- as.POSIXct(strptime(casttimes$upcast_start, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time.buffer
      casttimes$upcast_end <- as.POSIXct(strptime(casttimes$upcast_end, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time.buffer
      
      if(is.na(casttimes$downcast_start[1])) {
        stop(paste(dir.path, "/CastTimes.csv date format not recognized! (", ctime_2, " >>> ", casttimes$downcast_start[1], ")" , sep = ""))
      }
      
      if(is.na(corr_mk9hauls$ctime[1])) {
        stop(paste(dir.path, "/corr_MK9hauls.csv date format not recognized! (", ctime_1, " >>> ", corr_mk9hauls$ctime[1], ")" , sep = ""))
      }
      
      for(j in 1:nrow(casttimes)) {
        if(!silent) {
          print(paste("Cruise: ", casttimes$cruise[j], ", Vessel: ", casttimes$vessel[j], ", Haul: ", casttimes$haul[j], sep = ""))
        }
        vert <- trawllight:::tlu_vertical_profiles(light.data = corr_mk9hauls,
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
          
          lr.out <- trawllight:::light_proportion(filtered)
          
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
      return(list(loess_eval = loess_eval,
                  atten_values = atten_values,
                  light_ratios = light_ratios,
                  resid_fit = resid_fit))
    } else {
     print(paste("corr_MK9hauls.csv and/or CastTimes.csv file(s) not found in ", dir.path, "."))
      return(NULL)
    }

}

#' Loop surface_light
#'
#' Calculate surface light for upcasts and downcasts.
#'
#' @param dir.structure A character vector containing filepaths for all of the directories containing light measurements from surface/deck archival tags.
#' @param adjust.time Should trawllight::tlu_time_adjustments function be used to adjust surface times to match survey time.
#' @param survey RACE survey region as a character vector ("BS", "NBS", "AI", "GOA" or "SLOPE")
#' @param ... Additional arguments to be passed to surface_light or time_adjustments
#' @noRd

tlu_process_all_surface <- function(dir.structure, adjust.time = T, survey, ...) {

  surface.output <- NULL
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  
  for(t in 1:length(dir.structure)) {
    
    # Check for CastTimes
    if(!file.exists(paste(dir.structure[t], "/CastTimes.csv", sep = ""))) {
      print(paste("process_all_surface: CastTimes.csv not found in" , paste(dir.structure[t])))
    } else {
      
      # Import CastTImes
      print(paste("Processing", dir.structure[t]))
      cast.times <- read.csv(paste(dir.structure[t], "/CastTimes.csv", sep = ""))
      
      deck.exists <- file.exists(paste0(dir.structure[t], "/corr_deck.csv"))
      
      if(deck.exists) {
        print(paste0("tlu_process_all_surface: Loading ", dir.structure[t], "/corr_deck.csv"))
        deck.data <- read.csv(file = paste0(dir.structure[t], "/corr_deck.csv"))
        deck.data$ctime <- as.POSIXct(strptime(deck.data$ctime, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
      } else {
        # Find names of deck files
        deck.files <- list.files(path = dir.structure[t], pattern = "^deck.*\\.csv", full.names = T)
        
        # Check for CastTimes
        if(length(deck.files) < 1) {
          warning(paste("process_all_surface: Deck light measurements not found in" , paste(dir.structure[t])))
          next
        } else {
          
          #Import first deck file
          deck.data <- read.csv(file = deck.files[1], header = F, skip = 3)
          deck.data$ctime <- paste(deck.data[,1], deck.data[,2], sep = " ")
          
          # Import additional deck files if multiple exist in one directory
          if(length(deck.files) > 1) {
            for(b in 2:length(deck.files)) {
              add.deck <- read.csv(file = deck.files[b], header = F, skip = 3)
              add.deck$ctime <- paste(add.deck[,1], add.deck[,2], sep = " ")
              deck.data <- rbind(deck.data, add.deck)
            }
          }
          # Convert times into POSIXct
          deck.data$ctime <- as.POSIXct(strptime(deck.data$ctime, format = "%m/%d/%Y %H:%M:%S", tz = "America/Anchorage"))
          
          if(is.na(deck.data$ctime[1])) {
            stop(paste(dir.structure[t], "date format not recognized! (", deck.data$ctime, " >>> ", deck.data$ctime[1], ")" , sep = ""))
          }
          
        }
      }
          
          # Convert cast times to POSIXct format, add 30 second offset to each cast time to avoid truncating cast data
          cast.times$downcast_start <- as.POSIXct(strptime(cast.times$downcast_start,
                                                           format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
          cast.times$downcast_end <- as.POSIXct(strptime(cast.times$downcast_end,
                                                         format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
          cast.times$upcast_start <- as.POSIXct(strptime(cast.times$upcast_start,
                                                         format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
          cast.times$upcast_end <- as.POSIXct(strptime(cast.times$upcast_end,
                                                       format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
          
          if(adjust.time & !deck.exists) {
            # Correct cases where there is a mismatch between survey time and tag time
            deck.data <- trawllight:::tlu_time_adjustments(light.data = deck.data,
                                                           cast.data = cast.times,
                                                           survey = survey)
          }
          if(nrow(deck.data) > 0) {
            
            # Find surface measurements
            surface_profiles <- trawllight:::tlu_surface_light(light.data = deck.data,
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

#' Get haul data and write to RDS
#' 
#' Run query to retrieve haul data from RACEBASE
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @noRd

tlu_prep_haul_data <- function(channel = NULL, 
                               survey) {
  
  survey <- toupper(survey)
  
  channel <- trawllight:::get_connected(channel = channel)
  
  qry <- paste0("select vessel, cruise, haul, start_time, stationid, start_latitude, start_longitude, end_latitude, end_longitude, bottom_depth, performance, haul_type, region, stratum from racebase.haul where cruise > 200400 and region = '", survey, "'")
  
  haul_dat <- RODBC::sqlQuery(channel = channel,
                              query = qry)
  names(haul_dat) <- tolower(names(haul_dat))
  
  saveRDS(haul_dat, 
          file = here::here("data", paste0(survey, "_haul_data.rds")))
}

#' Write light data directories to CSV
#' 
#' Find all of the directories with light data for a region and save the list of directories to a csv file.
#' 
#' @param survey RACE survey code as a character vector.
#' @param light_data_root Root directory for light data as a character vector.
#' @param omit_string Optional character vector to filter directories to be removed.
#' @noRd

tlu_prep_dir_list <- function(survey, 
                              light_data_root = "G:/RACE_LIGHT/LightData/Data",
                              omit_string = "oldtags") {
  
  region_light <- c("ebs", "nbs", "goa", "ai")[match(survey, c("BS", "NBS", "GOA", "AI"))]
  
  region_dirs <- list.files(light_data_root, 
                            pattern = region_light, 
                            recursive = TRUE, 
                            include.dirs = TRUE, 
                            full.names = TRUE)
  
  if(region_light == "ebs") {
    print("Combining EBS and NBS directories")
    region_dirs <- c(region_dirs, 
                     list.files(light_data_root, 
               pattern = "nbs", 
               recursive = TRUE, 
               include.dirs = TRUE, 
               full.names = TRUE))
  }
  
  vessel_dirs <- character()
  
  for(ii in 1:length(region_dirs)) {
    
    vessel_dirs <- c(vessel_dirs, list.files(region_dirs[ii],
                                             pattern = "v_",
                                             recursive = TRUE,
                                             include.dirs = TRUE,
                                             full.names = TRUE))
  }
  
  if(!is.null(omit_string)) {
    vessel_dirs <- vessel_dirs[!grepl(pattern = omit_string, 
                                      x = vessel_dirs)]
  }
  
  if(!dir.exists(here::here("imports"))) {
    dir.create(here::here("imports"))
  }
  
  print(paste0("Saving directory paths to ", here::here("imports", "directories.csv")))
  write.csv(data.frame(path = rev(vessel_dirs)), # Backwards so newest survey start processing first 
            file = here::here("imports", "directories.csv"), 
            row.names = FALSE)
}

#' Subsets light measurements from upcasts/downcasts
#'
#' Assigns light measurements to upcast or downcast based on downcast start time and downcast end time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @noRd

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


#' Find surface light measurements during casts
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Uses cast start and end times to find concurrent measurements obtained by the deck-mounted archival tag, including a buffer. The buffer is added around the start and end times, so a 30 second buffer = one minute.
#'
#' @param light.data Light measurements from the surface/deck archival tag.
#' @param cast.data Haul event times indicating the start and end times for net deployment and retrieval.
#' @param time.buffer Time buffer before and after the start of the cast.
#' @param agg.fun Function applied to calculate central tendency metric for light measurements sampled during the cast time window. Default trawllight::geometric.mean
#' @param ... Additional arguments passed to trawllight::convert_light()
#' @noRd

tlu_surface_light <- function(light.data, cast.data, time.buffer = 30, agg.fun = trawllight::geometric.mean, ...) {
  
  # Select and rename light and time columns
  if(ncol(light.data) >= 6) {
    light.data <- light.data[,5:6]
  } else if(ncol(light.data) == 4) {
    light.data <- light.data[,3:4]
    }else {
    warning(paste("surface_light: Deck light measurements not found for" , cast.data$vessel[1], "-", cast.data$cruise[1]))
    return(NULL)
  }
  
  colnames(light.data) <- c("surf_llight", "ctime")
  light.data$surf_trans_llight <- trawllight::convert_light(light.data$surf_llight, ...)
  light.data$vessel <- rep(cast.data$vessel[1], nrow(light.data))
  light.data$cruise <- rep(cast.data$cruise[1], nrow(light.data))
  
  for(i in 1:nrow(cast.data)) {
    # Assign upcast or downcast to tag time
    light.data$updown[light.data$ctime > (cast.data$downcast_start[i] - time.buffer) &
                        light.data$ctime < (cast.data$downcast_start[i] + time.buffer)] <- "downcast"
    light.data$updown[light.data$ctime > (cast.data$upcast_start[i] - time.buffer) &
                        light.data$ctime < (cast.data$upcast_end[i] + time.buffer)] <- "upcast"
    light.data$haul[light.data$ctime > (cast.data$downcast_start[i] - time.buffer) &
                      light.data$ctime < (cast.data$upcast_end[i] + time.buffer)] <- cast.data$haul[i]
  }
  
  # Remove measurements outside of time window
  light.data <- subset(light.data, !is.na(updown))
  
  llight <- aggregate(surf_trans_llight ~ haul + updown + vessel + cruise, data = light.data, FUN = agg.fun)
  
  ctime <- aggregate(ctime ~ haul + updown + vessel + cruise, data = light.data, FUN = mean)
  
  ctime$ctime <- lubridate::with_tz(ctime$ctime, "America/Anchorage")
  light.data <- dplyr::inner_join(llight, ctime, by = c("haul", "updown", "vessel", "cruise"))
  
  return(light.data)
}



#' Correct tag time in cases where offsets were incorrect
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Make adjustments to correct inconsistencies between tag time and survey time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @param survey RACE survey region as a character vector ("BS", "NBS", "AI", "GOA" or "SLOPE")
#' @keywords internal
#' @noRd

tlu_time_adjustments <- function(light.data, cast.data, survey) {
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  
  if(cast.data$vessel[1] == 134 & cast.data$cruise[1] == 200601 & region_light == "ebs") {
    print("Correcting 134-200601")
    light.data$ctime[lubridate::month(light.data$ctime) >=7 & lubridate::day(light.data$ctime) > 8] <- light.data$ctime[lubridate::month(light.data$ctime) >=7 & lubridate::day(light.data$ctime) > 8] - 3600*12 # Time off by 1 hour
  }
  
  return(light.data)
}

#' Combine haul data with data from upcasts, downcasts, and surface
#' 
#' Combine haul data with light data from upcasts, downcasts, and surface.
#' 
#' @param haul.dat Haul data as a data frame. Default NULL searches for haul_data file in the output directory.
#' @param survey RACE survey region as a character vector ("BS", "NBS", "AI", "GOA" or "SLOPE")
#' @param surface Surface light data as a data frame. Default NULL searches for surface data in the output directory.
#' @param downcasts Downcast data as a list that includes the light_data data frame. Default NULL searches for downcast data in the output directory.
#' @param upcasts Upcast data as a list that includes the light_Data data frame. Default NULL searaches for upcast data in the output directory.
#' @noRd

tlu_combine_casts <- function(haul.dat = NULL,
                              survey,
                              surface = NULL,
                              downcasts = NULL,
                              upcasts = NULL) {
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  
  if(is.null(haul.dat)) {
    haul.dat <- list.files(path = here::here("data"), 
               pattern = paste0(survey, "_haul_data"), 
               full.names = TRUE)
    print(paste0("tlu_combine_casts: Loading ", haul.dat))
    haul.dat <- readRDS(file = haul.dat )
  }
  if(is.null(surface)) {
    surface <- list.files(path = here::here("output"), 
                          pattern = paste0(region_light, "_surface"), 
                          full.names = TRUE)
    print(paste0("tlu_combine_casts: Loading ", surface))
    surface <- readRDS(file = surface)
  }
  
  if(is.null(downcasts)) {
    downcasts <- list.files(path = here::here("output"), 
               pattern = paste0(region_light, "_downcast"), 
               full.names = TRUE)
    print(paste0("tlu_combine_casts: Loading ", downcasts))
    downcasts <- readRDS(file = downcasts)
  }
  
  if(is.null(upcasts)) {
    upcasts <- list.files(path = here::here("output"), 
                         pattern = paste0(region_light, "_upcast"), 
                         full.names = TRUE)
    print(paste0("tlu_combine_casts: Loading ", upcasts))
    upcasts <- readRDS(file = upcasts)
  }

  down <- merge(downcasts$light_ratios,
                dplyr::select(
                  haul.dat,
                  vessel,
                  cruise,
                  haul,
                  start_latitude,
                  start_longitude
                )
  )
  up <- merge(upcasts$light_ratios,
              dplyr::select(haul.dat, vessel, cruise, haul, end_latitude, end_longitude)
  )
  names(up)[which(names(up) == "end_longitude")] <- "longitude"
  names(up)[which(names(up) == "end_latitude")] <- "latitude"
  
  names(down)[which(names(down) == "start_longitude")] <-
    "longitude"
  names(down)[which(names(down) == "start_latitude")] <-
    "latitude"
  
  casts <- rbind(down, up)
  casts <- subset(casts, !is.na(latitude) & !is.na(longitude))
  casts <-
    merge(
      casts,
      dplyr::select(
        haul.dat,
        vessel,
        cruise,
        haul,
        bottom_depth,
        stationid,
        haul_type
      )
    )
  casts <-
    dplyr::full_join(casts,
                     dplyr::select(surface, vessel, cruise, haul, updown, surf_trans_llight))
  casts <- subset(casts, !is.na(surface_time))
  
  return(casts)
}


#' Wrapper function to calculate residuals
#' 
#' Wrapper function over trawllight::tag_residuals_indirect() and trawllight_tag_residuals_direct(), which calculate residuals for the trawllight algorithm based on RACE's data structure.
#' 
#' @param input Input data, i.e., the data frame output by trawllight::tlu_combine_casts()
#' @param depth.bins Depth bins to use for calculating residuals. Default c(1,3,5) described by Rohan et al. (2020) and Rohan et al. (2021).
#' @param write.rds If TRUE, saves output to output/tag_residuals.rds
#' @noRd

tlu_calc_resids <- function(input, 
                            depth.bins = c(1, 3, 5),
                            write.rds = FALSE) {
  indirect_res <- trawllight::tag_residuals_indirect(
    x = input,
    formula = log10(trans_llight) ~ s(log10(PAR +
                                              0.001), bs = "cr"),
    depth.bins = depth.bins,
    utc.offset = -8,
    dark.adjust = 0.001,
    lat.col = "latitude",
    lon.col = "longitude",
    time.col = "surface_time",
    light.col = "trans_llight"
  )
  
  # Version of tag_residuals that uses linear regression
  direct_res <-
    trawllight::tag_residuals_direct(
      x = subset(indirect_res$resid_df, !is.na(surf_trans_llight)),
      formula = log10(trans_llight) ~ log10(surf_trans_llight) + interaction(vessel, cruise),
      water.col = "trans_llight",
      surface.col = "surf_trans_llight",
      depth.col = "cdepth",
      depth.bins = depth.bins
    )
  
  direct_res$resid_df <-
    dplyr::full_join(indirect_res$resid_df, direct_res$resid_df)
  
  # Tag residuals based on ratio between surface and water column
  direct_res$resid_df$surf_ratio <-
    log10(direct_res$resid_df$trans_llight / direct_res$resid_df$surf_trans_llight)
  
  #  saveRDS(
  #   direct_res$resid_df,
  #   file = here::here("output", "tag_residuals.rds")
  # )
  
  return(direct_res)
}

#' Applies a function to all casts in a data frame
#' 
#' For use processing AOPs from AFSC/RACE/GAP data structure.
#'
#' @param x Data frame
#' @param id.col Identification column
#' @param FUN Function to apply to all id.col variables
#' @param min.rows Minimum number of rows necessary for the function to run
#' @param silent Should id variables be printed to console.
#' @param ... Optional arguments passed to FUN.
#' @noRd

tlu_cast_wrapper <- function(x, id.col, FUN, min.rows = 4, silent = TRUE, ...) {
  
  # Create a unique ID column over which to loop
  x$id.col <- eval(parse(text = paste0("paste(", paste(paste0("x$", id.col), collapse = ","),")")))
  unique_ids <- unique(x$id.col)
  output.df <- NULL
  rowind <- 0
  
  for(i in 1:length(unique_ids)) {
    
    if(!silent) {
      print(unique_ids[i])
    }
    
    # Apply function to cast
    EEE <- subset(x, id.col == unique_ids[i])
    if(nrow(EEE) > min.rows) {
      EEE.out <- do.call(FUN, args = list(x = EEE, ...))
      
      if(is.null(output.df)) {
        output.df <- EEE.out
        rowind <- nrow(output.df)+1
        output.df[rowind:(nrow(output.df)+(nrow(x)*1.2)),] <- NA
      } else {
        output.df[rowind:(rowind + nrow(EEE.out)-1),] <- EEE.out[1:nrow(EEE.out),]
        rowind <- rowind + nrow(EEE.out)
      }
      
      if(i %% 1000 == 0) {
        print(i)
      }
    }
  }
  output.df <- output.df[,-which(names(output.df) == "id.col")]
  output.df <- output.df[-c(rowind:nrow(output.df)),]
  return(output.df)
}

#' Flag orientation
#' 
#' Assign orientation using trawllight::threshold_select()
#' 
#' @param resids Residuals data frame
#' @param direct.threshold Threshold for direct method residual error assignment as a numeric vector or ("auto") to use the method described in Rohan et al. (2020, 2021)
#' @param indirect.threshold Threshold for indirect method residual error assignment as a numeric vector or ("auto") to use the method described in Rohan et al. (2020, 2021)
#' @param direct.col Name of direct residual column. Default = "direct_residual"
#' @param indirect.col Name of indirect residual column. Default = "indirect_residual"
#' @noRd

tlu_flag_orientation <- function(resids,
                             direct.threshold  = "auto",
                             indirect.threshold = "auto",
                             direct.col = "direct_residual",
                             indirect.col = "indirect_residual") {
  dbins <- unique(resids$resid_df$cdepth)
  dt <- vector(length = length(dbins))
  it <- vector(length = length(dbins))
  if (direct.threshold == "auto") {
    for (i in 1:length(dbins)) {
      print("Calculating direct threshold (auto mode)")
      pdf(file = here::here("output", paste0("threshold_direct_", dbins[i], ".pdf")))
      sel <-
        resids$resid_df[which(resids$resid_df$cdepth == dbins[i]), which(names(resids$resid_df) == direct.col)]
      sel <- sel[!is.na(sel)]
      dir.threshold <-
        trawllight::threshold_select(vals = sel,
                                     method = "kernel",
                                     bandwidth.select = "ICV")
      print(dir.threshold$omit.rate)
      dt[i] <- dir.threshold$threshold
      print(dt[i])
      dev.off()
    }
  } else {
    dt <- direct.threshold
  }
  
  if (indirect.threshold == "auto") {
    for (i in 1:length(dbins)) {
      print("Calculating indirect threshold (auto mode)")
      pdf(file = here::here("output", paste0("threshold_indirect_", dbins[i], ".pdf")))
      sel <-
        resids$resid_df[which(resids$resid_df$cdepth == dbins[i]), which(names(resids$resid_df) == indirect.col)]
      sel <- sel[!is.na(sel)]
      ind.threshold <-
        trawllight::threshold_select(vals = sel,
                                     method = "kernel",
                                     bandwidth.select = "ICV")
      print(ind.threshold$omit.rate)
      it[i] <- ind.threshold$threshold
      print(it[i])
      dev.off()
    }
  } else {
    it <- direct.threshold
  }
  
  orient <-
    dplyr::full_join(
      trawllight::assign_orientation(
        resids$resid_df,
        formula = vessel + cruise + haul + updown ~ cdepth,
        resid.col = indirect.col,
        resid.cut = it,
        output.name = "indirect_orientation"
      ),
      trawllight::assign_orientation(
        resids$resid_df,
        formula = vessel + cruise + haul + updown ~ cdepth,
        resid.col = direct.col,
        resid.cut = dt,
        output.name = "direct_orientation"
      )
    )
  orient$orientation <- NA
  orient$orientation[orient$indirect_orientation == "Bad"] <- "Bad"
  orient$orientation[orient$indirect_orientation == "Good"] <-
    "Good"
  orient$orientation[orient$direct_orientation == "Good"] <- "Good"
  orient$orientation[orient$direct_orientation == "Bad"] <- "Bad"
  
  return(orient)
}

#' Combine orientation and quality for upcasts and downcasts
#' 
#' Combine orientation and quality for upcasts and downcasts
#' 
#' @param input Input data frame
#' @noRd

tlu_use_casts <- function(input) {
  qq <-
    unique(dplyr::select(input, vessel, cruise, haul, updown, quality, orientation))
  qq$final_qa <- 0
  qq$final_qa[qq$quality == 1 & qq$orientation == "Good"] <- 1
  qq <-
    reshape::cast(
      data = qq,
      vessel + cruise + haul ~ updown,
      value = "final_qa",
      fill = 0
    )
  return(qq)
}

#' QA/QC: Interpolated too bright
#' 
#' Omit casts where estimated near-surface light exceeds maximum observed
#' 
#' @param input Input data frame
#' @noRd

tlu_use_casts2 <- function(input) {
  input <-
    subset(input, trans_llight < max(input$trans_llight[input$cdepth == 1 &
                                                          !input$estimate_ref]))
  flag <- rep(F, nrow(input))
  flag[input$updown == "downcast" & input$downcast == 1] <- TRUE
  flag[input$updown == "upcast" & input$upcast == 1] <- TRUE
  return(input[which(flag),])
}

#' Merge attenuation profiles from casts that passed QA/QC
#' 
#' Merge attenuation profiles from upcast and downcasts with a data frame where data passed QA/QC
#' 
#' @param downcasts Data frame containing downcast attenuation profiles
#' @param upcasts Data frame containing upcast attenuation profiles
#' @param keep Data frame with QA/QC flags
#' @noRd

tlu_find_atten_profiles <- function(downcasts, upcasts, keep) {
  downcasts$updown <- "downcast"
  upcasts$updown <- "upcast"
  bbb <- rbind(downcasts, upcasts)
  bbb <- merge(bbb, keep)
  return(bbb)
}


#' Calculate relative light level and linear attenuation coefficient
#'
#' \code{light_proportion} calculates the proportion of light and the linear attenuation coefficient for each depth bin, relative to the highest observed light level for a cast.
#'
#' @param x A data frame which contains light and depth measurements.
#' @param light.col Name of the column containing light measurements.
#' @param depth.col Name of the column containing depth measurements.
#' @return Returns a data frame containing input values, the proportion of light relative to the highest observed light measurement as \code{light_ratio}, and the attenuation coefficient of diffuse downwelling irradiance relative to the highest observed light measurement as \code{k_linear}.
#' @keywords internal
#' @export

light_proportion <- function(x, light.col = "trans_llight", depth.col = "cdepth", ...) {
  
  names(x)[names(x) == light.col] <- "trans_llight"
  names(x)[names(x) == depth.col] <-  "cdepth"
  # Light ratio relative to shallowest bin.
  x$light_ratio <- x$trans_llight / max(x$trans_llight, na.rm = T)
  
  # Linear attenuation coefficient relative to shallowest bin.
  x$k_linear <- log(x$trans_llight/max(x$trans_llight)) / (min(x$cdepth) - x$cdepth)
  
  # Whole column attenuation coefficient.
  x$k_column <- rep((log(min(x$trans_llight)/max(x$trans_llight))) / (min(x$cdepth) - max(x$cdepth)), nrow(x))
  names(x)[names(x) == "trans_llight"] <- light.col
  names(x)[names(x) == "cdepth"] <-  depth.col
  return(x)
}

#' Write AOPs
#' 
#' Calculate cast-level variables (near-bottom optical depth, Z10, Z1) from filtered and interpolated data and write data output/final_stn_vars.rds; and write profile variables (Kd(Z)) that pass QA/QC checks to output/final_profile_vars.rds.
#' 
#'  @param survey Survey region "BS", "NBS", "GOA", "AI", or "SLOPE".
#'  @noRd 

tlu_calc_summary <- function(survey) {
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  
  print("tlu_calc_summary: Summarizing near bottom optical depth and writing to output/final_nbod.rds")
  aggregate(data = dplyr::select(
    readRDS(here::here("output", paste0("temp_", region_light, "_filtered_huds.rds"))),
    vessel,
    cruise,
    haul,
    updown,
    bottom_depth,
    cdepth,
    latitude,
    longitude,
    stationid),
    cdepth ~ bottom_depth + vessel + cruise + haul + updown + latitude + longitude + stationid,
    FUN = max) |>
    dplyr::inner_join(readRDS(here::here("output", paste0("temp_", region_light, "_filtered_huds.rds")))) |>
    dplyr::rename(near_bottom_optical_depth = optical_depth) |>
    dplyr::select(-downcast, -upcast, - k_linear, -k_column, -light_ratio) |>
    saveRDS(here::here("output", paste0("temp_", region_light, "_nbod.rds")))
  
  print("tlu_calc_summary: Finding attenuation profiles from casts that passed QA/QC and writing to output/final_profile_vars.rds")
  trawllight:::tlu_find_atten_profiles(
    downcasts = readRDS(list.files("output", pattern = paste0(region_light, "_downcast"), full.names = TRUE))$atten_values,
    upcasts =  readRDS(list.files("output", pattern = paste0(region_light, "_upcast"), full.names = TRUE))$atten_values,
    keep = unique(dplyr::select(
      readRDS(here::here("output", paste0("temp_", region_light, "_nbod.rds"))),
      vessel,
      cruise,
      haul,
      updown,
      latitude,
      longitude))) |>
    saveRDS(here::here("output", paste0(region_light, "_final_profile_vars.rds")))
  
  od_df <- readRDS(list.files(here::here("output"), pattern = paste0("temp_", region_light, "_filtered_huds"), full.names = TRUE))
  
  print("tlu_calc_summary: Estimating z[10]")  
  Z10_df <- trawllight:::tlu_cast_wrapper(
    x = od_df,
    id.col = c("vessel", "cruise", "haul", "updown"),
    FUN = trawllight::find_optical_depth,
    with.input = TRUE,
    return.col = "Z10",
    target.od = -1 * log(0.1)) |>
    dplyr::select(vessel, cruise, haul, updown, latitude, longitude, haul_type, bottom_depth, Z10) |>
    unique()
  
  print("tlu_calc_summary: Estimating z[1]")  
  Z1_df <- trawllight:::tlu_cast_wrapper(
    x = od_df,
    id.col = c("vessel", "cruise", "haul", "updown"),
    FUN = trawllight::find_optical_depth,
    with.input = TRUE,
    return.col = "Z1",
    target.od = -1 * log(0.01)) |>
    dplyr::select(vessel, cruise, haul, updown, latitude, longitude, haul_type, bottom_depth, Z1) |>
    unique()
  
  print("tlu_calc_summary: Combining temp_nbod, z1, and z10")
  dplyr::full_join(Z10_df, 
                   Z1_df, 
                   by = c("vessel", "cruise", "haul", "updown", "latitude", "longitude", "haul_type", "bottom_depth")) |>
    dplyr::full_join(
      readRDS(here::here("output", paste0("temp_", region_light, "_nbod.rds"))), 
      by = c("vessel", "cruise", "haul", "updown", "latitude", "longitude", "haul_type", "bottom_depth")) |>
    saveRDS(here::here("output", paste0(region_light, "_final_stn_vars.rds")))
  
}


#' Load and combine data frames from multiple .rds files
#' 
#' @param sel_dir Directory containing files to combine
#' @param pattern Character vector to search for in rds file names.
#' @param n_batch Number of files to combine at a time.
#' @noRd

.combine_rds_df <- function(sel_dir = here::here("output"), 
                            pattern,
                            n_batch = 5) {
  
  sel_files <- list.files(sel_dir, 
                          pattern = pattern,
                          full.names = TRUE)
  dat_out <- data.frame()
  int_index <- 1
  kk_iter <- length(sel_files)
  
  for(kk in 1:kk_iter) {
    dat_in <- try(readRDS(sel_files[kk]), silent = TRUE)
    
    if(class(dat_in) != "try-error") {
      dat_out <- dplyr::bind_rows(dat_out, dat_in)
    }
    
    if(kk == n_batch || kk == kk_iter) {
      assign(x = paste0(pattern, "_", int_index), value = dat_out)
      dat_out <- data.frame()
      int_index <- int_index + 1
    }
  }
  
  return(do.call(dplyr::bind_rows, 
                 mget(objects(pattern = pattern))))
  
}