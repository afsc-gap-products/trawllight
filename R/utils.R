#' Applies a function to all casts in a data frame
#' 
#' For use processing AOPs from AFSC/RACE/GAP data structure.
#'
#' @param x Data frame
#' @param id.col Name of column containing unique identifier for variables (i.e. haul, vessel, cruise, up/downcast) apply function.
#' @param FUN Function to apply
#' @param min.rows Minimum number of rows necessary to apply functions.
#' @param silent Should messages be printed to console while the function is running?
#' @param ... Additional arguments passed to internal functions.
#' @export

cast_wrapper <- function(x, id.col, FUN, min.rows = 4, silent = TRUE, ...) {
  
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
        #output.df <- plyr::rbind.fill(output.df, EEE.out)
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

#' Open and bind rows of csv to data frame
#' 
#' For use with AFSC/RACE/GAP data structure.
#' 
#' @param directory Path to directory where csv files are stored.
#' @param string Unique pattern used in the csv files.
#' @export

csv_rbind <- function(directory, string) {
  file.list <- grep(pattern = string, x = dir(directory))
  
  if(substr(directory, nchar(directory), nchar(directory)) != "/") {
    if(substr(directory, nchar(directory), nchar(directory)) != "\\") {
      directory <- paste0(directory, "/")
    }
  }
  
  for(i in 1:length(file.list)) {
    if(i == 1) {
      out.df <- read.csv(file = paste0(directory, dir(directory)[file.list[i]]), stringsAsFactors = F)
      out.df$fname <- dir(directory)[file.list[i]]
    } else {
      out.comb <- read.csv(file = paste0(directory, dir(directory)[file.list[i]]), stringsAsFactors = F)
      out.comb$fname <- dir(directory)[file.list[i]]
      out.df <- rbind(out.df, out.comb)
    }
  }
  return(out.df)
}

#' Find and apply time offset to Mk9 light data
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure.
#'
#' @param light Data frame containing Mk9 data. Must include columns: ldate_time (POSIXct), ldepth (numeric)
#' @param mbt Data frame containing MBT data. Must include columns: date_time (POSIXct), depth (numeric)
#' @param try.offsets A vector of offsets to try. Default is seq(-8,8,0.5)
#' @param results.file Character vector specifying the filepath where information about the offset and correlation between Mk9 and MBT depths are stored.
#' @return The input light data frame with date_time adjusted according to the offset.
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

find_mk9_offset <- function(light, mbt, try.offsets = seq(-8,8,0.5), results.file = NULL) {
  
  if(!(("ldate_time" %in% names(light)) & ("ldepth" %in% names(light)))) {
    stop("find_mk9_offset: Columns named ldate_time and/or ldepth are missing from the light argument.")
  }
  
  if(!(("date_time" %in% names(mbt)) & ("depth" %in% names(mbt)))) {
    stop("find_mk9_offset: Columns named date_time and/or depth are missing from the mbt argument.")
  }
  # Initilize vector to store correlations from different offsets
  try.cor <- vector(length = length(try.offsets))
  
  # Loop through offsets
  for(i in 1:length(try.offsets)) {
    # Create offset to try
    light$ldate_time_offset <- light$ldate_time + try.offsets[i]*3600
    offset.df <- dplyr::inner_join(light, mbt, by = c("ldate_time_offset" = "date_time"))
    try.cor[i] <- cor(offset.df$ldepth, offset.df$depth)
  }
  
  print(paste0("Offset for Mk9 is " , try.offsets[which.max(try.cor)], " hrs, with correlation between Mk9 and MBT depth of ", try.cor[which.max(try.cor)], "."))
  
  # Remove try column
  light <- light[,-which(colnames(light) == "ldate_time_offset")]
  
  # Transform based on the best offset
  light$ldate_time <- light$ldate_time + try.offsets[which.max(try.cor)]*3600
  
  # Write offset and correlation to .txt file
  if(!is.null(results.file)) {
    fconn <- file(results.file)
    writeLines(c(as.character(Sys.Date()),
                 paste0("Offset: " , try.offsets[which.max(try.cor)], " hrs"),
                 paste0("Corr: ", try.cor[which.max(try.cor)])),fconn)
    close(fconn)
  }
  
  return(light)
  
}

#' Loop surface_light
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Iterates surface_light over directories containing deck**.csv files and CastTimes.csv files.
#'
#' @param dir.structure A character vector containing filepaths for all of the directories containing light measurements from surface/deck archival tags.
#' @param ... Additional arguments to be passed to surface_light or time_adjustments
#' @export

process_all_surface <- function(dir.structure, adjust.time = T, ...) {
  
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
      deck.data <- read.csv(file = deck.files[1], header = F)
      deck.data$ctime <- paste(deck.data[,1], deck.data[,2], sep = " ")
      
      # Import additional deck files if multiple exist in one directory
      if(length(deck.files) > 1) {
        for(b in 2:length(deck.files)) {
          deck.data <- rbind(deck.data, read.csv(file = deck.files[b], header = F))
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
        deck.data <- time_adjustments(light.data = deck.data,
                                      cast.data = cast.times)
      }
      if(nrow(deck.data) > 0) {
        
        # Find surface measurements
        surface_profiles <- surface_light(light.data = deck.data,
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

process_all <- function(dir.structure,
                        time.buffer = 20,
                        cast.dir = "Downcast",
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
        vert <- TLUtilities::vertical_profiles(light.data = corr_mk9hauls,
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
          
          if(cast.dir == "Upcast") {
            lr.out$surface_time <- casttimes$upcast_end[j]
          } else if(cast.dir == "Downcast") {
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

#' Raw light
#' 
#' For use processing AOPs from AFSC/RACE/GAP data structure. 
#' 
#' @param dir.structure Filepath to directory with corr_MK9hauls.csv and CastTimes.csv files
#' @param time.buffer Numeric vector 1L. Buffer around haul start time. 
#' @param silent Logical. Should messages be printed to the console during processing
#' @param coulumn Logical. Should light measurements from casts corr_MK9hauls.csv be retrieved?
#' @param adjust.time Should times be shifted to account for archival tags being set to the wrong time zones?
#' @param ... Additional arguments passed to internal functions.
#' @export

raw_light <- function(dir.structure,
                      time.buffer = 20,
                      silent = T,
                      column = T,
                      adjust.time = T,
                      ...) {
  
  # Step 1. Input directory and CastTimes and corr_MK9_hauls files
  ind <- 1
  out <- rep(-1, 1e7)
  
  for(i in 1:length(dir.structure)) {
    if(column) {
      if(file.exists(paste(dir.structure[i], "/corr_MK9hauls.csv", sep = "")) &
         file.exists(paste(dir.structure[i], "/CastTimes.csv", sep = ""))) {
        
        corr_mk9hauls <- read.csv(paste(dir.structure[i], "/corr_MK9hauls.csv", sep = ""))
        casttimes <- read.csv(paste(dir.structure[i], "/CastTimes.csv", sep = ""))
        
        corr_mk9hauls$ctime <- as.POSIXct(strptime(corr_mk9hauls$ctime,
                                                   format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
        casttimes$downcast_start <- as.POSIXct(strptime(casttimes$downcast_start,
                                                        format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time.buffer
        casttimes$downcast_end <- as.POSIXct(strptime(casttimes$downcast_end,
                                                      format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time.buffer
        casttimes$upcast_start <- as.POSIXct(strptime(casttimes$upcast_start,
                                                      format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time.buffer
        casttimes$upcast_end <- as.POSIXct(strptime(casttimes$upcast_end,
                                                    format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time.buffer
        
        for(j in 1:nrow(casttimes)) {
          if(!silent) {
            print(paste("Cruise: ", casttimes$cruise[j], ", Vessel: ", casttimes$vessel[j], ", Haul: ", casttimes$haul[j], sep = ""))
          }
          vert <- vertical_profiles(light.data = corr_mk9hauls,
                                    cast.data = subset(casttimes, haul == casttimes$haul[j]))
          
          
          if(nrow(vert) > 0) {
            out[ind:(ind+nrow(vert)-1)] <- vert$llight
            ind <- ind + nrow(vert)
          }
        }
      }
    } else {
      if(!file.exists(paste(dir.structure[i], "/CastTimes.csv", sep = ""))) {
        stop(paste("raw_light: CastTimes.csv not found in" , paste(dir.structure[i])))
      }
      
      # Import CastTImes
      print(paste("Processing", dir.structure[i]))
      casttimes <- read.csv(paste(dir.structure[i], "/CastTimes.csv", sep = ""))
      
      # Find names of deck files
      deck.files <- list.files(path = dir.structure[i], pattern = "^deck.*\\.csv", full.names = T)
      
      # Check for CastTimes
      if(length(deck.files) < 1) {
        warning(paste("raw_light: Deck light measurements not found in" , paste(dir.structure[i])))
      } else {
        
        #Import first deck file
        deck.data <- read.csv(file = deck.files[1], header = F)
        deck.data$ctime <- paste(deck.data[,1], deck.data[,2], sep = " ")
        
        # Import additional deck files if multiple exist in one directory
        if(length(deck.files) > 1) {
          for(b in 2:length(deck.files)) {
            deck.data <- rbind(deck.data, read.csv(file = deck.files[b], header = F))
          }
        }
        
        # Convert times into POSIXct
        deck.data$ctime <- as.POSIXct(strptime(deck.data$ctime, format = "%m/%d/%Y %H:%M:%S", tz = "America/Anchorage"))
        
        # Convert cast times to POSIXct format, add 30 second offset to each cast time to avoid truncating cast data
        casttimes$downcast_start <- as.POSIXct(strptime(casttimes$downcast_start,
                                                        format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
        casttimes$downcast_end <- as.POSIXct(strptime(casttimes$downcast_end,
                                                      format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
        casttimes$upcast_start <- as.POSIXct(strptime(casttimes$upcast_start,
                                                      format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
        casttimes$upcast_end <- as.POSIXct(strptime(casttimes$upcast_end,
                                                    format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
        
        if(adjust.time) {
          # Correct cases where there is a mismatch between survey time and tag time
          deck.data <- time_adjustments(light.data = deck.data,
                                        cast.data = casttimes)
        }
        if(nrow(deck.data) > 0) {
          
          
          if(ncol(deck.data) >= 6) {
            deck.data <- deck.data[,5:6]
            colnames(deck.data) <- c("surf_llight", "ctime")
            deck.data$vessel <- rep(casttimes$vessel[1], nrow(deck.data))
            deck.data$cruise <- rep(casttimes$cruise[1], nrow(deck.data))
            
            for(j in 1:nrow(casttimes)) {
              # Assign upcast or downcast to tag time
              deck.data$updown[deck.data$ctime > (casttimes$downcast_start[j] - time.buffer) &
                                 deck.data$ctime < (casttimes$downcast_start[j] + time.buffer)] <- "Downcast"
              deck.data$updown[deck.data$ctime > (casttimes$upcast_start[j] - time.buffer) &
                                 deck.data$ctime < (casttimes$upcast_end[j] + time.buffer)] <- "Upcast"
              deck.data$haul[deck.data$ctime > (casttimes$downcast_start[j] - time.buffer) &
                               deck.data$ctime < (casttimes$upcast_end[j] + time.buffer)] <- casttimes$haul[j]
            }
            
            # Remove measurements outside of time window
            deck.data <- subset(deck.data, !is.na(updown))
            out[ind:(ind+nrow(deck.data)-1)] <- deck.data$surf_llight
            ind <- ind + nrow(deck.data)
          }
        }
        
      }
      
    }
  }
  return(out[out > 0])
}

#' Find surface light measurements during casts
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Uses cast start and end times to find concurrent measurements obtained by the deck-mounted archival tag, including a buffer. The buffer is added around the start and end times, so a 30 second buffer = one minute.
#'
#' @param light.data Light measurements from the surface/deck archival tag.
#' @param cast.data Haul event times indicating the start and end times for net deployment and retrival.
#' @param time.buffer Time buffer before and after the start of the cast.
#' @param agg.fun Function applied to calculate central tendency metric for light measurements sampled during the cast time window. Default trawllight::geometric.mean
#' @param ... Additional arguments passed to internal functions.
#' @export


surface_light <- function(light.data, cast.data, time.buffer = 30, agg.fun = trawllight::geometric.mean, ...) {
  
  # Select and rename light and time columns
  if(ncol(light.data) >= 6) {
    light.data <- light.data[,5:6]
    colnames(light.data) <- c("surf_llight", "ctime")
    light.data$surf_trans_llight <- convert_light(light.data$surf_llight, ...)
    light.data$vessel <- rep(cast.data$vessel[1], nrow(light.data))
    light.data$cruise <- rep(cast.data$cruise[1], nrow(light.data))
    
    for(i in 1:nrow(cast.data)) {
      # Assign upcast or downcast to tag time
      light.data$updown[light.data$ctime > (cast.data$downcast_start[i] - time.buffer) &
                          light.data$ctime < (cast.data$downcast_start[i] + time.buffer)] <- "Downcast"
      light.data$updown[light.data$ctime > (cast.data$upcast_start[i] - time.buffer) &
                          light.data$ctime < (cast.data$upcast_end[i] + time.buffer)] <- "Upcast"
      light.data$haul[light.data$ctime > (cast.data$downcast_start[i] - time.buffer) &
                        light.data$ctime < (cast.data$upcast_end[i] + time.buffer)] <- cast.data$haul[i]
    }
    
    # Remove measurements outside of time window
    light.data <- subset(light.data, !is.na(updown))
    llight <- aggregate(surf_trans_llight ~ haul + updown + vessel + cruise, data = light.data, FUN = agg.fun)
    ctime <- aggregate(ctime ~ haul + updown + vessel + cruise, data = light.data, FUN = mean)
    ctime$ctime <- lubridate::with_tz(ctime$ctime, "America/Anchorage")
    light.data <- merge(llight, ctime)
    
    return(light.data)
  } else {
    warning(paste("surface_light: Deck light measurements not found for" , cast.data$vessel[1], "-", cast.data$cruise[1]))
    return(NULL)
  }
}

#' Correct tag time in cases where offsets were incorrect
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Make adjustments to correct inconsistencies between tag time and survey time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @export

time_adjustments <- function(light.data, cast.data) {
  # Add vessel/cruise combination corrections for processing.
  # Offsets for tags set to the wrong time zone
  if(cast.data$cruise[1] == 201601) {
    print("Correcting 2016")
    light.data$ctime <- light.data$ctime + 3600 # Time off by 1 hour
  }
  
  if(cast.data$vessel[1] == 94 & cast.data$cruise[1] == 201501) {
    print("Correcting 94-201501")
    light.data$ctime <- light.data$ctime - 3600
  }
  
  if(cast.data$vessel[1] == 94 & cast.data$cruise[1] == 201401) {
    print("Correcting 94-201401")
    light.data$ctime <- light.data$ctime - 3600
  }
  
  if(cast.data$vessel[1] == 162 & cast.data$cruise[1] == 201101) {
    print("Correcting 162-201101")
    light.data$ctime <- light.data$ctime - (3600*8)
  }
  
  if(cast.data$vessel[1] == 134 & cast.data$cruise[1] == 200601) {
    print("Correcting 134-200601")
    light.data$ctime[lubridate::month(light.data$ctime) >=7 & lubridate::day(light.data$ctime) > 8] <- light.data$ctime[lubridate::month(light.data$ctime) >=7 & lubridate::day(light.data$ctime) > 8] - 3600*12 # Time off by 1 hour
  }
  
  return(light.data)
}

#' Subsets light measurements from upcasts/downcasts
#'
#' Assigns light measurements to upcast or downcast based on downcast start time and downcast end time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @export

vertical_profiles <- function(light.data, cast.data) {
  
  # Remove surface data
  light.data <- subset(light.data, cdepth >= 0)
  
  # Create empty columns for cast direction (updown) and haul number
  light.data$updown <- rep(NA, nrow(light.data))
  light.data$haul <- rep(NA, nrow(light.data))
  
  # Assign upcast or downcast to times
  light.data$updown[light.data$ctime > cast.data$downcast_start & light.data$ctime < cast.data$downcast_end[1]] <- "Downcast"
  light.data$updown[light.data$ctime > cast.data$upcast_start & light.data$ctime < cast.data$upcast_end[1]] <- "Upcast"
  light.data$haul[light.data$ctime > cast.data$downcast_start & light.data$ctime < cast.data$upcast_end[1]] <- cast.data$haul[1]
  
  # Remove on-bottom and errant sampling not from casts
  light.data <- subset(light.data, is.na(updown) == F)
  return(light.data)
}

#' Convert radiometric energy to photon flux density
#' 
#' Convert radiometric energy in watts to micromoles of photons per meter squared per second.
#' 
#' @param x Numeric vector. Energy for a wavelength in units of watts.
#' @param wavelength Numeric vector.Wavelength in nanometers.
#' @return Numeric vector of energy in photon flux density.
#' @export

energy_to_quanta <- function(x, wavelength) {
  return(x/(1e-6*6.02214076*10^23*(3e8/(wavelength*10e-9))*6.63e-34))
  
  
}

#' Photon flux density to radiometric energy 
#' 
#' Convert quantum units of micromoles of photons per meter squared per second to radiometric energy (watts per meter squared per second)
#' 
#' Convert energy to quanta based on wavelength.
#' @param x Numeric vector. Energy for a wavelength in units of watts.
#' @param wavelength Numeric vector.Wavelength in nanometers.
#' @return Numeric vector of energy in radiometric energy in watts.
#' @export

quanta_to_energy <- function(x, wavelength) {
  return(x*1e-6*6.02214076*10^23*(3e8/(wavelength*10e-9))*6.63e-34)
}