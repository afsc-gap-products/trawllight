#' Find and apply time offset to deck tag data
#' 
#' Calculates correlations between predicted and observed irradiance to calculate an archival tag offset for cases where tag time doesn't match survey time. Iterates through surface irradiance predictions from fishmethods::astrocalc4r based on location data and offset time.
#'
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @param results_file Character vector indicating a filepath where offset results should be stored.
#' @return Writes time-shifted surface values to corr_deck.csv
#' @export

mk9_find_surface_offset <- function(channel = NULL, survey, vessel, cruise, try_offsets = seq(-8,8,1), results_file = NULL) {
  
  channel <- trawllight:::get_connected(channel)
  
  light_dir <- here::here("data", "mk9", survey, cruise, vessel)
  region_light <- c("ebs", "nbs", "goa", "ai")[match(survey, c("BS", "NBS", "GOA", "AI"))]
  
  offset_cor <- numeric(length = length(try_offsets))
  
  if(survey == "NBS") {
    survey <- "BS"
  }
  
  haul_dat <- RODBC::sqlQuery(channel = channel,
                              query = paste0("select vessel, cruise, haul, stationid, start_latitude, start_longitude, end_latitude, end_longitude, region from racebase.haul where cruise > 200400 and region = '", survey, "'")) |>
    dplyr::mutate(LATITUDE = (START_LATITUDE + END_LATITUDE)/2,
                  LONGITUDE = (START_LONGITUDE + END_LONGITUDE)/2) |>
    dplyr::select(VESSEL, CRUISE, HAUL, STATIONID, LATITUDE, LONGITUDE)
  
  names(haul_dat) <- tolower(names(haul_dat))
  
  casttimes <- read.csv(file = paste0(light_dir, "/CastTimes.csv"))
  
  casttimes <- dplyr::inner_join(casttimes, haul_dat, by = c("vessel", "cruise", "haul"))
  
  casttimes$downcast_start <- as.POSIXct(strptime(casttimes$downcast_start,
                                                  format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
  casttimes$upcast_end <- as.POSIXct(strptime(casttimes$upcast_end,
                                              format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
  
  # Find names of deck files
  deck_files <- list.files(path = light_dir, pattern = "^deck.*\\.csv", full.names = T)
  
  # Check for CastTimes
  if(length(deck_files) < 1) {
    warning(paste("mk9_find_surface_offset: Deck light measurements not found in" , light_dir))
  } else {
    
    #Import first deck file
    deck_dat <- read.csv(file = deck_files[1], header = F, skip = 3)
    deck_dat$ctime <- paste(deck_dat[,1], deck_dat[,2], sep = " ")
    
    # Import additional deck files if multiple exist in one directory
    if(length(deck_files) > 1) {
      for(b in 2:length(deck_files)) {
        add_deck <- read.csv(file = deck_files[b], header = F, skip = 3)
        add_deck$ctime <- paste(add_deck[,1], add_deck[,2], sep = " ")
        deck_dat <- rbind(deck_dat, add_deck)
      }
    }
    
    if(ncol(deck_dat) == 6) {
      deck_dat <- deck_dat[,3:6]
    } else if(ncol(deck_dat) == 7){
      deck_dat <- deck_dat[,c(3:5, 7)] 
    } else {
      stop(paste0("mk9_find_surface_offset: Unexpected number of columns (", ncol(deck_dat), "). in deck data from cruise ", casttimes$cruise[1], " vessel ", casttimes$vessel, ". Expected: 6. May need to add processing option for the deck file configuration."))
    }
    
    names(deck_dat) <- c("ldepth", "ltemp", "llight", "ctime")
    
    # Convert times into POSIXct
    deck_dat$ctime <- as.POSIXct(strptime(deck_dat$ctime, format = "%m/%d/%Y %H:%M:%S", tz = "America/Anchorage"))
    
    if(casttimes$vessel[1] == 134 & casttimes$cruise[1] == 200601) {
      print("mk9_find_surface_offset: Additional correction for vessel 134, cruise 200601; -12 hours for data after July 8, 2006")
      deck_dat$ctime[lubridate::month(deck_dat$ctime) >=7 & lubridate::day(deck_dat$ctime) > 8] <- 
        deck_dat$ctime[lubridate::month(deck_dat$ctime) >=7 & lubridate::day(deck_dat$ctime) > 8] - 3600*12
    }
    
    sample_surf_light <- data.frame()
    
    for(jj in 1:nrow(casttimes)) {
      
      sub_light <- dplyr::filter(deck_dat, ctime >= casttimes$downcast_start[jj] & ctime <= casttimes$upcast_end[jj])
      
      if(nrow(sub_light) > 0) {
        if(nrow(sub_light) >= 30) {
          sub_light <- sub_light[sample(1:nrow(sub_light), size = 30, replace = FALSE), ]
        }
        
        sub_light$latitude <- casttimes$latitude[jj]
        sub_light$longitude <- casttimes$longitude[jj]
        sub_light$vessel <- casttimes$vessel[jj]
        sub_light$cruise <- casttimes$cruise[jj]
        sub_light$haul <- casttimes$haul[jj]
        
        sample_surf_light <- dplyr::bind_rows(sample_surf_light, sub_light)
        
      }
      
    }
    
    if(nrow(sample_surf_light) > 0) {
      for(kk in 1:length(try_offsets)) {
        
        sample_surf_light$ctime_adj_utc <- lubridate::with_tz(sample_surf_light$ctime + try_offsets[kk] * 3600, 
                                                              tzone = "UTC")
        
        sample_surf_light$ac_PAR <- fishmethods::astrocalc4r(day = lubridate::day(sample_surf_light$ctime_adj_utc),
                                                             month = lubridate::month(sample_surf_light$ctime_adj_utc),
                                                             year = lubridate::year(sample_surf_light$ctime_adj_utc),
                                                             hour = lubridate::hour(sample_surf_light$ctime_adj_utc) + 
                                                               lubridate::minute(sample_surf_light$ctime_adj_utc)/60,
                                                             timezone = rep(0, nrow(sample_surf_light)),
                                                             lat = sample_surf_light$latitude,
                                                             lon = sample_surf_light$longitude, 
                                                             withinput = FALSE)$PAR
        
        offset_cor[kk] <- cor(sample_surf_light$ac_PAR, 
                              trawllight::convert_light(sample_surf_light$llight, 
                                                        convert.method = "wc"))
        
      }
      
      best_offset <- try_offsets[which.max(offset_cor)]
      best_cor <- offset_cor[which.max(offset_cor)]
      
      pdf(file = paste0(light_dir, "/plot_surface_offset.pdf"))
      plot(try_offsets, offset_cor, xlab = "Offset (hrs)", 
           ylab = "Correlation", 
           ylim = c(0,1),
           main = paste0("Cruise: ", casttimes$cruise[1], ", Vessel: ", casttimes$vessel[1]))
      abline(v = best_offset, lty = 2)
      points(best_offset, best_cor, col = "red", pch = 16)
      dev.off()
      
      # Write offset and correlation to .txt file
      if(is.null(results_file)) {
        results_file <- paste0(light_dir, "/log_surface_offset.txt")
      }
      fconn <- file(results_file)
      writeLines(c(best_offset,
                   as.character(Sys.Date()),
                   paste0("Offset: " , best_offset, " hrs"),
                   paste0("Corr: ", best_cor )),fconn)
      close(fconn)
      
      print(paste0("mk9_find_surface_offset: Applying ", best_offset, " hr offset to deck data from cruise ", casttimes$cruise[1], ", vessel ", casttimes$vessel[1], ". Correlation: ", round(best_cor, 3)))
      deck_dat$ctime <- deck_dat$ctime + 3600 * best_offset
      
      # Plot best offset
      
      sample_surf_light$best_ctime_adj_utc <- lubridate::with_tz(sample_surf_light$ctime + 3600 * best_offset, 
                                                                 tzone = "UTC")
      
      sample_surf_light$llight_converted <- trawllight::convert_light(sample_surf_light$llight, 
                                                                      convert.method = "wc")
      
      sample_surf_light$best_ac_PAR <- fishmethods::astrocalc4r(day = lubridate::day(sample_surf_light$best_ctime_adj_utc),
                                                                month = lubridate::month(sample_surf_light$best_ctime_adj_utc),
                                                                year = lubridate::year(sample_surf_light$best_ctime_adj_utc),
                                                                hour = lubridate::hour(sample_surf_light$best_ctime_adj_utc) + 
                                                                  lubridate::minute(sample_surf_light$best_ctime_adj_utc)/60,
                                                                timezone = rep(0, nrow(sample_surf_light)),
                                                                lat = sample_surf_light$latitude,
                                                                lon = sample_surf_light$longitude, 
                                                                withinput = FALSE)$PAR
      
      best_offset_lm <- lm(log10(llight_converted+1e-5) ~ log10(best_ac_PAR+1e-5), data = sample_surf_light)
      
      sample_surf_light$llight_resid <- residuals(best_offset_lm)
      
      pdf(paste0(light_dir, "/plot_surface_residuals_by_day.pdf"))
      print(ggplot() +
              geom_boxplot(data = sample_surf_light,
                         mapping = aes(x = lubridate::day(ctime + 3600 * best_offset),
                                       y = llight_resid,
                                       group = lubridate::day(ctime + 3600 * best_offset))) +
              geom_smooth(data = sample_surf_light,
                          mapping = aes(x = lubridate::day(ctime + 3600 * best_offset),
                                        y = llight_resid)) +
              geom_hline(yintercept = 0, linetype = 2) +
              facet_wrap(~lubridate::month(best_ctime_adj_utc, 
                                           label = TRUE), scales = "free_x", nrow = 3) +
              scale_x_continuous(name = "Day") +
              scale_y_continuous(name = "Residual"))
      dev.off()
      
      print(paste0("mk9_find_surface_offset: Writing time-shifted surface light to ", light_dir, "/corr_deck.csv"))
      write.csv(x = deck_dat, paste0(light_dir, "/corr_deck.csv"), row.names = FALSE)
    } else {
      stop("mk9_find_surface_offset: No surface data found during cast times. Check that deck csv files contain data for the chosen vessel/cruise.")
    }
    
  }
  
}