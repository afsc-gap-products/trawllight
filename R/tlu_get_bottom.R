#' Retrieve bottom light data
#' 
#' Get bottom light data
#' 
#' @param directory_structure File path to csv file that lists directories to be processed.
#' @param survey RACE survey region as a character vector ("BS", "NBS", "AI", "GOA" or "SLOPE")
#' @param time_buffer Time buffer in seconds to add/subtract from upcast and downcast times. Default = 20
#' @param agg_fun Function applied to calculate central tendency metric for light measurements sampled during the cast time window. Default = mean
#' @export

tlu_get_bottom <- function(directory_structure, agg_fun = mean, time_buffer = 30, survey, ...) {
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  
  out_path <- here::here("output", paste0(region_light, "_bottom_light.rds"))
  out_path_mean <- here::here("output", paste0(region_light, "_mean_bottom.rds"))
  
  for(ii in 1:nrow(directory_structure)) {
    
    # Read-in directory and CastTimes and corr_MK9_hauls files, select near-bottom light data
    
    if(!file.exists(here::here("output", paste0("temp_bottom_light_", ii, ".rds")))) {  
    
    light_path = paste0(directory_structure[ii, ], "/corr_MK9hauls.csv")
    casttimes_path = paste0(directory_structure[ii, ], "/CastTimes.csv")
    
    # Import light files from directory
    print(paste0("tlu_bottom_light: Loading ", light_path))
    light_data <- read.csv(light_path)
    
    if("lcond" %in% names(light_data)) {
      if(class(light_data$lcond) != "integer") {
        light_data$lcond <- as.integer(light_data$lcond)
      }
    }
    
    print(paste0("tlu_bottom_light: Loading ", casttimes_path))
    cast_data <- read.csv(casttimes_path)
    
    # Convert cast times to POSIXct format, add time buffer.
    cast_data$downcast_end <- as.POSIXct(strptime(cast_data$downcast_end, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) + time_buffer
    cast_data$upcast_start <- as.POSIXct(strptime(cast_data$upcast_start, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage")) - time_buffer 
    
    # Convert light data to POSIXct format
    light_data$ctime <- as.POSIXct(strptime(light_data$ctime, format = "%Y-%m-%d %H:%M:%S", tz = "America/Anchorage"))
    
    # Create empty columns for cast direction and haul number
    light_data$bottom <- rep(NA, nrow(light_data))
    light_data$haul <- rep(NA, nrow(light_data))

    
    for(jj in 1:nrow(cast_data)) {
      # Assign upcast or downcast to times
      light_data$bottom[light_data$ctime > cast_data$downcast_end[jj] & light_data$ctime < cast_data$upcast_start[jj]] <- "bottom"
      light_data$haul[light_data$ctime > cast_data$downcast_end[jj] & light_data$ctime < cast_data$upcast_start[jj]] <- cast_data$haul[jj]
      
    }
    
    # Remove profile data
    bottom_light_df <- subset(light_data, is.na(bottom) == F)
    
    # Calculate geometric mean near-bottom light
    bottom_light_df$trans_llight <- trawllight::convert_light(bottom_light_df$llight, ...)
  
    if(!is.null(bottom_light_df)) {
      saveRDS(bottom_light_df, file = here::here("output", paste0("temp_bottom_light_", ii, ".rds")))
    }
    
    }
    
  }
  
  output_bottom_light <- trawllight:::.combine_rds_df(pattern = "temp_bottom_light_", n_batch = 5)
  
  loc_dat <- readRDS(file = here::here("output", paste0(region_light, "_tag_residuals.rds"))) |>
    dplyr::select(vessel, cruise, haul, latitude, longitude, bottom_depth, stationid) |>
    unique() |>
    dplyr::group_by(vessel, cruise, haul, bottom_depth, stationid) |>
    dplyr::summarise(latitude = mean(latitude),
                     longitude = mean(longitude))
  
  print("tlu_get_bottom: combining bottom light with orientation flags.")
  output_bottom_light <- readRDS(file = here::here("output", paste0(region_light, "_tag_residuals.rds"))) |>
    dplyr::select(vessel, cruise, haul, orientation) |>
    unique() |>
    dplyr::mutate(orientation = as.numeric(as.character(factor(orientation, levels = c("Bad", "Good"), labels = c(0,1))))) |>
    dplyr::group_by(vessel, cruise, haul) |>
    dplyr::summarise(orientation = mean(orientation), .groups = "keep") |>
    dplyr::ungroup() |>
    dplyr::right_join(output_bottom_light |> 
                        dplyr::select(-ltemp, -ldepth, -lcond),
                      by = c("vessel", "cruise", "haul")) |>
    dplyr::inner_join(loc_dat, 
                      by = c("vessel", "cruise", "haul"))
  
  print(paste0("tlu_get_bottom: writing output to ",  out_path))
  saveRDS(output_bottom_light, file = out_path)
  
  print(paste0("tlu_get_bottom: writing mean bottom light to ", out_path_mean))
  output_bottom_light |>
    dplyr::group_by(vessel, cruise, haul) |>
    dplyr::summarise(tag_bottom_light = agg_fun(llight),
                     bottom_light = agg_fun(trans_llight),
                     cdepth = mean(cdepth),
                     orientation = mean(orientation)) |>
    dplyr::ungroup() |>
    dplyr::inner_join(loc_dat, by = c("vessel", "cruise", "haul")) |>
    saveRDS(file = out_path_mean)
    
    print("tlu_get_bottom: removing temporary rds files from output")
    file.remove(list.files(here::here("output"), pattern = "temp_bottom_light", full.names = TRUE))
    
}