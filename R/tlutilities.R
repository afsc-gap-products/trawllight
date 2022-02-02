#' Copy MK9 data and setup directory
#' 
#' Setup directory and retrieve data for processing MK9 data.
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param light_data_root Filepath to the data source
#' @export

tlu_setup_dir <- function(channel = NULL, survey, light_data_root = "G:/RACE_LIGHT/LightData/Data") {
  
  if(!(survey %in% c("BS", "NBS", "GOA", "AI"))) {
    stop(paste0("survey selection, ", survey, " invalid. Must be 'BS', 'NBS', 'GOA', or 'AI'."))
  }

  channel <- mk9process::get_connected(channel = channel)
  
  if(!dir.exists(here::here("data"))) {
    dir.create(here::here("data"))
  } else {
    warning(paste0("'data' directory already exists and was not overwritten."))
  }
    
  if(!dir.exists(here::here("output"))) {
    dir.create(here::here("output"))
  } else {
    warning(paste0("'output' directory already exists and was not overwritten."))
  }
  
  if(!dir.exists(here::here("imports"))) {
    dir.create(here::here("imports"))
  } else {
    warning(paste0("'imports' directory already exists and was not overwritten."))
  }
  
  trawllight::tlu_prep_haul_data(channel = channel,
                     survey = survey)
  
  trawllight::tlu_prep_dir_list(survey = survey,
                                light_data_root = light_data_root)
  
}

#' Get haul data and write to RDS
#' 
#' Run query to retrieve haul data from RACEBASE
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @export

tlu_prep_haul_data <- function(channel = NULL, 
                               survey) {
  
  survey <- toupper(survey)
  
  channel <- mk9process::get_connected(channel = channel)
  
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
#' @export

tlu_prep_dir_list <- function(survey, 
                              light_data_root = "G:/RACE_LIGHT/LightData/Data") {
  
  region_light <- c("ebs", "nbs", "goa", "ai")[match(survey, c("BS", "NBS", "GOA", "AI"))]
  
  region_dirs <- list.files(light_data_root, 
                            pattern = region_light, 
                            recursive = TRUE, 
                            include.dirs = TRUE, 
                            full.names = TRUE)
  
  vessel_dirs <- character()
  
  for(ii in 1:length(region_dirs)) {
    
    vessel_dirs <- c(vessel_dirs, list.files(region_dirs[ii],
                                             pattern = "v_",
                                             recursive = TRUE,
                                             include.dirs = TRUE,
                                             full.names = TRUE))
  }
  
  if(!dir.exists(here::here("imports"))) {
    dir.create(here::here("imports"))
  }
  
  print(paste0("Saving directory paths to ", here::here("imports", "directories.csv")))
  write.csv(data.frame(path = rev(vessel_dirs)), # Backwards so newest survey start processing first 
            file = here::here("imports", "directories.csv"), 
            row.names = FALSE)
}


#' Wrapper function to retrieve cast data and aggregate
#' 
#' Description goes here... writes binned cast data to an rds file after interpolating missing surface values.
#' 
#' @param directory_structure File path to csv file that lists directories to be processed.
#' @param survey RACE survey region as a character vector ("BS", "NBS", "AI", "GOA" or "SLOPE".
#' @param cast.dir Cast direction, either "upcast" or "downcast"
#' @param time.buffer Time buffer in seconds to add/subtract from upcast and downcast times. Default = 20
#' @param bin.size Passed to trawllight::filter_stepwise. Depth bin size for aggregating light measurements.  Default = 2
#' @param bin.gap Passed to trawllight::filter_stepwise. Maximum allowable gap in observations before a profile is considered to not meet continuity standards. Default = 6 (i.e., three bins if bin.size = 2)
#' @param agg.fun Function used to calculate summary statistic for a depth bin. Default = geometric.mean
#' @param silent Passed to... 
#' @param ... Optional arguments passed to filter_stepwise or calculate_attenuation.
#' @export

tlu_get_casts <- function(directory_structure,
                      survey,
                      cast.dir = "downcast",
                      time.buffer = 20,
                      bin.size = 2,
                      bin.gap = 6,
                      agg.fun = geometric.mean,
                      silent = TRUE,
                      ...) {
  
  cast.dir <- tolower(cast.dir)
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  out_path <- here::here("output", paste0(region_light, "_", cast.dir, ".rds"))
  
  
  # Batch load rds
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

  for(jj in 1:nrow(directory_structure)) {
    
    # Condition in case processing stops in the middle due to missing files in source path
    if(!file.exists(here::here("output", paste0("temp_resid_", jj, ".rds")))) {
      
      cast_dat <- trawllight::tlu_process_all(
        dir.structure = directory_structure[jj,],
        cast.dir = cast.dir,
        time.buffer = time.buffer,
        silent = silent,
        binsize = bin.size,
        bin.gap = bin.gap,
        agg.fun = agg.fun)
      
      saveRDS(cast_dat$light_ratios, file = here::here("output", paste0("temp_od_", jj, ".rds")))
      saveRDS(cast_dat$atten_values, file = here::here("output", paste0("temp_kd_", jj, ".rds")))
      saveRDS(cast_dat$loess_eval, file = here::here("output", paste0("temp_loess_", jj, ".rds")))
      saveRDS(cast_dat$resid_fit, file = here::here("output", paste0("temp_resid_", jj, ".rds")))
    }
    
  }
  
  print("tlu_get_casts: combining temp rds files from output")
  out_list <- list(light_ratios = .combine_rds_df(pattern = "temp_od", n_batch = 5),
                   atten_values = .combine_rds_df(pattern = "temp_kd", n_batch = 5),
                   loess_eval = .combine_rds_df(pattern = "temp_loess", n_batch = 5),
                   resid_fit = .combine_rds_df(pattern = "temp_resid", n_batch = 5))
  
  print(paste0("tlu_get_casts: writing output to ",  out_path))
  saveRDS(out_list, file = out_path)
  
  print("tlu_get_casts: removing temporary rds files from output")
  file.remove(list.files(here::here("output"), pattern = "temp", full.names = TRUE))
  
}


#' Wrapper function to retrieve cast data and aggregate
#' 
#' Description goes here... writes binned cast data to an rds file after interpolating missing surface values.
#' 
#' @param directory_structure File path to csv file that lists directories to be processed.
#' @param survey Survey as a character vector ("BS", "NBS", "GOA", "AI", or "SLOPE") 
#' @export

tlu_get_surface <- function(directory_structure, 
                            survey,
                            ...) {
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  out_path <- here::here("output", paste0(region_light, "_surface", ".rds"))
  
  surf <- trawllight::tlu_process_all_surface(dir.structure = directory_structure[,1], ...)

  saveRDS(surf, file = out_path)
  return(surf)
}

#' Wrapper function to run trawllight
#' 
#' Runs all of the trawllight functions in order to produces data products from upcast, downcast, and surface data.
#' 
#' @param rm.temp Should temporary files produced during data processing be removed?
#' @export

tlu_run_trawllight <- function(rm.temp = TRUE) {
  
  
  if(!file.exists(here::here("output", "temp_combined_huds.rds"))) {
    print("Combining haul, upcast, downcast, and surface (HUDS) data and writing to output/combined_huds.rds.")
    trawllight::tlu_combine_casts() |>
      saveRDS( here::here("output", "temp_combined_huds.rds"))
  }
  
  
  
  if(!file.exists(here::here("output", "temp_tag_residuals.rds"))) {
    print("Calculate residuals that are used for flagging archival tag orientation errors.")
    trawllight::tlu_calc_resids(input = 
                                  readRDS(here::here("output", "temp_combined_huds.rds"))) |>
      trawllight::tlu_flag_orientation(direct.threshold = "auto",
                                       indirect.threshold = "auto",
                                       direct.col = "surf_ratio",
                                       indirect.col = "indirect_residual") |>
      saveRDS(file = here::here("output", "temp_tag_residuals.rds"))
  }
  
  if(!file.exists(here::here("output", "temp_interp_huds.rds"))) {
    print("Estimating irradiance when missing/omitted from reference depth bin.")
    trawllight::tlu_cast_wrapper(
      x = readRDS(here::here("output", "temp_combined_huds.rds")),
      id.col = c(
        "vessel",
        "cruise",
        "haul",
        "updown",
        "quality",
        "k_column",
        "surf_trans_llight",
        "haul_type",
        "stationid",
        "bottom_depth",
        "orientation"
      ), silent = TRUE, FUN = "estimate_surface"
    ) |> 
      saveRDS(here::here("output", "temp_interp_huds.rds"))
  }
  
  if(!file.exists(here::here("output", "temp_filtered_huds.rds"))) {
    print("Combining residual and interpolated HUDS, seting QA/QC flags, and recalculating optical depth for interpolated data.")
    readRDS(here::here("output", "temp_tag_residuals.rds")) |>
      trawllight::tlu_use_casts() |> 
      merge(readRDS(here::here("output", "temp_interp_huds.rds"))) |>
      trawllight::tlu_use_casts2() |>
      cast_wrapper(id.col = c("vessel", "cruise", "haul", "updown"),
                   FUN = light_proportion) |>
      dplyr::mutate(optical_depth = log(1) - log(light_ratio)) |>
      saveRDS(here::here("output", "temp_filtered_huds.rds"))
  }
  
  print("Summarizing near bottom optical depth and writing to output/final_nbod.rds")
  aggregate(data = dplyr::select(
    readRDS(here::here("output", "temp_filtered_huds.rds")),
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
    dplyr::inner_join(readRDS(here::here("output", "temp_filtered_huds.rds"))) |>
    dplyr::rename(near_bottom_optical_depth = optical_depth) |>
    dplyr::select(-downcast, -upcast, - k_linear, -k_column, -light_ratio) |>
    saveRDS(here::here("output", "final_nbod.rds"))
  
  print("Finding attenuation profiles from casts that passed QA/QC and writing to output/final_atten.rds")
  tlu_find_atten_profiles(
    downcasts = readRDS(list.files("output", pattern = "downcast", full.names = TRUE))$atten_values,
    upcasts =  readRDS(list.files("output", pattern = "upcast", full.names = TRUE))$atten_values,
    keep = unique(dplyr::select(
      readRDS(here::here("output", "final_nbod.rds")),
      vessel,
      cruise,
      haul,
      updown,
      latitude,
      longitude))) |>
    saveRDS(here::here("output", "final_atten.rds"))
  
  if(rm.temp) {
    print("Removing temporary files.")
    file.remove(list.files(path = here::here("output"), pattern = "temp", full.names = TRUE))
  }
  
}

#' Correct tag time in cases where offsets were incorrect
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Make adjustments to correct inconsistencies between tag time and survey time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @export

tlu_time_adjustments <- function(light.data, cast.data) {
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

#' Combine haul data with data from upcasts, downcasts, and surface
#' 
#' Combine haul data with light data from upcasts, downcasts, and surface.
#' 
#' @param haul.dat Haul data as a data frame. Default NULL searches for haul_data file in the output directory.
#' @param surface Surface light data as a data frame. Default NULL searches for surface data in the output directory.
#' @param downcasts Downcast data as a list that includes the light_data data frame. Default NULL searches for downcast data in the output directory.
#' @param upcasts Upcast data as a list that includes the light_Data data frame. Default NULL searaches for upcast data in the output directory.
#' @export

tlu_combine_casts <- function(haul.dat = NULL,
                              surface = NULL,
                              downcasts = NULL,
                              upcasts = NULL) {
  
  if(is.null(haul.dat)) {
    haul.dat <- readRDS(file = list.files(path = here::here("data"), pattern = "haul_data", full.names = TRUE))
  }
  
  if(is.null(surface)) {
    surface <- readRDS(file = list.files(path = here::here("output"), pattern = "surface", full.names = TRUE))
  }
  
  if(is.null(downcasts)) {
    downcasts <- readRDS(file = list.files(path = here::here("output"), pattern = "downcast", full.names = TRUE))
  }
  
  if(is.null(upcasts)) {
    upcasts <- readRDS(file = list.files(path = here::here("output"), pattern = "upcast", full.names = TRUE))
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
#' @export

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
#' @export

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
#' @export

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
#' @export

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
#' @export

tlu_use_casts2 <- function(input) {
  input <-
    subset(input, trans_llight < max(input$trans_llight[input$cdepth == 1 &
                                                          !input$estimate_ref]))
  flag <- rep(F, nrow(input))
  flag[input$updown == "downcast" & input$downcast == 1] <- T
  flag[input$updown == "upcast" & input$upcast == 1] <- T
  return(input[which(flag),])
}

#' Merge attenuation profiles from casts that passed QA/QC
#' 
#' Merge attenuation profiles from upcast and downcasts with a data frame where data passed QA/QC
#' 
#' @param downcasts Data frame containing downcast attenuation profiles
#' @param upcasts Data frame containing upcast attenuation profiles
#' @param keep Data frame with QA/QC flags
#' @export

tlu_find_atten_profiles <- function(downcasts, upcasts, keep) {
  downcasts$updown <- "downcast"
  upcasts$updown <- "upcast"
  bbb <- rbind(downcasts, upcasts)
  bbb <- merge(bbb, keep)
  return(bbb)
}
