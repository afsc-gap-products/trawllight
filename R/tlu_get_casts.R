#' Wrapper function to retrieve cast data and aggregate
#' 
#' Description goes here... writes binned cast data to an rds file after interpolating missing surface values.
#' 
#' @param directory_structure File path to csv file that lists directories to be processed.
#' @param survey RACE survey region as a character vector ("BS", "NBS", "AI", "GOA" or "SLOPE")
#' @param cast.dir Cast direction, either "upcast" or "downcast"
#' @param time.buffer Time buffer in seconds to add/subtract from upcast and downcast times. Default = 20
#' @param bin.size Passed to trawllight::filter_stepwise. Depth bin size for aggregating light measurements.  Default = 2
#' @param bin.gap Passed to trawllight::filter_stepwise. Maximum allowable gap in observations before a profile is considered to not meet continuity standards. Default = 6 (i.e., three bins if bin.size = 2)
#' @param agg.fun Function used to calculate summary statistic for a depth bin. Default = geometric.mean
#' @param silent Passed to... 
#' @param ... Optional arguments passed to filter_stepwise or calculate_attenuation.
#' @export

tlu_get_casts <- function(directory_structure = NULL,
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
  
  if(is.null(directory_structure)) {
    directory_structure <- read.csv(file = here::here("imports", "directories.csv"))
  }
  
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
      
      cast_dat <- trawllight:::tlu_process_all(
        dir.path = directory_structure[jj,],
        cast.dir = cast.dir,
        time.buffer = time.buffer,
        silent = silent,
        binsize = bin.size,
        bin.gap = bin.gap,
        agg.fun = agg.fun)
      
      
      if(!is.null(cast_dat)) {
        saveRDS(cast_dat$light_ratios, file = here::here("output", paste0("temp_od_", jj, ".rds")))
        saveRDS(cast_dat$atten_values, file = here::here("output", paste0("temp_kd_", jj, ".rds")))
        saveRDS(cast_dat$loess_eval, file = here::here("output", paste0("temp_loess_", jj, ".rds")))
        saveRDS(cast_dat$resid_fit, file = here::here("output", paste0("temp_resid_", jj, ".rds")))
      }
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