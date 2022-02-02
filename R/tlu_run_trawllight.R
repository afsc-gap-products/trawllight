#' Wrapper function to run trawllight
#' 
#' Runs all of the trawllight functions in order to produces data products from upcast, downcast, and surface data.
#' 
#' @param rm.temp Should temporary files produced during data processing be removed?
#' @export

tlu_run_trawllight <- function(rm.temp = TRUE) {
  
  
  if(!file.exists(here::here("output", "temp_combined_huds.rds"))) {
    print("Combining haul, upcast, downcast, and surface (HUDS) data and writing to output/combined_huds.rds.")
    trawllight:::tlu_combine_casts() |>
      saveRDS( here::here("output", "temp_combined_huds.rds"))
  }
  
  
  
  if(!file.exists(here::here("output", "temp_tag_residuals.rds"))) {
    print("Calculate residuals that are used for flagging archival tag orientation errors.")
    trawllight:::tlu_calc_resids(input = 
                                   readRDS(here::here("output", "temp_combined_huds.rds"))) |>
      trawllight:::tlu_flag_orientation(direct.threshold = "auto",
                                        indirect.threshold = "auto",
                                        direct.col = "surf_ratio",
                                        indirect.col = "indirect_residual") |>
      saveRDS(file = here::here("output", "temp_tag_residuals.rds"))
  }
  
  if(!file.exists(here::here("output", "temp_interp_huds.rds"))) {
    print("Estimating irradiance when missing/omitted from reference depth bin.")
    trawllight:::tlu_cast_wrapper(
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
      trawllight:::tlu_use_casts() |> 
      merge(readRDS(here::here("output", "temp_interp_huds.rds"))) |>
      trawllight:::tlu_use_casts2() |>
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
  trawllight:::tlu_find_atten_profiles(
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