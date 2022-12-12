#' Wrapper function to run trawllight
#' 
#' Runs all of the trawllight functions in order to produces data products from upcast, downcast, and surface data.
#' 
#' @param rm.temp Should temporary files produced during data processing be removed?
#' @param survey Survey region "BS", "NBS", "GOA", "AI", or "SLOPE".
#' @export

tlu_run_trawllight <- function(rm.temp = TRUE, survey) {
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  
    print(paste0("tlu_run_trawllight: Combining haul, upcast, downcast, and surface (HUDS) data. Writing temp file to output/", 
                 region_light, 
                 "_temp_combined_huds.rds."))
    trawllight:::tlu_combine_casts(survey = survey) |>
      saveRDS(here::here("output", paste0("temp_", region_light, "_combined_huds.rds")))
  
    print("tlu_run_trawllight: Calculating residuals for flagging archival tag orientation errors.")
    trawllight:::tlu_calc_resids(input = 
                                   readRDS(here::here("output", paste0("temp_", region_light, "_combined_huds.rds")))) |>
      trawllight:::tlu_flag_orientation(direct.threshold = "auto",
                                        indirect.threshold = "auto",
                                        direct.col = "surf_ratio",
                                        indirect.col = "indirect_residual") |>
      saveRDS(file = here::here("output", paste0(region_light, "_tag_residuals.rds")))
  
    print("tlu_run_trawllight: Estimating reference irradiance if missing/omitted.")
    trawllight:::tlu_cast_wrapper(
      x = readRDS(here::here("output", paste0("temp_", region_light, "_combined_huds.rds"))),
      id.col = c(
        "vessel",
        "cruise",
        "haul",
        "path",
        "updown",
        "quality",
        "k_column",
        "surf_trans_llight",
        "surf_llight",
        "haul_type",
        "stationid",
        "bottom_depth",
        "orientation"
      ), silent = TRUE, FUN = "estimate_surface"
    ) |> 
      saveRDS(here::here("output", paste0("temp_", region_light, "_interp_huds.rds")))
  
    print("tlu_run_trawllight: Combining residuals and interpolated HUDS, filtering using QA/QC criteria, and recalculating optical depth for interpolated data.")
    readRDS(here::here("output", paste0(region_light, "_tag_residuals.rds"))) |>
      trawllight:::tlu_use_casts() |> 
      merge(readRDS(here::here("output", paste0("temp_", region_light, "_interp_huds.rds")))) |>
      trawllight:::tlu_use_casts2() |>
      trawllight:::tlu_cast_wrapper(id.col = c("vessel", "cruise", "haul", "updown", "path"),
                   FUN = light_proportion) |>
      dplyr::mutate(optical_depth = log(1) - log(light_ratio)) |>
      saveRDS(here::here("output", paste0("temp_", region_light, "_filtered_huds.rds")))
  
  print("tlu_run_trawllight: Generating cast and profile variables.")
  trawllight:::tlu_calc_summary(survey = survey)
  
  if(rm.temp) {
    print("Removing temporary files.")
    file.remove(list.files(path = here::here("output"), pattern = "temp", full.names = TRUE))
  }
  
}