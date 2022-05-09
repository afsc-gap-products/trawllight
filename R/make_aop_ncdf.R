#' Make netCDF files
#' 
#' Bundle apparent optical properties data product into a netCDF file for the specified year and region.
#' 
#' @param region Survey region as a character vector ("BS" = EBS and NBS, "GOA" = Gulf of Alaska, "AI" = Aleutian Islands)
#' @param year Year or years for which data product should be generated.
#' @param global_attributes List of global attributes that is passed to gapctd::df_to_netcdf(global_attributes).
#' @export

make_aop_ncdf <- function(region, year, global_attributes = list(title = "Apparent Optical Properties from AFSC 2021 EBS Shelf and NBS Bottom Trawl Surveys",
                                                                 references = "Rohan, S.K., Kotwicki, S., Kearney, K.A., Schulien, J.A., Laman, E.A., Cokelet, E.D., Beauchamp, D.A., Britt, L.L., Aydin, K.Y., & Zador, S.G. (2021). Using bottom trawls to monitor subsurface water clarity in marine ecosystems. Progress in Oceanography, 194, 102554. https://doi.org/10.1016/j.pocean.2021.102554",
                                                                 id = "https://doi.org/10.5281/zenodo.3688864",
                                                                 cdm_data_type = "Point",
                                                                 cruise = "2021 Eastern Bering Sea Continental Shelf and Northern Bering Sea Bottom-Trawl Surveys",
                                                                 institution = "NOAA Alaska Fisheries Science Center",
                                                                 contributor_name = "Rebecca Haehn, Ned Laman",
                                                                 creator_name = "Sean Rohan",
                                                                 creator_institution = "NOAA Alaska Fisheries Science Center",
                                                                 creator_email = "sean.rohan@noaa.gov",
                                                                 publisher = "NOAA Alaska Fisheries Science Center",
                                                                 publisher_type = "institution",
                                                                 publisher_url = "https://www.fisheries.noaa.gov/about/alaska-fisheries-science-center",
                                                                 geospatial_bounds_crs = "EPSG:4326",
                                                                 license = "http://www.usa.gov/publicdomain/label/1.0/",
                                                                 instrument = "Archival Tag",
                                                                 Conventions = "CF-1.8",
                                                                 standard_name_vocabulary = "CF Standard Name Table v79",
                                                                 source = paste0("Data processed using trawllight ", packageVersion(pkg = "trawllight")))) {
  
  req_attributes <- c("title", "references", "id", "cdm_data_type", "cruise", "institution", "contributor_name",
                      "creator_name", "creator_institution", "creator_email","publisher", "publisher_type", 
                      "publisher_url", "geospatial_bounds_crs", "license", "metadata_link", "instrument", "standard_name_vocabulary", "Conventions", "source")
  
  sel_year <- year
  
  sel_region <- c("ebs", "ebs", "goa", "ai")[match(tolower(region), c("bs", "ebs", "goa", "ai"))]
  
  profile_vars <- readRDS(file = here::here("output", paste0(sel_region, "_final_profile_vars.rds")))
  bottom_vars <- readRDS(file = here::here("output", paste0(sel_region, "_mean_bottom.rds")))
  station_vars <- readRDS(file = here::here("output", paste0(sel_region, "_final_stn_vars.rds")))
  
  profile_vars$year <- floor(profile_vars$cruise/100)
  bottom_vars$year <- floor(bottom_vars$cruise/100)
  station_vars$year<- floor(station_vars$cruise/100)
  
  
  
  profile_vars <- dplyr::filter(profile_vars, year %in% sel_year) |>
    dplyr::select(-longitude, -latitude)
  bottom_vars  <- dplyr::filter(bottom_vars , year %in% sel_year) |>
    dplyr::select(vessel, cruise, haul, tag_bottom_light, year)
  station_vars <- dplyr::filter(station_vars, year %in% sel_year)
  
  unique_profiles <- dplyr::select(profile_vars,
                                   vessel,
                                   cruise, 
                                   haul,
                                   year,
                                   updown) |>
    unique()
  
  interp_profiles <- data.frame()
  
  message(paste0("Making binned diffuse attenuation profiles."))
  
  for(ii in 1:nrow(unique_profiles)) {
    
    profile_out <- unique_profiles[ii,]
    
    sel_profile <- profile_vars |>
      dplyr::inner_join(profile_out,
                        by = c("vessel", "cruise", "haul", "updown", "year")) |>
      dplyr::arrange(depth)
    
    drange <- floor(min(sel_profile$depth)):ceiling(max(sel_profile$depth))
    
    kdz_out <- data.frame(depth = drange,
                          
                          kdz = oce::oce.approx(x = sel_profile$depth, 
                                                y = sel_profile$kdz, 
                                                xout = drange, 
                                                method = "unesco"))
    
    if(is.na(kdz_out$kdz[kdz_out$depth == max(kdz_out$depth)])) {
      kdz_out$kdz[kdz_out$depth == max(kdz_out$depth)] <- sel_profile$kdz[sel_profile$depth == max(sel_profile$depth)]
    }
    
    profile_out <- dplyr::full_join(profile_out, kdz_out, by = character())
    
    interp_profiles <- dplyr::bind_rows(interp_profiles,
                                        profile_out)
    
  }
  
  
  all_profiles <- interp_profiles |>
    dplyr::inner_join(station_vars, by = c("vessel", "cruise", "haul", "year", "updown")) |>
    dplyr::inner_join(bottom_vars, by = c("vessel", "cruise", "haul", "year")) |>
    dplyr::mutate(surface_time = surface_time) |>
    dplyr::rename(haul_depth = bottom_depth,
                  sea_floor_downwelling_irradiance = trans_llight,
                  sea_floor_archival_tag_irradiance = tag_bottom_light,
                  surface_downwelling_irradiance = surf_trans_llight) |>
    dplyr::mutate(quality_flag = ifelse(estimate_ref, 1, 0))
  
  # Define temporal coverage
  time_coverage <- paste0(as.character(range(lubridate::with_tz(all_profiles$surface_time, tz = "UTC"))), " UTC")
  
  all_profiles$surface_time <- as.character(all_profiles$surface_time)
  
  
  # Define spatial extent of data set using WKT polygon
  geospatial_bounds <- cbind(
    c(
      min(all_profiles$latitude),
      min(all_profiles$latitude),
      max(all_profiles$latitude),
      max(all_profiles$latitude),
      min(all_profiles$latitude)),
    c(
      min(all_profiles$longitude), 
      max(all_profiles$longitude),
      max(all_profiles$longitude),
      min(all_profiles$longitude),
      min(all_profiles$longitude))
  )
  
  geospatial_bounds <- paste0("POLYGON ((",
                              paste(apply(X = geospatial_bounds, MARGIN = 1, FUN = paste, collapse = " "), collapse = ", "),
                              "))")
  
  g_attributes <- list(references = global_attributes$references,
                       id = global_attributes$id,
                       cruise =  global_attributes$cruise,
                       institution = global_attributes$institution,
                       contributor_name = global_attributes$contributor_name,
                       creator_name = global_attributes$creator_name,
                       creator_institution = global_attributes$creator_institution,
                       creator_email = global_attributes$creator_email,
                       publisher = global_attributes$publisher,
                       publisher_type = global_attributes$publisher_type,
                       publisher_url = global_attributes$publisher_url,
                       geospatial_bounds = geospatial_bounds,
                       geospatial_bounds_crs = global_attributes$geospatial_bounds_crs,
                       license = global_attributes$license,
                       date_created = as.character(Sys.Date()),
                       instrument = global_attributes$instrument,
                       Conventions = global_attributes$Conventions,
                       standard_name_vocabulary = global_attributes$standard_name_vocabulary,
                       cdm_data_type = global_attributes$cdm_data_type,
                       time_coverage_start = time_coverage[1],
                       time_coverage_end = time_coverage[2],
                       source = global_attributes$source)
  
  
  gapctd::df_to_ncdf(x = all_profiles,
                     output_filename = paste0("AOPLIGHT_", region, paste0(unique(range(c(year))), collapse = "_"), ".nc"),
                     dim_names_2d = c("surface_time", "longitude", "latitude"),
                     dim_long_names_2d = NULL,
                     dim_units_2d = c("time", "degrees_north", "degrees_east"),
                     var_names_2d = c("haul_depth", "near_bottom_optical_depth", "Z10", "Z1", "sea_floor_downwelling_irradiance", "sea_floor_archival_tag_irradiance", "surface_downwelling_irradiance", "quality_flag"),
                     var_long_names_2d = c("Mean towed depth of instrument during haul", "near bottom optical depth", "Depth where optical depth equals 2.3", "Depth where optical depth equals 4.6", "Downwelling irradiance at the sea floor", "Downwelling irradiance at the sea floor in archival tag units", "Downwelling irradiance at the sea surface", "Data Quality Assurance Flag"),
                     var_flag_values_2d = list("quality_flag" = c(0, 1)),
                     var_flag_meanings_2d = list("quality_flag" = c("Good quality", "Surface irradiance at 1 m missing or omitted, reference depth estimated for optical depth calculations.")),
                     var_units_2d = c("m", "1", "m", "m", "micromol m-2 s-1", "1", "micromol m-2 s-1", "1"),
                     dim_names_3d = "depth",
                     dim_long_names_3d = "Depth",
                     dim_positive_3d = list("depth" = "down"),
                     dim_units_3d = "m",
                     dim_sort_3d = TRUE,
                     var_names_3d = c("kdz"),
                     var_long_names_3d = c("Downwelling diffuse attenuation coefficient"),
                     var_units_3d = c("m-1"),
                     var_flag_values_3d = NULL,
                     var_flag_meanings_3d = NULL,
                     instrument_attributes = c("make_model"),
                     instrument_values = list(make_model = "Wildlife Computers TDR-MK9 Archival Tag"),
                     global_attributes = g_attributes)
  
}
