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
  
  trawllight:::tlu_prep_haul_data(channel = channel,
                                 survey = survey)
  
  trawllight:::tlu_prep_dir_list(survey = survey,
                                light_data_root = light_data_root)
  
}