#' Copy MK9 data and setup directory
#' 
#' @param source_path Filepath to the data source
#' @param overwrite Should the local data directory be overwritten?
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @export

mk9_setup_dir <- function(source_path, survey, cruise, vessel, overwrite = FALSE) {
  
  flag <- FALSE
  
  if(!dir.exists(here::here("data", "mk9", survey, cruise, vessel))) {
    flag <- TRUE
  } else if(dir.exists(here::here("data", "mk9", survey, cruise, vessel)) & overwrite) {
    flag <- TRUE
  } else {
    warning(paste0("data directory already exists! No files copied from ", source_path))
  }
  
  if(flag) {
    dir.create(here::here("data", "mk9", survey, cruise, vessel), recursive = TRUE)
    
    file.copy(list.files(source_path, full.names = TRUE), 
              to = here::here("data", "mk9", survey, cruise, vessel))
  }

}