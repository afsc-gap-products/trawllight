#' Wrapper function to retrieve cast data and aggregate
#' 
#' Description goes here... writes binned cast data to an rds file after interpolating missing surface values.
#' 
#' @param directory_structure File path to csv file that lists directories to be processed.
#' @param survey Survey as a character vector ("BS", "NBS", "GOA", "AI", or "SLOPE") 
#' @param ... Additional arguments passed to tlu_surface_light
#' @export

tlu_get_surface <- function(directory_structure = NULL, 
                            survey,
                            ...) {
  
  if(is.null(directory_structure)) {
    directory_structure <- read.csv(file = here::here("imports", "directories.csv"))
  }
  
  region_light <- c("ebs", "nbs", "goa", "ai", "slope")[match(survey, c("BS", "NBS", "GOA", "AI", "SLOPE"))]
  out_path <- here::here("output", paste0(region_light, "_surface", ".rds"))
  
  surf <- trawllight:::tlu_process_all_surface(dir.structure = directory_structure[,1], survey = survey, ...)
  
  saveRDS(surf, file = out_path)
  return(surf)
}