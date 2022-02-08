#' Convert WCDAP archive csv to trawllight csv format
#' 
#' Converts an archive CSV file produced by WC-DAP or WC Instrument helper to a CSV format that is compatible with trawllight. The CSV format for trawllight is based on a format that had been produced by Wildlife Computers HexDecoder software prior to 2017
#' 
#' @param dir_path Path to the directory containing the Archive.csv file
#' @param file_name Name of the Archive.csv file. If NULL, searches for the file name in dir_path.
#' @param plot_check Should a plot (.png) of datetime versus depth be written to dir_path?
#' @param type User specified type. Detected automatically if type = NULL (default). Otherwise, either "trwl" for trawl-mounted tag or "deck" for deck-mounted tag.
#' @export

mk9_convert_wcdap_archive <- function(dir_path = NULL,
                              file_name = NULL,
                              plot_check = TRUE,
                              type = NULL) {
  
  if(is.null(dir_path)) {
    dir_path <- here::here()
  }
  
  if(is.null(file_name)) {
    file_name <- list.files(dir_path, 
                            pattern = "Archive.csv", 
                            full.names = TRUE)
  }
  
  file.copy("process_mk9.Rmd", system.file(pacakge = "mk9process"))
  
  for(ii in 1:length(file_name)) {  
    
    tag_id <- stringr::str_split_fixed(file_name[ii], "_", n = 2)[,1]
    
    # Read data and skip headers
    print(paste0("Reading ", here::here(dir_path, file_name[ii])))
    dat <- read.csv(file = here::here(dir_path, file_name[ii]), 
                    skip = 2,
                    header = FALSE)
    
    # Select wet/dry sensor values >0
    dat <- dat[which(dat[,5] > 0),]
    
    # Remove unused columns
    dat <- dat[,-((ncol(dat)-4):ncol(dat))]
    
    # Assign column names
    names(dat) <- c("ldatetime", "Depth", "Temperature", "Light Level")
    
    # Split datetime
    dat <- dplyr::mutate(dat, 
                         Time = stringr::str_split_fixed(ldatetime, " ", 2)[,2],
                         Date = stringr::str_split_fixed(ldatetime, " ", 2)[,2])
    
    if(is.null(type)) {
      if(diff(range(dat$Depth, na.rm = TRUE)) >= 20) {
        type <- "trwl"
      } else {
        type <- "deck"
      }
    } else {
      if(!(type %in% c("trwl", "deck"))) {
        stop(paste0("Chosen type (", type, ") is invalid. Must be 'trwl' or 'deck'"))  
      }
    }
    
    if(plot_check) {
      png(file = here::here(dir_path, paste0("plotdat_", tag_id, ".png")), 
          width = 8, 
          height = 6, 
          units = "in", 
          res = 150)
      print(
        plot(x = as.POSIXct(dat$ldatetime), 
             y = -1*dat$Depth, 
             type = 'l', 
             xlab = "Datetime", 
             ylab = "Depth",
             main = paste0(type, " ", tag_id))
      )
      dev.off()
    }
    
    # Reorder columns for output
    dat <- data.frame(Date = dat$Date,
                      Time = dat$Time,
                      Depth = dat$Depth,
                      Temperature = dat$Temperature,
                      `Light Level` = dat$`Light Level`)
    
    print("Writing output to .csv")
    write.csv(dat, here::here(dir_path, paste0(type, "_", tag_id, ".csv")), row.names = FALSE)
  }
}
