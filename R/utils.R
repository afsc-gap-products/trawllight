#' Open and bind rows of csv to data frame
#' 
#' For use with AFSC/RACE/GAP data structure.
#' 
#' @param directory Path to directory where csv files are stored.
#' @param string Unique pattern used in the csv files.
#' @export

csv_rbind <- function(directory, string) {
  file.list <- grep(pattern = string, x = dir(directory))
  
  if(substr(directory, nchar(directory), nchar(directory)) != "/") {
    if(substr(directory, nchar(directory), nchar(directory)) != "\\") {
      directory <- paste0(directory, "/")
    }
  }
  
  for(i in 1:length(file.list)) {
    if(i == 1) {
      out.df <- read.csv(file = paste0(directory, dir(directory)[file.list[i]]), stringsAsFactors = F)
      out.df$fname <- dir(directory)[file.list[i]]
    } else {
      out.comb <- read.csv(file = paste0(directory, dir(directory)[file.list[i]]), stringsAsFactors = F)
      out.comb$fname <- dir(directory)[file.list[i]]
      out.df <- rbind(out.df, out.comb)
    }
  }
  return(out.df)
}


#' Correct tag time in cases where offsets were incorrect
#'
#' For use processing AOPs from AFSC/RACE/GAP data structure. Make adjustments to correct inconsistencies between tag time and survey time.
#' 
#' @param light.data Data frame with light data
#' @param cast.data Data frame containing case data.
#' @export

time_adjustments <- function(light.data, cast.data) {
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

#' Convert radiometric energy to photon flux density
#' 
#' Convert radiometric energy in watts to micromoles of photons per meter squared per second.
#' 
#' @param x Numeric vector. Energy for a wavelength in units of watts.
#' @param wavelength Numeric vector.Wavelength in nanometers.
#' @return Numeric vector of energy in photon flux density.
#' @export

energy_to_quanta <- function(x, wavelength) {
  return(x/(1e-6*6.02214076*10^23*(3e8/(wavelength*10e-9))*6.63e-34))
  
  
}

#' Photon flux density to radiometric energy 
#' 
#' Convert quantum units of micromoles of photons per meter squared per second to radiometric energy (watts per meter squared per second)
#' 
#' Convert energy to quanta based on wavelength.
#' @param x Numeric vector. Energy for a wavelength in units of watts.
#' @param wavelength Numeric vector.Wavelength in nanometers.
#' @return Numeric vector of energy in radiometric energy in watts.
#' @export

quanta_to_energy <- function(x, wavelength) {
  return(x*1e-6*6.02214076*10^23*(3e8/(wavelength*10e-9))*6.63e-34)
}