#' Create a database connection using RODBC
#'
#' A function that accepts a data source name, username, and password to establish returns an Oracle DataBase Connection (ODBC) as an RODBC class in R.
#'
#' @param schema Data source name (DSN) as a character vector.
#' @return An RODBC class ODBC connection.
#' @noRd

get_connected <- function(channel = NULL, schema = NA){
  if(is.null(channel)) {
  (echo = FALSE)
  if(is.na(schema)) {
    schema <- getPass::getPass(msg = "Enter ORACLE schema: ")
  }
  username <- getPass::getPass(msg = "Enter your ORACLE Username: ")
  password <- getPass::getPass(msg = "Enter your ORACLE Password: ")
  channel  <- RODBC::odbcConnect(dsn = paste(schema),
                                 uid = paste(username),
                                 pwd = paste(password),
                                 believeNRows = FALSE)
  }
  return(channel)
}



#' Find and apply time offset to Mk9 light data
#'
#' @param light Data frame containing Mk9 data. Must include columns: ldate_time (POSIXct), ldepth (numeric)
#' @param mbt Data frame containing MBT data. Must include columns: date_time (POSIXct), depth (numeric)
#' @param try.offsets A vector of offsets to try. Default is seq(-8,8,0.5)
#' @param results.file Character vector specifying the filepath where information about the offset and correlation between Mk9 and MBT depths are stored.
#' @return The input light data frame with date_time adjusted according to the offset.
#' @noRd

mk9_find_offset <- function(light, mbt, try.offsets = seq(-8,8,0.5), results.file = NULL) {
  
  if(!(("ldate_time" %in% names(light)) & ("ldepth" %in% names(light)))) {
    stop("find_mk9_offset: Columns named ldate_time and/or ldepth are missing from the light argument.")
  }
  
  if(!(("date_time" %in% names(mbt)) & ("depth" %in% names(mbt)))) {
    stop("find_mk9_offset: Columns named date_time and/or depth are missing from the mbt argument.")
  }
  # Initilize vector to store correlations from different offsets
  try.cor <- vector(length = length(try.offsets))
  
  # Loop through offsets
  for(i in 1:length(try.offsets)) {
    # Create offset to try
    light$ldate_time_offset <- light$ldate_time + try.offsets[i]*3600
    offset.df <- dplyr::inner_join(light, mbt, by = c("ldate_time_offset" = "date_time"))
    cor_i <- try(cor(offset.df$ldepth, offset.df$depth, use = "complete.obs"), silent = TRUE)
    
    if(class(cor_i) == "try-error") {
      cor_i <- 0
    }
    
    try.cor[i] <- cor_i
  }
  
  # Remove try column
  light <- light[,-which(colnames(light) == "ldate_time_offset")]
  
  # Transform based on the best offset
  
  
  if(max(try.cor) > 0.5) {
    best_offset <- try.offsets[which.max(try.cor)]
    light$ldate_time <- light$ldate_time + try.offsets[which.max(try.cor)]*3600
  } else {
    best_offset <- 0
  }

  print(paste0("Offset for Mk9 is " , best_offset, " hrs, with correlation between Mk9 and MBT depth of ", try.cor[try.offsets == best_offset], "."))
  
  # Write offset and correlation to .txt file
  if(!is.null(results.file)) {
    fconn <- file(results.file)
    writeLines(c(best_offset,
                 as.character(Sys.Date()),
                 paste0("Offset: " , best_offset, " hrs"),
                 paste0("Corr: ", try.cor[which.max(try.cor)])),fconn)
    close(fconn)
  }
  
  return(light)
  
}



#' Retrieve haul metadata from RACEBASE.HAUL
#' 
#' This function runs a SQLPLUS query against RACEBASE.HAUL returning Vessel, Cruise, Haul, Performance, Haul_type, Wire_length, and Bottom_depth for a specific Region (survey), vessel and cruise.
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @return Returns data frame to object that includes haul metadata for the specified survey, vessel, cruise combination.
#' @noRd


mk9_get_haul_list <- function(channel = NA, survey, vessel, cruise){
  
  # load library if not supplied in function call above
  
  channel <- trawllight:::get_connected(channel = channel)
  
  ## cruise year
  year = floor(cruise/100)
  
  ## Converting RACE_DATA survey definitions back to RACEBase REGIONs
  ## All incoming "survey" designations are equivalent to RACEBase REGIONs except for the Slope
  if(survey == "SLOPE" | survey == "NBS") survey <- "BS"
  
  ## data vintage not tested here because only grainy data streams appear in multiple Oracle schema by date of collection
  ## RACEBASE.HAUL does not utilize RACE_DATA survey.name so translate incoming "survey" to RACEBase REGION
  hauls = RODBC::sqlQuery(channel, paste("select vessel,cruise,haul,performance,haul_type,wire_length,bottom_depth from racebase.haul 
		where vessel = ",vessel," and cruise = ",cruise," and region = '",survey,"' order by haul", sep = ""), 
                          errors = TRUE,  rows_at_time = 1)
  
  names(hauls) <- casefold(names(hauls))
  
  return(hauls)
  
}



#' Retrieve MBT data from RACE_DATA
#' 
#' Retrieves MBT data from the selected survey, vessel, and cruise combination.
#'
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @noRd

mk9_get_mbt_data <- function(survey, vessel, cruise, channel = NULL){
  
  channel <- trawllight:::get_connected(channel = channel)
  survey <- toupper(survey)
  
  ## cruise year
  year = floor(cruise/100)
  
  ## RACE_DATA survey definitions
  if(survey == "GOA") survey.name <- "Gulf of Alaska Bottom Trawl Survey"
  if(survey == "AI") survey.name <- "Aleutian Islands Bottom Trawl Survey"
  if(survey == "BS" | survey == "EBS") survey.name <- "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"
  if(survey == "SLOPE") survey.name <- "Eastern Bering Sea Slope Bottom Trawl Survey"
  if(survey == "NBS") survey.name <- "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"  #"Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey - Triennial Extension"
  
  if(year == 2004){
    
    mbt = RODBC::sqlQuery(channel, paste("select vessel, cruise, haul, date_time, depth
			from race_edit.rb2_btd where vessel = ",vessel,"  and cruise = ",cruise," 
			and datum_code <> 0 order by haul, date_time"), errors = TRUE,  rows_at_time = 1)
    
  }else{
    
    # get data from RACE_DATA edit tables
    bt.data1 <- RODBC::sqlQuery(channel, query = paste("select a.vessel_id vessel, a.cruise cruise, 
			b.haul haul, d.edit_date_time date_time, d.edit_depth depth 
			from race_data.cruises a, race_data.edit_hauls b, race_data.edit_bathythermic_headers c, 
			race_data.edit_bathythermics d, race_data.datum_codes e, race_data.surveys f, 
			race_data.survey_definitions g where g.survey_definition_id = f.survey_definition_id 
			and year = ", year, " and g.survey_name = '", survey.name, "' and f.survey_id = a.survey_id 
			and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and a.vessel_id = ", vessel,
                                                       " and d.bathythermic_header_id = c.bathythermic_header_id and e.datum_code = d.datum_code 
			and e.use_in_analysis = 'Y' order by vessel, haul", sep = ""))
    
    # get data from RACE_DATA final tables
    bt.data2 <- RODBC::sqlQuery(channel, query = paste("select a.vessel_id vessel, a.cruise cruise, b.haul haul, 
			d.edit_date_time date_time, d.edit_depth depth 
			from race_data.cruises a, race_data.hauls b, race_data.bathythermic_headers c, 
			race_data.bathythermics d, race_data.datum_codes e, race_data.surveys f, 
			race_data.survey_definitions g where g.survey_definition_id = f.survey_definition_id 
			and year = ", year, " and g.survey_name = '", survey.name, "' and f.survey_id = a.survey_id 
			and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and a.vessel_id = ", vessel,
                                                       " and d.bathythermic_header_id = c.bathythermic_header_id and e.datum_code = d.datum_code 
			and e.use_in_analysis = 'Y' order by vessel, haul", sep = ""))
    
    # test if data come from the edit or final tables
    if(length(bt.data2[,1]) != 0) {
      mbt <- bt.data2
    }else{
      mbt <- bt.data1
    }
    
  }
  
  names(mbt) <- casefold(names(mbt))
  
  return(mbt)
  
}



#' Retrieve MBT data from RACE_DATA
#' 
#' Retrieves MBT data from the selected survey, vessel, and cruise combination.
#'
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @noRd

mk9_get_mbt_data <- function(survey, vessel, cruise, channel = NULL){
  
  channel <- trawllight:::get_connected(channel = channel)
  survey <- toupper(survey)
  
  ## cruise year
  year = floor(cruise/100)
  
  ## RACE_DATA survey definitions
  if(survey == "GOA") survey.name <- "Gulf of Alaska Bottom Trawl Survey"
  if(survey == "AI") survey.name <- "Aleutian Islands Bottom Trawl Survey"
  if(survey == "BS" | survey == "EBS") survey.name <- "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"
  if(survey == "SLOPE") survey.name <- "Eastern Bering Sea Slope Bottom Trawl Survey"
  if(survey == "NBS") survey.name <- "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"  #"Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey - Triennial Extension"
  
  if(year == 2004){
    
    mbt = RODBC::sqlQuery(channel, paste("select vessel, cruise, haul, date_time, depth
			from race_edit.rb2_btd where vessel = ",vessel,"  and cruise = ",cruise," 
			and datum_code <> 0 order by haul, date_time"), errors = TRUE,  rows_at_time = 1)
    
  }else{
    
    # get data from RACE_DATA edit tables
    bt.data1 <- RODBC::sqlQuery(channel, query = paste("select a.vessel_id vessel, a.cruise cruise, 
			b.haul haul, d.edit_date_time date_time, d.edit_depth depth 
			from race_data.cruises a, race_data.edit_hauls b, race_data.edit_bathythermic_headers c, 
			race_data.edit_bathythermics d, race_data.datum_codes e, race_data.surveys f, 
			race_data.survey_definitions g where g.survey_definition_id = f.survey_definition_id 
			and year = ", year, " and g.survey_name = '", survey.name, "' and f.survey_id = a.survey_id 
			and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and a.vessel_id = ", vessel,
                                                       " and d.bathythermic_header_id = c.bathythermic_header_id and e.datum_code = d.datum_code 
			and e.use_in_analysis = 'Y' order by vessel, haul", sep = ""))
    
    # get data from RACE_DATA final tables
    bt.data2 <- RODBC::sqlQuery(channel, query = paste("select a.vessel_id vessel, a.cruise cruise, b.haul haul, 
			d.edit_date_time date_time, d.edit_depth depth 
			from race_data.cruises a, race_data.hauls b, race_data.bathythermic_headers c, 
			race_data.bathythermics d, race_data.datum_codes e, race_data.surveys f, 
			race_data.survey_definitions g where g.survey_definition_id = f.survey_definition_id 
			and year = ", year, " and g.survey_name = '", survey.name, "' and f.survey_id = a.survey_id 
			and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and a.vessel_id = ", vessel,
                                                       " and d.bathythermic_header_id = c.bathythermic_header_id and e.datum_code = d.datum_code 
			and e.use_in_analysis = 'Y' order by vessel, haul", sep = ""))
    
    # test if data come from the edit or final tables
    if(length(bt.data2[,1]) != 0) {
      mbt <- bt.data2
    }else{
      mbt <- bt.data1
    }
    
  }
  
  names(mbt) <- casefold(names(mbt))
  
  return(mbt)
  
}



#' Get SGT data
#' 
#' Retrieves SGT data
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @noRd

mk9_get_sgt_data <- function(survey, vessel, cruise, channel = NULL){
  
  channel <- trawllight:::get_connected(channel = channel)
  
  # make sure survey in uppercase
  survey <- toupper(survey)
  
  ## cruise year
  year = floor(cruise/100)
  
  ## RACE_DATA survey definitions
  if(survey == "GOA") survey.name <- "Gulf of Alaska Bottom Trawl Survey"
  if(survey == "AI") survey.name <- "Aleutian Islands Bottom Trawl Survey"
  if(survey == "BS" | survey == "EBS") survey.name <- "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"
  if(survey == "SLOPE") survey.name <- "Eastern Bering Sea Slope Bottom Trawl Survey"
  if(survey == "NBS") survey.name <-  "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension" #"Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey - Triennial Extension"
  
  
  print(paste0("Retrieving data for ", cruise, " ", survey.name, " vessel ", vessel))
  
  if(year == 2004){
    
    sgt = RODBC::sqlQuery(channel, paste("select * from race_edit.rb2_sgt where vessel = ",vessel,"
			and cruise = ",cruise," and time_flag in(1,3,5,6,7) order by date_time"),
                          errors = TRUE,  rows_at_time = 1)
    
  }else{
    
    # get data from RACE_DATA edit tables
    start.times1 <- RODBC::sqlQuery(channel, query = paste("select a.vessel_id vessel, a.cruise cruise, b.haul haul, 
			c.edit_date_time date_time, c.event_type_id time_flag, c.edit_latitude latitude, c.edit_longitude longitude 
			from race_data.cruises a, race_data.edit_hauls b, race_data.edit_events c, race_data.surveys d, 
			race_data.survey_definitions e, race_data.datum_codes f where e.survey_definition_id = d.survey_definition_id 
			and d.survey_id = a.survey_id and year = ", year, " and e.survey_name = '", survey.name, "' and 
			b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and a.vessel_id = ", vessel, " and 
			c.event_type_id in (1,3,5,6,7,15) and c.datum_code = f.datum_code 
			and f.use_in_analysis = 'Y' order by c.edit_date_time", sep = ""))
    
    # get data from RACE_DATA final tables
    start.times2 <- RODBC::sqlQuery(channel, query = paste("select a.vessel_id vessel, a.cruise cruise, b.haul haul, 
			c.edit_date_time date_time, c.event_type_id time_flag, c.latitude, c.longitude from race_data.cruises a, 
			race_data.hauls b, race_data.events c, race_data.surveys d, race_data.survey_definitions e, race_data.datum_codes f
			where e.survey_definition_id = d.survey_definition_id and d.survey_id = a.survey_id and year = ", year, 
                                                           " and e.survey_name = '", survey.name, "' and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and a.vessel_id = ",
                                                           vessel, " and c.event_type_id in (1,3,5,6,7,15) and c.datum_code = f.datum_code and f.use_in_analysis = 'Y' 
			order by c.edit_date_time", sep = ""))
    
    # test if start.times are available in RACE_DATA edit or final tables
    if (length(start.times2[,1]) != 0) {
      sgt <- start.times2
    }else{
      sgt <- start.times1
    }
    
  }
  
  names(sgt) <- casefold(names(sgt))
  
  return(sgt)
}

#' Retrieve regression coefficients for cast times correction
#' 
#' Retrieve offset correction coefficients
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @noRd

mk9_get_regr <- function(survey, vessel, cruise, channel = NULL){
  
  channel <- trawllight:::get_connected(channel = channel)
  
  year <- floor(cruise/100)
  
  ## RACE_DATA survey definitions
  if(survey == "GOA") survey.name <- "Gulf of Alaska Bottom Trawl Survey"
  if(survey == "AI") survey.name <- "Aleutian Islands Bottom Trawl Survey"
  if(survey == "BS") survey.name <- "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"
  if(survey == "SLOPE") survey.name <- "Eastern Bering Sea Slope Bottom Trawl Survey"
  if(survey == "NBS") survey.name <- "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension" #"Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey - Triennial Extension"
  
  ## test data vintage to direct sql query to correct Oracle schema
  if(year == 2004){
    # RB2 population tables in RACE_EDIT
    dat <- RODBC::sqlQuery(channel, query = paste("select a.region, a.vessel, a.cruise, a.haul, (24*60*(b.date_time - c.date_time)) duration,
		wire_out from race_edit.rb2_hpden a, race_edit.rb2_sgt b, race_edit.rb2_sgt c where a.vessel = b.vessel 
		and a.cruise = b.cruise and a.haul = b.haul and b.vessel = c.vessel and b.cruise = c.cruise 
		and b.haul = c.haul and b.time_flag = 8 and c.time_flag = 6 and a.region = '",survey,"' and a.vessel = ",vessel," 
		and a.cruise = ",cruise, sep = ""), errors = TRUE,  rows_at_time = 1)
    
  }else{
    
    # Final tables in RACE_DATA
    dat <- RODBC::sqlQuery(channel, query = paste("select survey_name, vessel, cruise, haul, wire_out, (24*60*(doors_up - hb)) duration
		from (select e.survey_name, a.vessel_id vessel, a.cruise cruise, b.haul haul, b.wire_out, c.edit_date_time hb 
		from race_data.cruises a, race_data.hauls b, race_data.events c, race_data.surveys d, race_data.survey_definitions e, 
		race_data.datum_codes f where e.survey_definition_id = d.survey_definition_id and d.survey_id = a.survey_id 
		and e.survey_name = '",survey.name,"' and year = ",year," and a.vessel_id = ",vessel," and b.cruise_id = a.cruise_id 
		and c.haul_id = b.haul_id and c.event_type_id = 6 and c.datum_code = f.datum_code and f.use_in_analysis = 'Y' 
		and b.performance >= 0 and b.haul_type = 3) right outer join (select e.survey_name, a.vessel_id vessel, a.cruise cruise, 
		b.haul haul, b.wire_out, c.edit_date_time doors_up from race_data.cruises a, race_data.hauls b, race_data.events c, 
		race_data.surveys d, race_data.survey_definitions e, race_data.datum_codes f where e.survey_definition_id = d.survey_definition_id 
		and d.survey_id = a.survey_id and e.survey_name = '",survey.name,"' and year = ",year," and a.vessel_id = ",vessel," 
		and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and c.event_type_id = 8 and c.datum_code = f.datum_code 
		and f.use_in_analysis = 'Y' and b.performance >= 0 and b.haul_type = 3) using (survey_name, vessel,cruise,haul,wire_out) 
		order by survey_name,cruise,vessel,haul", sep = ""), errors = TRUE,  rows_at_time = 1)
    
    if(length(dat[,1]) == 0){
      dat <- RODBC::sqlQuery(channel, query = paste("select survey_name, vessel, cruise, haul, wire_out, (24*60*(doors_up - hb)) duration
			from (select e.survey_name, a.vessel_id vessel, a.cruise cruise, b.haul haul, b.edit_wire_out wire_out, c.edit_date_time hb 
			from race_data.cruises a, race_data.edit_hauls b, race_data.edit_events c, race_data.surveys d, race_data.survey_definitions e, 
			race_data.datum_codes f where e.survey_definition_id = d.survey_definition_id and d.survey_id = a.survey_id 
			and e.survey_name = '",survey.name,"' and year = ",year," and a.vessel_id = ",vessel," and b.cruise_id = a.cruise_id 
			and c.haul_id = b.haul_id and c.event_type_id = 6 and c.datum_code = f.datum_code and f.use_in_analysis = 'Y' 
			and b.performance >= 0 and b.haul_type = 3) right outer join (select e.survey_name, a.vessel_id vessel, a.cruise cruise, 
			b.haul haul, b.edit_wire_out wire_out, c.edit_date_time doors_up from race_data.cruises a, race_data.edit_hauls b,
			race_data.edit_events c, race_data.surveys d, race_data.survey_definitions e, race_data.datum_codes f 
			where e.survey_definition_id = d.survey_definition_id 
			and d.survey_id = a.survey_id and e.survey_name = '",survey.name,"' and year = ",year," and a.vessel_id = ",vessel," 
			and b.cruise_id = a.cruise_id and c.haul_id = b.haul_id and c.event_type_id = 8 and c.datum_code = f.datum_code 
			and f.use_in_analysis = 'Y' and b.performance >= 0 and b.haul_type = 3) using (survey_name, vessel,cruise,haul,wire_out) 
			order by survey_name,cruise,vessel,haul", sep = ""), errors = TRUE,  rows_at_time = 1)
    }
  }
  
  names(dat) <- casefold(names(dat))
  
  l.r <- lm(duration ~ wire_out, data = dat)
  coef <- l.r[1]
  
  return(coef)
  
}


#' Internal. Low pass filter malfunctioning MK9
#' 
#' Lowpass filter to fix malfunctioning 1990030 pressure 
#' 
#' @param x Light data.frame
#' @param vessel Vessel
#' @param cruise
#' @noRd

mk9_lowpass_filter <- function(x, vessel, cruise) {
    
    lp_filter <- function(var, tc, freq_n) {
      
      n_var <- length(var)
      aa <- 1 / (1 + 2 * tc * (1/freq_n))
      bb <- (1 - 2 * tc * (1/freq_n)) * aa
      new_var <- numeric(length = n_var)
      new_var[1] <- var[1]
      
      for(jj in 2:n_var) {
        new_var[jj] <- aa*(var[jj]+var[jj-1]) - bb * new_var[jj-1]
      }
      
      return(new_var)
    }
    
    in_var <- x
    
    pass_1 <- lp_filter(var = in_var,
                        tc = 4,
                        freq_n = 1)
    pass_2 <- lp_filter(var = rev(pass_1),
                        tc = 4,
                        freq_n = 1)
    
    x <- rev(pass_2)
  
  return(x)
}
