#' Write time-corrected data to csv
#' 
#' Use estimated offsets to time-correct MK9 data.
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @param make.plots Should cast plots be written to pdf?
#' @export

mk9_extinction <- function(channel, survey, vessel, cruise, make.plots = TRUE){
  
  channel <- trawllight:::get_connected(channel = channel)

	## establish region names and light data location
	survey.names <- as.data.frame(cbind(cap_name = c("AI","BS","GOA","SLOPE","NBS"), 
	                                    dir_name = c("ai","ebs","goa","slope","nbs")))
	yy <- as.character(substr(cruise,3,4))
	sur <- survey.names$dir_name[survey.names$cap_name == survey]
	light.loc <-  here::here("data", "mk9", survey, cruise, vessel)
	# extinct.loc <- "G:\\RACE_LIGHT\\LightData\\Data\\ExtinctionCoefficients\\"

	## get data
	hauls = trawllight:::mk9_get_haul_list(survey = survey ,
	                      vessel = vessel, 
	                      cruise = cruise, 
	                      channel = channel)
	casts <- read.csv(file = paste(light.loc, "/CastTimes.csv", sep = ""), header = T)
	perf <- read.csv(file = paste(light.loc, "/MK9.performance.csv", sep = ""), header = T)
	sgt <- trawllight:::mk9_get_sgt_data(survey = survey ,
	                    vessel = vessel, 
	                    cruise = cruise, 
	                    channel = channel)
	mbt <- trawllight:::mk9_get_mbt_data(survey = survey, 
	                                     vessel = vessel, 
	                                     cruise = cruise, 
	                                     channel = channel)
	
	## getting values from regression to correct light meter times
	coef <- trawllight:::mk9_get_regr(survey = survey ,
	                 vessel = vessel, 
	                 cruise = cruise, 
	                 channel = channel)

	# get trawl-mounted MK9 good performance casts
	good.perf <- perf[perf$datum_code == 0, c("vessel","cruise","haul")]

	## get Light (MK9) data
	print("reading trawl-mounted MK9 file...")

	files.in.dir <- list.files(path = light.loc, pattern = ".csv")
	file.idx1 <- which(substr(files.in.dir, 1, 4) == "trwl")

	## read in trawl-mounted light meter data (there could be more than one file)
	for (file in file.idx1) {
		assign(paste("trwl", file, sep = ""),
		read.csv(paste(light.loc, "/", files.in.dir[file], sep= "" ), header = F, skip = 1))
		}

	if(cruise <= 201801){
		trwl.dat <- data.frame(V1 = numeric(), 
		                       V2 = numeric(), 
		                       V3 = numeric(), 
		                       V4 = numeric(), 
		                       V5 = numeric(), 
		stringsAsFactors = FALSE)
	}else{
		trwl.dat <- data.frame(V1 = numeric(), 
		                       V2 = numeric(), 
		                       V3 = numeric(), 
		                       V4 = numeric(), 
		                       V5 = numeric(), 
		                       V6 = numeric(), 
		                       stringsAsFactors = FALSE)
		}

	## put multiple files together
	obj <- objects(pattern = "trwl")
	for(i in 1:length(obj)){
		trwl.dat <- rbind(trwl.dat, get(obj[i]))
		}

	## make MK9 date and time comparable with sgt and mbt date format
	trwl.dat$date_time <- as.POSIXct(paste(trwl.dat$V1, trwl.dat$V2), format = "%m/%d/%Y %H:%M:%S")

	light <- trwl.dat

	# empty data set to receive raw MK9 data
	## 2004 and 2005 MK9 data have a different format than subsequent years; here's code to deal with it
	if(ncol(light) == 6){ #cruise <= 201801
		light <- light[,3:6]
		colnames(light) <- c("ldepth","ltemp","llight","ldate_time")

		## create null database to rbind records into from the loop below
		light.df <- data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), cdepth = numeric(), ctime = character(), 
			ldepth = numeric(), ltemp = numeric(), llight = numeric(), stringsAsFactors = FALSE)
			df.names <- names(light.df)
	} else if(ncol(light) == 7){
		light <- light[,3:7]
		colnames(light) <- c("ldepth","ltemp","llight","lcond","ldate_time")

		## create null database to rbind records into from the loop below
		light.df <- data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), cdepth = numeric(), ctime = character(), 
			ldepth = numeric(), ltemp = numeric(), llight = numeric(), lcond = numeric(), stringsAsFactors = FALSE)
			df.names <- names(light.df)
	} else {
		  stop(paste0("extinction: Unexpected number of columns, (", ncol(light), ") in light data csv. Expected 6 or 7."))
	}
	
	light <- trawllight:::mk9_find_offset(light = light, 
	                                      mbt = mbt, 
	                                      try.offsets = c(seq(-8,8,0.5),0.25,-0.25), 
	                                      results.file = paste0(light.loc, "/offset_step1_log.txt"))

	## identify offsets file and read in (there should be just one)
	file.idx2 <- which(substr(files.in.dir, 1, 7) == "offsets")
	assign("offsets", read.csv(paste(light.loc, "/", files.in.dir[file.idx2], sep = ""), header = TRUE))
	# eliminate hauls where an offset time could not be calculated

	offset.names <- names(offsets)
	## substitute smoothed offset times for NAs and for offsets where R-sqr < 0.97
	mixed.offsets <- ifelse(is.na(offsets$offset),round(offsets$smooth_offset,0),
		ifelse(offsets$max.rsqr < 0.97,round(offsets$smooth_offset,0),
		offsets$offset))
	## substitute the mixed_offsets vector for the original offsets column
	offsets <- cbind(offsets[,c("vessel","cruise","haul")], mixed.offsets, offsets[,c("max.rsqr","a","b","eq100m","smooth_offset")])
	names(offsets) <- offset.names
	## eliminate any remaining records where offset was NA; these will be cases where the whole record was NA
	## presumably because of no irradiance data
	offsets <- offsets[!is.na(offsets$offset),]

	## Depending on data vintage start time can be a 1 or a 15
	start_time <- sgt[sgt$time_flag == 1 | sgt$time_flag == 15, c("vessel","cruise","haul","date_time")]
	names(start_time) <- c("vessel","cruise","haul","starttime")
	## haul back time event remains 6
	hb.time <- sgt[sgt$time_flag == 6, c("vessel","cruise","haul","date_time")]
	## merge sgt and hauls to get hauls aligned
	hb.time <- merge(hauls, hb.time)

	## linear regression of duration from haul back to doors up + 50% more
	wire.time <- ceiling(2*60*(coef$coefficients[1] + (coef$coefficients[2]*hb.time$wire_length)))

	############# potential for different length haul lists needs to be handled:  use one master haul list #####################3

	## add predicted haul end to haulback time
	haul.end <- hb.time$date_time + wire.time
	end_time <- cbind(hb.time[, c("vessel","cruise","haul")], haul.end)
	## append start and end times to haul list
	haullist <- merge(offsets, start_time)
	haullist <- merge(haullist, end_time)
	haullist <- merge(haullist, good.perf)

	extinct.df <- data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), extinction = numeric(), longitude = numeric(), 
		latitude = numeric(), stringsAsFactors = FALSE)
		e.names <- names(extinct.df)

	## assign haul numbers, time offsets, and depth corrections to MK9 data to produce corrected data set
	## also compute extinction coefficients for down and upcasts and assign them to global positioning coordinates
	for(h in sort(haullist$haul)) {

		# print(h)
		## subset data by haul
		haultime <- haullist[haullist$haul == h,c("starttime","haul.end","a","b","offset")]
		## get smaller subset of data between start and end times
		light.haul <- with(light, subset(light, light$ldate_time >= haultime$starttime & light$ldate_time <= haultime$haul.end)) 

		if(length(haultime[,1]) == 0 | length(light.haul[,1]) == 0) next

		## correct 
		L <- length(light.haul[,1])
		vch <- as.data.frame(cbind(rep(vessel, L), rep(cruise, L), rep(h, L)))
		c.depth <- ((haultime$b * light.haul$ldepth) + haultime$a)
		c.time <- (light.haul$ldate_time + haultime$offset)

		if(ncol(light) == 4){
			c.light.dat <- cbind(vch[1], vch[2], vch[3], c.depth, c.time, light.haul$ldepth, light.haul$ltemp, light.haul$llight)
			names(c.light.dat) <- df.names
			light.df <- rbind(c.light.dat, light.df)
			names(light.df) <- df.names
		} else if(ncol(light) == 5){
			c.light.dat <- cbind(vch[1], vch[2], vch[3], c.depth, c.time, light.haul$ldepth, light.haul$ltemp, light.haul$llight, 
			light.haul$lcond)
			names(c.light.dat) <- df.names
			light.df <- rbind(c.light.dat, light.df)
			names(light.df) <- df.names
			}

		# get single cast event times
		cast <- casts[casts$haul == h, ]
		# convert character date times to POSIX
		dn.st <- as.POSIXct(cast$downcast_start)
		dn.nd <- as.POSIXct(cast$downcast_end)
		up.st <- as.POSIXct(cast$upcast_start)
		up.nd <- as.POSIXct(cast$upcast_end)
		# extract MK9 casts
		down <- light.df[light.df$haul == h & light.df$ctime >= dn.st & light.df$ctime <= dn.nd,c("llight","cdepth")]
		up <- light.df[light.df$haul == h & light.df$ctime >= up.st & light.df$ctime <= up.nd,c("llight","cdepth")]

		# get global position of OB and FB events
		## note that using other events here run risk of tapping multiple coordinate formats stored in RACE_DATA
		ob.pos <- sgt[sgt$haul == h & sgt$time_flag == 3, c("longitude","latitude")]
		fb.pos <- sgt[sgt$haul == h & sgt$time_flag == 7, c("longitude","latitude")]

		if(length(down[,1]) == 0){
			dn.cast.extinct = NA
		} else{
			# get light intercept from light regressed on depth
			down.coef <- lm(llight ~ cdepth, data = down)
			down.0 <- down.coef$coefficients[1]
			# there could be multiple records at max depth so take average
			ob.depth <- max(down$cdepth)
			ob.light <- mean(down[down$cdepth == ob.depth, "llight"])
			# take natural log of the change in light between surface and bottom over bottom depth
			dn.cast.extinct <- -(log(ob.light/down.0)/ob.depth)
			}

		if(length(up[,1]) == 0){
			up.cast.extinct = NA
		}else{
			up.coef <- lm(llight ~ cdepth, data = up)
			up.0 <- up.coef$coefficients[1]
			fb.depth <- max(up$cdepth)
			fb.light <- mean(up[up$cdepth == fb.depth, "llight"])
			up.cast.extinct <- -(log(fb.light/up.0)/fb.depth)
			}

		r1 <- cbind(unique(vch[1]), unique(vch[2]), unique(vch[3]), extinction = dn.cast.extinct, ob.pos)
		r2 <- cbind(unique(vch[1]), unique(vch[2]), unique(vch[3]), extinction = up.cast.extinct, fb.pos)
		## for rbind to work need to keep all of the field names the same throughout concatenation process
		r1r2 <- rbind(r1,r2); names(r1r2) <- e.names

		extinct.df <- rbind(extinct.df, r1r2); names(extinct.df) <- e.names

		}

	names(extinct.df) <- e.names
	names(light.df) <- df.names

	## place corrected haul-specific light meter data into the vessel's data directory
	write.csv(light.df, paste(light.loc, "/corr_MK9hauls.csv", sep= ""), row.names = F)
	
	if(make.plots) {
	  dat_to_plot <- light.df |>
	    dplyr::select(vessel, cruise, haul, cdepth, ctime, ltemp, llight)
	  
	  dat_to_plot$ctime <- as.POSIXct(dat_to_plot$ctime)
	  
	  head(dat_to_plot)
	  
	  unique_hauls <- dat_to_plot |> 
	    dplyr::select(vessel, cruise, haul) |>
	    unique()
	  
	  nrow(unique_hauls)
	  
	  print(paste0("Writing haul plots to ", light.loc, "/plots_cast_dat.pdf"))
	  pdf(file = paste0(light.loc, "/plots_cast_dat.pdf"), onefile = TRUE)
	  for(ii in 1:nrow(unique_hauls)) {
	    haul_sub <- dat_to_plot |>
	      dplyr::filter(vessel == unique_hauls$vessel[ii], 
	                    cruise == unique_hauls$cruise[ii], 
	                    haul == unique_hauls$haul[ii]) 
	    
	    print(
	      ggplot2::ggplot() +
	        ggplot2::geom_line(data =  haul_sub |> 
	                            tidyr::pivot_longer(cols = c("ltemp", "llight", "cdepth")),
	                          ggplot2::aes(x = ctime, 
	                                      y = value, 
	                                      color = name)) +
	        ggplot2::facet_wrap(~name, scales = "free_y", 
	                           nrow = 3) +
	        ggplot2::ggtitle(paste0("Vessel: ", unique_hauls$vessel[ii], ", Cruise: ", unique_hauls$cruise[ii], ", Haul: ", unique_hauls$haul[ii])) + 
	        ggplot2::theme_bw()
	    )
	  }
	  dev.off()
	}

	}
