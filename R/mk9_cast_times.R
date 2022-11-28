#' Get cast times
#' 
#' This cast.times function retrieves data from Oracle and offset tables to generate a table of event times defining a trawl event. These events are: downcast start, downcast end, on bottom, haul back, off bottom, upcast start and upcast end; they bracket he portion of the trawl haul event that is of interest to us for the use of MK9 irradiance data.  Within this function MK9 light meter depths and times are adjusted to using the info provided in the offsets table so that MK9 depth is calibrated to the survey MBT and MK9 clock drift is adjusted to synch up with the MBT/GPS time that is the basis for the survey clock.
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @return Objects created: MK9.performance.csv - trawl-mounted light meter deployment performance code, CastTimes.csv - catalog of events during trawl haul operation. MK9 performance codes: 0 = good performance, 1 = bad performance, 99 = no data between on and off bottom events.
#' @export

mk9_cast_times <- function(channel = NULL, survey, vessel, cruise){

	channel <- trawllight:::get_connected(channel = channel)

	## cruise year
	year = floor(cruise/100)

	## this code points to the text file named to utilize functions stored there
	## RACE data (either in RACEEDIT or RACEDATA)
	print("getting supporting data...")

	mbt = trawllight:::mk9_get_mbt_data(survey  = survey, 
	                   vessel = vessel, 
	                   cruise = cruise, 
	                   channel = channel)
	
	sgt = trawllight:::mk9_get_sgt_data(survey  = survey, 
	                   vessel = vessel, 
	                   cruise = cruise, 
	                   channel = channel)

	hauls = trawllight:::mk9_get_haul_list(survey = survey, 
	                      vessel = vessel, 
	                      cruise = cruise, 
	                      channel = channel)
	
	## getting values from regression to correct light meter times
	coef <- trawllight:::mk9_get_regr(survey = survey, 
	                             vessel = vessel, 
	                             cruise = cruise, 
	                             channel = channel)
	
	## get list of hauls where need to use light meter for depth data
	no.bt <- hauls[!(hauls$haul %in% unique(mbt$haul)),] 
	no.bt.names <- names(no.bt)

	## Light (MK9) data
	survey.names <- as.data.frame(cbind(cap_name = c("AI","BS","GOA","SLOPE","NBS"), dir_name = c("ai","ebs","goa","slope","nbs")))
	yy <- as.character(substr(cruise,3,4))
	sur <- survey.names$dir_name[survey.names$cap_name == survey]
	light.loc <- here::here("data", "mk9", survey, cruise, vessel)
	files.in.dir <- list.files(path = light.loc, pattern = ".csv")
	file.idx1 <- which(substr(files.in.dir, 1, 4) == "trwl")

	## identify offsets file and read in (there should be just one)
	print("reading time offsets...")

	file.idx2 <- which(substr(files.in.dir, 1, 7) == "offsets")
	# if file.idx2 doesn't return anything then handle that as an error with a note about what to do to get the offs file
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

	## read in trawl-mounted light meter data (there could be more than one file)
	print("reading trawl-mounted MK9 file...")

	for (file in file.idx1) {
		assign(paste("trwl", file, sep = ""),
		read.csv(paste(light.loc, "/", files.in.dir[file], sep= "" ), header = F, skip = 1))
		}

	# empty data set to receive raw MK9 data
	light <- data.frame(V1 = character(), V2 = character(), V3 = numeric(), V4 = numeric(),
	V5 = numeric(), V6 = numeric())

	## put multiple files together
	obj <- objects(pattern = "trwl")
	
	for(i in 1:length(obj)){
		light <- rbind(light, get(obj[i]))
		}

	## make MK9 date and time comparable with sgt and mbt date format
	light$date_time <- as.POSIXct(paste(light$V1, light$V2), format = "%m/%d/%Y %H:%M:%S")

	## 2004 and 2005 MK9 data have a different format than subsequent years; here's code to deal with it
	version_flag <- ncol(light) == 6
	if(cruise <= 200501 & version_flag){
		light <- light[,3:6]
		colnames(light) <- c("ldepth","ltemp","llight","ldate_time")

		## create null database to rbind records into from the loop below
		light.df <- data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), cdepth = numeric(), ctime = character(), 
		ldepth = numeric(), ltemp = numeric(), llight = numeric(), stringsAsFactors = FALSE)
		df.names <- names(light.df)

	}else{
		light <- light[,3:7]
		colnames(light) <- c("ldepth","ltemp","llight","lcond","ldate_time")

		## create null database to rbind records into from the loop below
		light.df <- data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), cdepth = numeric(), ctime = character(), 
		ldepth = numeric(), ltemp = numeric(), llight = numeric(), lcond = numeric(), stringsAsFactors = FALSE)
		df.names <- names(light.df)
	}
	
	light <- trawllight:::mk9_lowpass_filter(x = light, vessel = vessel, cruise = cruise)
	
	light <- trawllight:::mk9_find_offset(light = light, 
	                                      mbt = mbt, 
	                                      try.offsets = c(seq(-8,8,0.5),0.25,-0.25), 
	                                      results.file = paste0(light.loc, "/offset_step1_log.txt"))

	print("assigning hauls to MK9 data...")
	for (h in sort(haullist$haul)){
		## subset data by haul
		haultime <- haullist[haullist$haul == h,c("starttime","haul.end","a","b","offset")]
		## get smaller subset of data between start and end times
		light.haul <- with(light, subset(light, light$ldate_time >= haultime$starttime & light$ldate_time <= haultime$haul.end))
		 
		if(length(haultime[,1]) == 0 | length(light.haul[,1]) == 0) next

		L <- length(light.haul[,1])
		vch <- as.data.frame(cbind(rep(vessel, L), rep(cruise, L), rep(h, L)))
		c.depth <- ((haultime$b * light.haul$ldepth) + haultime$a)
		c.time <- (light.haul$ldate_time + haultime$offset)

		if(cruise <= 200501 & version_flag){
			c.light.dat <- cbind(vch[1], vch[2], vch[3], c.depth, c.time, light.haul$ldepth, light.haul$ltemp, light.haul$llight)
			light.df <- rbind(c.light.dat, light.df)
		}else{
			c.light.dat <- cbind(vch[1], vch[2], vch[3], c.depth, c.time, light.haul$ldepth, light.haul$ltemp, light.haul$llight, light.haul$lcond)
			light.df <- rbind(c.light.dat, light.df)
			}

		}

	names(light.df) <- df.names

	## substituting MK9 data to determine start and end times when MBT data are missing
	print("substituting MK9 depth data where MBT data are missing...")
	if(length(no.bt[,1]) == 0){
		mbt <- mbt
	}else{
		## glue onto MBT the MK9 depth/time data for hauls with no BT data
		for(h in sort(no.bt$haul)){
			mk9.dat <- light.df[light.df$haul == h, c("vessel","cruise","haul","ctime","cdepth")]
			names(mk9.dat) <- names(mbt)
			mbt <- rbind(mbt, mk9.dat)
			}
		}

	## get trawl mounted light meter performance 
	bottom.depth <- hauls[,c("vessel","cruise","haul","bottom_depth")]
	ob.time <- sgt[sgt$time_flag == 3, c("vessel","cruise","haul","date_time")]
	names(ob.time) <- c("vessel","cruise","haul","ob_time")
	fb.time <- sgt[sgt$time_flag == 7, c("vessel","cruise","haul","date_time")]
	names(fb.time) <- c("vessel","cruise","haul","fb_time")
	ob.int <- merge(ob.time, fb.time)

	## create null data frame to be populated in loop below
	perf.df = data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), datum_code = numeric(), stringsAsFactors = FALSE)
	perf.names <- names(perf.df)

	print("detecting trawl-mounted MK9 performance...")
	for(h in sort(ob.int$haul)){

		light.sub <- light.df[light.df$haul == h, ]
		ob.sub <- ob.int[ob.int$haul == h, ]
		ob.light <- with(light.sub, subset(light.sub, light.sub$ctime >= ob.sub$ob_time & light.sub$ctime <= ob.sub$fb_time)) 
		median.light.depth <- median(ob.light$cdepth, na.rm = T)

		## if there are no MK9 data between start and end then performance = -99, otherwise if the MK9 depth is <= 5 
		## then light meter performance is 1, if the MK9 depth is > 5 then it is considered 'good' deployment (performance = 0)
		ifelse(is.na(median.light.depth), perf.df <- rbind(perf.df, c(ob.sub$vessel, ob.sub$cruise, h, -99)),
			(if(!is.na(median.light.depth) <= 5){
			perf = 1
			perf.df = rbind(perf.df, c(ob.sub$vessel, ob.sub$cruise, h, perf))
		}else{
			perf = 0
			perf.df = rbind(perf.df, c(ob.sub$vessel, ob.sub$cruise, h, perf))
			}))
		}

	names(perf.df) <- perf.names

	write.csv(perf.df, paste(light.loc, "/MK9.performance.csv", sep = ""), row.names = F)

	## create CastTimes table as null data frame
	CastTimes <- data.frame(vessel = numeric(), cruise = numeric(), haul = numeric(), downcast_start = character(),
	downcast_end = character(), ob = character(), fb = character(), upcast_start = character(), upcast_end = character(),
	stringsAsFactors = FALSE)
	cast.names <- names(CastTimes)

	## eliminate hauls from ob.int where MK9 performance != 0
	mk9.good.hauls <- ob.int[ob.int$haul %in% (perf.df[perf.df$datum_code == 0,"haul"]), ]

	print("creating CastTimes table...")
	for(h in sort(mk9.good.hauls$haul)){

		mbt.sub <- mbt[mbt$haul == h, c("date_time","depth")]
		ob.sub <- ob.int[ob.int$haul == h, ]
		light.sub <- light.df[light.df$haul == h, ]
		if(length(light.sub[,1]) == 0) {next}

		## if there are no MBT or MK9 data then skip to next haul number
		if(length(mbt.sub[,1]) == 0 & length(light.sub[,1]) == 0) next

		## get downcast start  ##
		dn.cast <- with(mbt.sub, subset(mbt.sub, mbt.sub$date_time < ob.sub$ob_time))
		## if there is MK9 data but not MBT data, then dn.cast data gets MK9 data
		if(length(dn.cast[,1]) == 0){
			dn.cast <- with(light.sub, subset(light.sub, light.sub$ctime < ob.sub$ob_time))
			mbt.down.sub <- dn.cast[,c("ctime","cdepth")]
			names(mbt.down.sub) <- c("date_time","depth")
		}else{
			mbt.down.sub <- dn.cast
			}
		## order data by decreasing time 
		down.order <- mbt.down.sub[order(mbt.down.sub$date_time, decreasing = TRUE), ]
		## get row number where conditions are met
		downcast.idx <- which(down.order$depth >= 1 & down.order$depth <= 2)[1]
		## capture an event time between 0-5 m if there isn't one between 1-2
		downcast.idx <- ifelse(is.na(downcast.idx), which(down.order$depth >= 0 & down.order$depth <= 5)[1], downcast.idx)
		## if there is no observation between 1-2 m or 0-5 m, then take the min time of the downcast as the start
		if(is.na(down.order$date_time[downcast.idx])){
			downcast.start <- down.order$date_time[down.order$date_time == min(down.order$date_time) & down.order$depth >= 0]
		}else{
			downcast.start <- down.order$date_time[downcast.idx]
			}
		## downcast ends 1 second before on bottom time
		downcast.end <- ob.sub$ob_time - 1

		## get upcast end  ##
		up.cast <- with(mbt.sub, subset(mbt.sub, mbt.sub$date_time > ob.sub$fb_time))
		## if there is MK9 data but not MBT data, then dn.cast data gets MK9 data
		if(length(up.cast[,1]) == 0){
			up.cast <- with(light.sub, subset(light.sub, light.sub$ctime > ob.sub$fb_time))
			mbt.up.sub <- up.cast[,c("ctime","cdepth")]
			names(mbt.up.sub) <- c("date_time","depth")
		}else{
			mbt.up.sub <- up.cast
			}
		## order data by increasing time
		up.order <- mbt.up.sub[order(mbt.up.sub$date_time), ]
		## get row number where conditions are met
		upcast.idx <- which(up.order$depth >= 1 & up.order$depth <= 2)[1]
		## capture an event time between 0-5 m if there isn't one between 1-2
		upcast.idx <- ifelse(is.na(upcast.idx), which(up.order$depth >= 0 & up.order$depth <= 5)[1], upcast.idx)
		## if there is no observation between 1-2 m or 0-5 m, then take the max time of the upcast as the end
		if(is.na(up.order$date_time[upcast.idx])){
			upcast.end <- up.order$date_time[up.order$date_time == max(up.order$date_time) & up.order$depth >= 0]
		}else{
			upcast.end <- up.order$date_time[upcast.idx]
			}
		## upcast starts 1 second after off bottom time
		upcast.start <-  ob.sub$fb_time + 1

		casts <- cbind(ob.sub$vessel, ob.sub$cruise, ob.sub$haul, as.character(downcast.start), 
		as.character(downcast.end), as.character(ob.sub$ob_time), as.character(ob.sub$fb_time),
		as.character(upcast.start), as.character(upcast.end))
		CastTimes <- rbind(CastTimes, casts)

		}

	names(CastTimes) <- cast.names

	write.csv(CastTimes, paste(light.loc, "/CastTimes.csv", sep = ""), row.names = F)

	}
