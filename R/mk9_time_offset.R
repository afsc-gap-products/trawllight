#' Estimate offset between MK9 and MBT times
#' 
#' This function determines the offset in seconds between the date_times of the MK9 and those of the MBT by aligning depths during the downcasts.
#' 
#' @param channel Oracle connection as an RODBC connection object.
#' @param survey RACE survey code as a character vector.
#' @param vessel RACE vessel code a numeric vector.
#' @param cruise RACE cruise code as a numeric vector.
#' @export

mk9_time_offset <- function(vessel, cruise, survey, channel = NULL) {

  channel <- trawllight:::get_connected(channel = channel)
  
	# get BT and event data from Oracle
	mbt <- trawllight:::mk9_get_mbt_data(survey = survey, 
	                                vessel = vessel, 
	                                cruise = cruise, 
	                                channel = channel)
	sgt <- trawllight:::mk9_get_sgt_data(survey = survey, 
	                                vessel = vessel, 
	                                cruise = cruise, 
	                                channel = channel)
	
	if(!(nrow(mbt) > 0)) {
	  stop("mk9_time_offset: no data in mbt from mk9_get_mbt_data.")
	}
	
	if(!(nrow(sgt) > 0)) {
	  stop("mk9_time_offset: no data in sgt from mk9_get_sgt_data.")
	}

	# translating input parameters into variables to be used in eqns
	survey.names <- as.data.frame(cbind(cap_name=c("AI","BS","GOA","SLOPE","NBS"),
		dir_name = c("ai","ebs","goa","slope","nbs")))
	## the mathematical approach yields a numeric value but does not keep the leading zero
	## for years prior to 2010. Since the only place yy shows up after it is populated is in
	## the paste statement, I think the substr() approach should hold up.
	# yy <- floor(cruise/100) - floor(cruise/10000) * 100
	yy <- substr(floor(cruise/100), 3, 4)
	sur <- survey.names$dir_name[survey.names$cap_name == survey]

	# path to data and storage directory
	light.loc <- here::here("data", "mk9", survey, cruise, vessel)
	
	# identify the MK9 data that came from the trawl-mounted meter
	# this is a naming convention applied when writing the CSV out from HexDecode
	files.in.dir <- list.files(path = light.loc, pattern = ".csv")
	file.idx <- which(substr(files.in.dir, 1, 4) == "trwl")

	# there can be multiple files on a single vessel and these will be numbered sequentially
	# this loop reads them all and assigns them to an internal variable name
	for (file in file.idx) {

		# assign a value to a name in an environment
		assign(paste("trwl", file, sep = ""), read.csv(paste(light.loc, "/", files.in.dir[file], 
			sep = ""), header = F, skip = 1))
			
		}

	# calls the rbind function and the vector of names assigned in the loop above
	# basically appends all of the rows in however many trwl files were translated
	# from the MK9 to the light.loc
	light <- do.call("rbind", c(mget(objects(pattern = "trwl"))))

# TROUBLESHOOT
# When getting error message that light$date_time ... [,3:7]
# typically has to do with having the wrong number of fields extracted
# at the HexDecode step
	
	# format MK9 time to same as found in SGT and MBT
	light$date_time <- as.POSIXct(paste(light$V1, light$V2), format = "%m/%d/%Y %H:%M:%S")
	

	
	print(paste0("mk9_time_offset: Imported ", nrow(light), " rows of data in ", ncol(light), " variable columns."))
	
	#######changed from 7 or 6
	if(ncol(light) == 7) {
	  light <- light[ ,3:7]
	  colnames(light) = c("ldepth", "ltemp", "llight", "lcond", "ldate_time")
	} else {
	  stop("Check number of columns in mk9_corr and create new value in time_offset()")
	}

	# ldepth = MK9 depth
	# ltemp = MK9 temperature
	# llight = MK9 raw photon count number
	# lcond = MK9 conductivity
	# ldate_time = light$date_time from above

	print(paste0("mk9_time_offset: ", "Creating haul list from event files where there is an On Bottom event."))
	# creates haul list from event files where there is an On Bottom event
	haullist <- sgt$haul[which(sgt$time_flag == 3)]
	
	light <- trawllight:::mk9_find_offset(light = light, 
	                                     mbt = mbt, 
	                                     try.offsets = c(seq(-8,8,0.5),0.25,-0.25), 
	                                     results.file = paste0(light.loc, "/offset_step1_log.txt"))

	# create objects to be populated later
	rsqr = as.data.frame(matrix(data = NA, nrow = length(haullist), ncol = 361))
	aa = as.data.frame(matrix(data = NA, nrow = length(haullist), ncol = 361))
	bb = as.data.frame(matrix(data = NA, nrow = length(haullist), ncol = 361))
	a = vector(length = length(haullist))
	b = vector(length = length(haullist))

	# by haul loop to adjust ldepth to MBT DEPTH
	i = 2
	for (i in 1:length(haullist)) {
	
	  sub.mbt = subset(mbt, haul == haullist[i] & depth > 5)
	  sub.sgt = subset(sgt, haul == haullist[i])
	  min.mbt =  min(sub.mbt$date_time)
	  max.mbt =  max(sub.mbt$date_time)
	  sub.light = subset(light, ldate_time > min.mbt - 300 & ldate_time < max.mbt + 300)
	  colnames(sub.light)[5] = "date_time"
	  sub.light2 = sub.light
	  print(i)

	  for (j in 1:361) {
	  
		# sensitivity of fitting routine is 3 minutes
		sub.light2$date_time = sub.light$date_time + (-181 + j)
		if(length(sub.light2$date_time) < 11) {
			rsqr[i,j] = NA
			aa[i,j] = NA
			bb[i,j] = NA
		}else{
			merged = merge(sub.mbt, sub.light2)
			if(length(merged$date_time) < 11) {
				rsqr[i,j] = NA
				aa[i,j] = NA
				bb[i,j] = NA
			}else{
				# aligning light meter ldepth with MBT DEPTH with a linear model
				rsqr[i,j] = summary(lm(depth ~ ldepth, data = merged))$r.squared
				aa[i,j] = lm(depth ~ ldepth, data = merged)$coefficients[1]
				bb[i,j] = lm(depth ~ ldepth, data = merged)$coefficients[2]
				}
			}
		}
	}

	offset = vector(length = length(haullist))
	max.rsqr = vector(length = length(haullist))

	# by haul loop to generate time offsets
	i = 1
	print("mk9_time_offset: Making plots")
	pdf(file = paste0(light.loc, "/plot_offset_step_2.pdf"), onefile = TRUE)
	for (i in 1:length(rsqr[,1])) {
		if(is.na(rsqr[i,180])) {
			offset[i] = NA
			a[i] = NA
			b[i] = NA
		}else{
			offset[i] = which.max(rsqr[i,]) - 181
			a[i] = aa[i,][rsqr[i,] == max(rsqr[i,])]
			b[i] = bb[i,][rsqr[i,] == max(rsqr[i,])]
			max.rsqr[i] = max(rsqr[i,])
			
			if(length(aa[i,][rsqr[i,] == max(rsqr[i,])]) != 1) {
			  warning(paste0("mk9_time_offset: length of 'length(aa[i,][rsqr[i,] == max(rsqr[i,])])' is equal to ", length(aa[i,][rsqr[i,] == max(rsqr[i,])])))
			}
			
			if(length(bb[i,][rsqr[i,] == max(rsqr[i,])]) != 1) {
			  warning(paste0("mk9_time_offset: length of 'length(bb[i,][rsqr[i,] == max(rsqr[i,])])' is equal to ", length(bb[i,][rsqr[i,] == max(rsqr[i,])])))
			}
			
			if(any(is.na(c(a[i], b[i])))) {
			  warning("mk9_time_offset: max.rsqr for ", haullist[i], " is NA")
			}
			
			plot(c(-180:180), rsqr[i,], ylim = c(min(rsqr[i,], na.rm = TRUE), max(rsqr[i,], na.rm = TRUE)), main = paste("Haul No.", haullist[i], sep = " "),
				ylab = "iterative fit (r-sqr)", xlab = "3-minute range for fitting (seconds)")
			}
	}
	dev.off()

	haul.time = sgt$date_time[sgt$time_flag == 3]
	
	pdf(file = paste0(light.loc, "/plot_offset_by_haul.pdf"), onefile = TRUE)
	# plot offsets by haul
	plot(offset ~ haul.time, col = "darkblue", main = "MK9 time offset from MBT \nacross survey", 
		xlab = "time (haul start)", ylab = "offset (seconds)", ylim = c(0.5, 1))
	# use local polynomial smoothing across all hauls to get smoothed trend line for offsets over cruise
	# weights are keyed to max.rsqr, span controls the degree of smoothing
	smooth <- try(predict(loess(offset ~ as.numeric(haul.time), weights = ifelse(max.rsqr > .85, max.rsqr^10, 0), 
		span = 20/length(offset)), seq(min(haul.time), max(haul.time), I(60*30))), silent = TRUE)
	
	if(class(smooth)[1] == "try-error") {
	  smooth <- predict(loess(offset ~ as.numeric(haul.time), weights = ifelse(max.rsqr > .85, max.rsqr^10, 0), 
	                    span = 20/length(offset)), seq(min(haul.time), max(haul.time), I(60*30)))
	}
	lines(seq(min(haul.time), max(haul.time), I(60*30)), smooth, col = 'red')
	dev.off()
	# smooth2 is the loess predicted offset value for each haul
	smooth2 = predict(loess(offset ~ as.numeric(haul.time), weights = ifelse(max.rsqr > .85, max.rsqr^10,0), 
		span = 0.75), haul.time)
	offsets = cbind(vessel, cruise, haul = haullist, offset, max.rsqr, a, b,
		eq100m = a + b * 100, smooth_offset = smooth2)
	
	# CSV contains Vessel, Cruise, Haul, MK9 time offset (seconds), and parameters of linear model
	write.csv(offsets, paste(light.loc, "/offsets.csv", sep = ""), row.names = F)

	}
