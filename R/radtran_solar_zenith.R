#' Calculate solar zenith angle
#' 
#' Calculate solar zenith angle for a latitude and longitude based on position, date, time, and timezone. From the function fishmethods::astrocalc4r.
#' 
#' @param lat Latitude in decimal degrees.
#' @param lon Longitude in decimal degrees.
#' @param month Month. Integer.
#' @param day Day of month. Integer.
#' @param hour Hour (decimal hour). Numeric.
#' @param year Year. Numeric.
#' @param timezone Local time relative to GMT (Alaska Daylight Time is -8, i.e., GMT-8)
#' @return Returns a numeric vector of solar zenith angle (degrees).
#' @export

radtran_solar_zenith <- function(lat, lon, day, month, year, hour, timezone) {
  deg2rad <- pi/180
  
  is.leap <- function(x) return((((x%%4 == 0) & (x%%100 != 
                                                   0))) | (x%%400 == 0))
  date.list <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 
                 31)
  
  JulianDay <- function(xday, xmonth, xyear) {
    mm <- xmonth
    xmonth[mm <= 2] <- xmonth[mm <= 2] + 12
    xyear[mm <= 2] <- xyear[mm <= 2] - 1
    xa <- floor(xyear/100)
    xb <- 2 - xa + floor(xa/4)
    jd <- floor(365.25 * (xyear + 4716)) + floor(30.6001 * 
                                                   (xmonth + 1)) + xday + xb - 1524.5
    return(jd)
  }
  
  daymonth <- function(mth, yr) {
    day[is.leap(yr)] <- c(31, 29, 31, 30, 31, 30, 31, 31, 
                          30, 31, 30, 31)[mth[is.leap(yr)]]
    day[!is.leap(yr)] <- c(31, 28, 31, 30, 31, 30, 31, 31, 
                           30, 31, 30, 31)[mth[!is.leap(yr)]]
    return(day)
  }
  
  hourtemp <- hour - timezone
  hour <- ifelse(hourtemp > 24, hourtemp - 24, hourtemp)
  change_day <- !(hour == hourtemp)
  dm <- daymonth(month, year)
  daytemp <- day
  daytemp[change_day] <- ifelse((day[change_day] < dm[change_day]), 
                                day[change_day] + 1, 1)
  change_month <- abs(day - daytemp) > 1
  monthtemp <- month
  monthtemp[change_month] <- ifelse(month[change_month] < 
                                      12, month[change_month] + 1, 1)
  change_year <- abs(month - monthtemp) > 1
  
  yeartemp <- year
  yeartemp[change_year] <- year[change_year] + 1
  
  xy <- yeartemp
  xm <- monthtemp
  xd <- daytemp + hourtemp/24
  jd <- JulianDay(xd, xm, xy) * 100/100
  
  jc <- (jd - 2451545)/36525
  xx <- 280.46646 + 36000.76983 * jc + 0.0003032 * jc^2
  gmls <- xx%%360
  xx <- 357.52911 + 35999.05029 * jc - 0.0001537 * jc^2
  gmas <- xx%%360
  
  eeo <- 0.016708634 - 4.2037e-05 * jc - 1.267e-07 * jc^2
  scx <- (1.914602 - 0.004817 * jc - 1.4e-05 * jc^2) * sin(gmas * 
                                                             deg2rad) + (0.019993 - 0.000101 * jc) * sin(2 * gmas * 
                                                                                                           deg2rad) + 0.000289 * sin(3 * gmas * deg2rad)
  Stl <- gmls + scx
  Sta <- gmas + scx
  srv <- 1.000001018 * (1 - eeo^2)/(1 + eeo * cos(Sta * deg2rad))
  omega <- 125.04 - 1934.136 * jc
  lambda <- Stl - 0.00569 - 0.00478 * sin(omega * deg2rad)
  epsilon <- (23 + 26/60 + 21.448/60^2) - (46.815/60^2) * 
    jc - (0.00059/60^2) * jc^2 + (0.001813/60^2) * jc^3
  oblx <- 0.00256 * cos(omega * deg2rad)
  epsilon <- epsilon + oblx
  alpha <- atan2(cos(epsilon * deg2rad) * sin(lambda * deg2rad), 
                 cos(lambda * deg2rad))/deg2rad
  declin <- asin(sin(epsilon * deg2rad) * sin(lambda * deg2rad))/deg2rad
  y <- tan(epsilon * deg2rad/2)^2
  eqtime <- (y * sin(2 * gmls * deg2rad) - 2 * eeo * sin(gmas * 
                                                           deg2rad) + 4 * eeo * y * sin(gmas * deg2rad) * cos(2 * 
                                                                                                                gmls * deg2rad) - y^2 * sin(4 * gmls * deg2rad)/2 - 
               5/4 * eeo^2 * sin(2 * gmas * deg2rad))/deg2rad * 4
  
  tst <- (hourtemp * 60 + eqtime + 4 * lon)%%1440
  tsa <- ifelse(tst < 0, tst/4 + 180, tst/4 - 180)
  
  zenith <- 90 - asin(sin(lat * deg2rad) * sin(declin * deg2rad) + 
                        cos(lat * deg2rad) * cos(declin * deg2rad) * cos(tsa * 
                                                                           deg2rad))/deg2rad
  azimuth <- acos((sin(lat * deg2rad) * sin((90 - zenith) * 
                                              deg2rad) - sin(declin * deg2rad))/cos(lat * deg2rad)/cos((90 - 
                                                                                                          zenith) * deg2rad))/deg2rad + 180
  azimuth <- ifelse(tsa > 0, azimuth%%360, 360 - azimuth%%360)
  
  return(zenith)
}
