#' SPCTRAL2 Solar Irradiance Model
#' 
#' An R-based implementation of the Bird and Riordan (1984, 1986) SPCTRAL2 model. Produces estimates direct and diffuse spectral irradiance on horizontal and titled surfaces at the Earth's surface for 300-4000 nm wavelengths. This R implementation is based on code for the FORTRAN program SPCTRLA2, including corrections to code and equations described by Bird and Riordan (1986). As of July 2021, FORTRAN and C implementations of SPCTRAL2 using the original Bird and Riordan (1984, 1986) formulation is available from the National Renewable Energy Laboratory: [https://www.nrel.gov/grid/solar-resource/spectral.html](https://www.nrel.gov/grid/solar-resource/spectral.html).
#' 
#' @param latitude Latitude in decimal degrees (East is positive, West is negative)
#' @param longitude Longitude in decimal degrees (North is positive, South is negative)
#' @param doy Day of year (numeric vector of length 1 with a values of 1-365)
#' @param hour Local hour
#' @param minute Local minute
#' @param alpha Power on angstrom turbidity expression (default = 1.14 for rural areas
#' @param angle_of_incidence Angle of incidence of direct beam on flat surface (in degrees)
#' @param ozone Ozone (O_3) in atmospheres per centimeter
#' @param albedo Ground albedo. Default is a generic value of 0.06 for sea surface
#' @param surface_tilt Tilt angle of ground surface from horizontal (degrees)
#' @param surface_azimuth Aziumuth of surface 
#' @param surface_pressure Surface pressure (milibars)
#' @param aod aerosol optical depth at 0.5 microns (base e)
#' @param water_vapor precipitable water vapor (cm)
#' @param asym Aerosol scattering asymmetry factor (forward to total scattering ratio)
#' @param omega Single scattering albedo at 400 nm. Default = 0.945
#' @param omega_p Single scattering wavelength variation factor. Default = 0.095
#' @param cloud_modification Logical. Should the cloud cover modification from Bird et al. (1987) be used. Default = TRUE.
#' 
#' @return A data frame containing spectral direct and diffuse irradiance by wavelength at the earth's surface.
#' 
#' @references Bird, R., Riordan, C. 1984. Simple solar spectral model for direct and diffuse irradiance on horizontal and titled planes at the Earth's surface for cloudless atmospheres. SERI/TR-2145-2436. Solar Energy Research Institute. 37 pp.
#' @references Bird, R.E., Riordan, C., 1986. Simple solar spectral model for direct and diffuse irradiance on horizontal and tilted planes at the Earth’s surface for cloudless atmospheres. J. Clim. Appl. Meteorol. 25, 87–97. https://doi.org/10.1175/1520-0450(1986)025<0087:SSSMFD>2.0.CO;2

spctral2 <- function(latitude, longitude, doy, hour, minute, alpha, angle_of_incidence, ozone = NA, albedo, surface_tilt, surface_azimuth = 180, surface_pressure, aod, water_vapor, asym = 0.65, omega = 0.945, omega_p = 0.095, cloud_modification = FALSE) {
  
  # Longitude flag
  if(longitude > 0) {
    lambda <- 1
  } else {
    lambda <- -1
  }
  
  # Photon flux constants
  CONST <- (1/(6.626176e-34*299792458))*10^-10
  
  surface_azimuth <- (surface_azimuth-180)/(180/3.14159)
  slope_tilt <- surface_tilt*pi/180
  
  # Calculate Solar Zenith Angle
  day_angle <- 6.283185*(doy-1)/365
  declination <- (0.006918-0.399912*cos(day_angle)+0.070257*sin(day_angle)-0.006758*cos(2*day_angle) +0.000907*sin(2*day_angle)-0.002697*cos(3*day_angle)+0.00148*sin(3*day_angle))*(180/3.14159)
  eqntime <- (0.000075+0.001868*cos(day_angle)-0.032077*sin(day_angle)-0.014615*cos(2*day_angle) -0.040849*sin(2*day_angle))*(229.18)
  hour_angle <- 15*(hour+minute/60+eqntime/60+((as.integer(floor(longitude/15))*15-longitude)*4)/60+12)-360
  zenith_degrees <- 180/pi*acos(cos(declination*pi/180)*cos(latitude*pi/180)*cos(hour_angle*pi/180)+sin(declination*pi/180)*sin(latitude*pi/180))
  
  # Calculate solar azimuth angle
  azimuth_degrees <- 180-(180/pi*asin(-1*cos(pi/180*declination)*sin(pi/180*hour_angle)/cos(pi/180*(90-zenith_degrees))))
  
  # Angle of Incidence (for a fixed surface)
  cos_aoi <- sin(declination*pi/180)*sin(latitude*pi/180)*cos(slope_tilt)-sin(declination*pi/180)*cos(latitude*pi/180)*sin(slope_tilt)*cos(surface_azimuth) + 
    cos(declination*pi/180)*cos(latitude*pi/180)*cos(slope_tilt)*cos(hour_angle*pi/180) + 
    cos(declination*pi/180)*sin(latitude*pi/180)*sin(slope_tilt)*cos(surface_azimuth)*cos(hour_angle*pi/180) + 
    cos(declination*pi/180)*sin(slope_tilt)*sin(surface_azimuth)*sin(hour_angle*pi/180)
  
  aoi <- 180/pi*acos(cos_aoi)
  
  # Earth-sun distance factor (Spencer 1971)
  rad_vec <- 1.00011+0.034221*cos(day_angle)+0.00128*sin(day_angle)+0.000719*cos(day_angle) +0.000077*sin(day_angle)
  
  # Effective airmass of ozone 
  amoz <- (1+22/6370)/(cos(zenith_degrees*pi/180)^2+2*(22/6370))^0.5
  
  # Geometric airmass path length
  geometric_airmass <- 1/(cos(pi/180*(zenith_degrees))+0.15*(93.885-zenith_degrees)^(-1.253)) 
  
  # Pressure-corrected airmass
  pressure_corr_airmass <- geometric_airmass*(surface_pressure/1013)
  
  # Total Ozone Thickness (Van Heuklon)
  if(latitude > 0) {
    oz_par <- c(0.15, 0.04, -30, 3, 1.28)
  } else {
    oz_par <- c(0.1, 0.03, 152.625, 2, 1.25)
  }
  
  if(lambda > 0 & latitude > 0) {
    oz_par <- c(oz_par, 20)
  } else if(lambda > 0 & latitude <= 0) {
    oz_par <- c(oz_par, -75)
  } else {
    oz_par <- c(oz_par, 0)
  }
  
  if(is.na(ozone)) {
    total_ozone <- 0.235+(oz_par[1] + oz_par[2]*sin(0.9865*(doy+oz_par[3])*0.017453)+0.02*sin(0.017453*(oz_par[4])*(oz_par[6])))*(sin(oz_par[5]*0.017453*latitude))^2    
  } else {
    total_ozone <- ozone
  }

  
  w_v <- c(0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.36, 0.37, 0.38, 
           0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 
           0.54, 0.55, 0.57, 0.593, 0.61, 0.63, 0.656, 0.668, 0.69, 0.71, 0.718, 0.724, 0.74, 0.753, 
           0.758, 0.763, 0.768, 0.78, 0.8, 0.816, 0.824, 0.832, 0.84, 0.86, 0.88, 0.905, 0.915, 
           0.925, 0.93, 0.937, 0.948, 0.965, 0.98, 0.994, 1.04, 1.07, 1.1, 1.12, 1.13, 1.145, 1.161, 
           1.17, 1.2, 1.24, 1.27, 1.29, 1.32, 1.35, 1.395, 1.443, 1.463, 1.477, 1.497, 1.52, 1.539, 
           1.558, 1.578, 1.592, 1.61, 1.63, 1.646, 1.678, 1.74, 1.8, 1.86, 1.92, 1.96, 1.985, 2.005, 
           2.035, 2.065, 2.1, 2.148, 2.198, 2.27, 2.36, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 
           3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4)
  e_v <- 1e-3*c(535.9, 558.3, 622, 692.7, 715.1, 832.9, 961.9, 931.9, 900.6, 911.3, 975.5, 975.9, 1119.9, 
           1103.8, 1033.8, 1479.1, 1701.3, 1740.4, 1587.2, 1837, 2005, 2043, 1987, 2027, 1896, 1909, 
           1927, 1831, 1891, 1898, 1892, 1840, 1768, 1728, 1658, 1524, 1531, 1420, 1399, 1374, 1373, 
           1298, 1269, 1245, 1223, 1205, 1183, 1148, 1091, 1062, 1038, 1022, 998.7, 947.2, 893.2, 
           868.2, 829.7, 830.3, 814, 786.9, 768.3, 767, 757.6, 688.1, 640.7, 606.2, 585.9, 570.2, 
           564.1, 544.2, 533.4, 501.6, 477.5, 442.7, 440, 416.8, 391.4, 358.9, 327.5, 317.5, 307.3, 
           300.4, 292.8, 275.5, 272.1, 259.3, 246.9, 244, 243.5, 234.8, 220.5, 190.8, 171.1, 144.5, 
           135.7, 123, 123.8, 113, 108.5, 97.5, 92.4, 82.4, 74.6, 68.3, 63.8, 49.5, 48.5, 38.6, 
           36.6, 32, 28.1, 24.8, 22.1, 19.6, 17.5, 15.7, 14.1, 12.7, 11.5, 10.4, 9.5, 8.6)
  absorption_water <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0.075, 0, 0, 0, 0, 0.016, 0.0125, 1.8, 2.5, 0.061, 
                        8e-04, 1e-04, 1e-05, 1e-05, 6e-04, 0.036, 1.6, 2.5, 0.5, 0.155, 1e-05, 
                        0.0026, 7, 5, 5, 27, 55, 45, 4, 1.48, 0.1, 1e-05, 0.001, 3.2, 115, 70, 75, 
                        10, 5, 2, 0.002, 0.002, 0.1, 4, 200, 1000, 185, 80, 80, 12, 0.16, 0.002, 
                        5e-04, 1e-04, 1e-05, 1e-04, 0.001, 0.01, 0.036, 1.1, 130, 1000, 500, 100, 4, 
                        2.9, 1, 0.4, 0.22, 0.25, 0.33, 0.5, 4, 80, 310, 15000, 22000, 8000, 650, 
                        240, 230, 100, 120, 19.5, 3.6, 3.1, 2.5, 1.4, 0.17, 0.0045)
  absorption_ozone <- c(10, 4.8, 2.7, 1.35, 0.8, 0.38, 0.16, 0.075, 0.04, 0.019, 0.007, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0.003, 0.006, 0.009, 0.014, 0.021, 0.03, 0.04, 0.048, 0.063, 
                        0.075, 0.085, 0.12, 0.119, 0.12, 0.09, 0.065, 0.051, 0.028, 0.018, 0.015, 
                        0.012, 0.01, 0.008, 0.007, 0.006, 0.005, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  absorption_umg <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0, 0, 0, 0, 0, 0, 4, 0.35, 0, 0, 0, 0, 
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.3, 
                      0.02, 2e-04, 0.00011, 1e-05, 0.05, 0.011, 0.005, 6e-04, 0, 0.005, 0.13, 0.04, 
                      0.06, 0.13, 0.001, 0.0014, 1e-04, 1e-05, 1e-05, 1e-04, 0.001, 4.3, 0.2, 21, 
                      0.13, 1, 0.08, 0.001, 0.00038, 0.001, 5e-04, 0.00015, 0.00014, 0.00066, 100, 
                      150, 0.13, 0.0095, 0.001, 0.8, 1.9, 1.3, 0.075, 0.01, 0.00195, 0.004, 0.29, 
                      0.025)
  
  # Transmission
  omegl <- omega*exp(-1*omega_p*(log(w_v/0.4))^2) # Single scattering albedo
  transmission_rayleigh <- exp(-1*pressure_corr_airmass/(w_v^4*(115.6406-1.335/w_v^2))) # Rayleigh
  transmission_toz <- exp(-1*absorption_ozone*amoz*total_ozone) # Total ozone
  transmission_umg <- exp(-1.41*absorption_umg*pressure_corr_airmass/(1+118.3*absorption_umg*pressure_corr_airmass)^0.45) # Uniform mixed gasses using 118.3 (Bird and Riordan 1986)
  transmission_water <- exp(-0.2385*absorption_water*water_vapor*geometric_airmass/(1+20.07*absorption_water*absorption_water*water_vapor*geometric_airmass)^0.45) # Water vapor
  aero1 <- aod*((w_v/0.5)^(-1*alpha))
  transmission_aerosol <- exp(-1*geometric_airmass*aero1) # Aerosol
  
  # Direct
  direct_wm2 <- e_v*rad_vec*transmission_rayleigh*transmission_toz*transmission_umg*transmission_water*transmission_aerosol
  
  # Diffuse scattering
  alg <- log(1-asym)
  afs <- alg*(1.459+alg*(0.1595+alg*0.4129))
  bfs <- alg*(0.0783+alg*(-0.3824-alg*0.5874))
  fsp <- 1-0.5*exp((afs+bfs/1.8)/1.8)
  
  total_as <- exp(-omegl*aero1*geometric_airmass)
  total_aa <- exp(-(1-omegl)*aero1*geometric_airmass)
  
  diff_ray <- e_v*cos(zenith_degrees*pi/180)*transmission_toz*transmission_umg*transmission_water*total_aa*(1-transmission_rayleigh^0.95)*0.5*rad_vec
  diff_aer <- e_v*cos(zenith_degrees*pi/180)*transmission_toz*transmission_umg*transmission_water*total_aa*transmission_rayleigh^1.5*(1-total_as)*(1-0.5*exp((afs+bfs*cos(zenith_degrees*pi/180))*cos(zenith_degrees*pi/180)))*rad_vec
  
  
  trp <- exp(-1.8/(w_v^4*(115.6406-1.335/w_v^2)))
  twp <- exp(-0.2385*absorption_water*water_vapor*1.8/((1+20.07*absorption_water*water_vapor*1.8)^0.45))
  tup <- exp(-1.41*absorption_umg*1.8/((1+118.93*absorption_umg*1.8)^0.45))
  tasp <- exp(-omegl*aero1*1.8)
  taap <- exp(-(1-omegl)*aero1*1.8)
  rhoa <- tup*twp*taap*(0.5*(1-trp)+(1-fsp)*trp*(1-tasp))
  drgd <- (direct_wm2*cos(zenith_degrees*pi/180)+(diff_ray+diff_aer))*albedo*rhoa/(1-albedo*rhoa)
  
  hz_diff <- numeric(length = length(w_v))
  hz_diff[w_v <= 0.45] <- (diff_ray[w_v <= 0.45]+diff_aer[w_v <= 0.45]+drgd[w_v <= 0.45])*(w_v[w_v <= 0.45]+0.55)^1.8
  hz_diff[w_v > 0.45] <- (diff_ray[w_v > 0.45]+diff_aer[w_v > 0.45]+drgd[w_v > 0.45])
  
  ref <- (direct_wm2*cos(zenith_degrees*pi/180)+hz_diff)*albedo*(1-cos(slope_tilt))/2
  difsc <- hz_diff*(direct_wm2/e_v)*cos_aoi/cos(zenith_degrees*pi/180)
  difsi <- 0.5*hz_diff*(1-(direct_wm2/e_v))*(1+cos(surface_tilt*pi/180))
  
  diffuse_wm2 <- ref+difsc+difsi
  
  # Cloud cover modification (Bird et al., 1987) *In development*
  # diffuse_wm2[w_v <= 0.55] <- diffuse_wm2[w_v <= 0.55] * (w_v[w_v <= 0.55]+0.45)^-1
  
  # Total energy
  total_wm2 <- direct_wm2*cos_aoi+diffuse_wm2
  
  # Photon flux density
  pfd_dni <- w_v * direct_wm2 * CONST
  pfd_dif <- w_v * diffuse_wm2 * CONST
  pfd_tot <- w_v * total_wm2 * CONST
  
  return(data.frame(wavelength = w_v,
              direct_pfd = pfd_dni,
              diffuse_pfd = pfd_dif,
              total_pfd = pfd_tot,
              direct_wm2_nm = direct_wm2,
              diffuse_wm2_nm = diffuse_wm2,
              total_wm2_nm = total_wm2))
}