#' RADTRAN Solar Irradiance Model with delta-Eddington cloud transmission
#' 
#' An R-based implementation of the Gregg and Carder (1990) model for cloudless maritime atmospheres, along with an option to use the delta-Eddington approximation of two-stream approach for spectral cloud transmission following Slingo (1989). Produces estimates of direct and diffuse spectral irradiance at the sea surface and just below the sea surface at 1 nm resolution for 350-700 nm wavelengths. Atmospheric path length calculation follows Kasten and Young (1989) instead of Kasten (1965), as suggested by Gregg (2002). Ozone scale height calculation 
#' 
#' @param latitude Latitude in decimal degrees (East is positive, West is negative)
#' @param longitude Longitude in decimal degrees (North is positive, South is negative)
#' @param doy Day of year (numeric vector of length 1 with a values of 1-365)
#' @param hour Local hour
#' @param minute Local minute
#' @param alpha Power on angstrom turbidity expression 
#' @param surface_pressure Surface pressure (milibars)
#' @param water_vapor precipitable water vapor (cm)
#' @param aerosol_scale_height Aerosol scale height in kilometers. Only used if cloud_modification = FALSE.
#' @param visual_range Visibility in kilometers. Only used if cloud_modification = FALSE.
#' @param air_mass Air mass type (1-10). Typically 1-10, where 1 correspondswith marine aerosols, 10 with continental aerosols. 
#' @param relative_humidity Relative humidity (percent).
#' @param wind_speed Wind speed in meters per second. Must be greater than or equal to zero. Default = 0 (calm)
#' @param cloud_modification Logical. Use the cloud cover modification from Gregg (2002), based on Slingo (1989)? Default = TRUE.
#' @param below_surface Logical. TRUE = estimate irradiance just below the sea surface. FALSE = estimate irradiance just above the sea surface.
#' @param water_refraction Refractive index of water. Default for seawater = 1.341.
#' @param foam_spectrum Logical. TRUE = use Gregg (2002) wavelength-dependent correction for sea-foam reflectance-based on Frouin et al. (1996). FALSE = original implementation with no spectral dependence (Gregg and Carder 1990)
#' @param cloud_droplet_radius Cloud droplet radius in micrometers. Used for cloud_modification. Default 11.8 micrometers = mean oceanic value from Han et al. (1994).
#' @param liquid_water_path Liquid water path used for cloud_modification.
#' @param ozone Either a 1L numeric vector indicating the total ozone value or a character vector indicating an algorithmic approximation method. The only approximation method currently available is "vanheuklon" based on Van Heuklon (1979).
#' 
#' @return A data frame containing spectral direct and diffuse irradiance by wavelength at the sea surface and just below the sea surface.
#' 
#' @references Gregg, W.W., 2002. A coupled ocean-atmosphere radiative model for global ocean biogeochemical models. NASA Global Modeling and Assimilation Series, M. Suarez (ed.), NASA Technical Memorandum 2002-104606, Volume 22.
#' @references Gregg, W.W., Carder, K.L., 1990. A simple spectral solar irradiance model for cloudless maritime atmospheres. Limnol. Oceanogr. 35, 1657–1675. https://doi.org/10.4319/lo.1990.35.8.1657
#' @references Kasten, F., 1965. A new table and approximation formula for the relative optial air mass. Arch. für Meteorol. Geophys. und Bioklimatologie Ser. B 14, 206–223. https://doi.org/10.1007/BF02248840
#' @references Kasten, F., Young, A.T., 1989. Revised optical air mass tables and approximation formula. Appl. Opt. 28, 4735–4738. https://doi.org/10.1364/AO.28.004735
#' @references Slingo, A., 1989. A GCM parameterization for the shortwave radiative properties of water clouds. J. Atmos. Sci. 46, 1419–1427. https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2
#' @references Van Heuklon, T.K., 1979. Estimating atmospheric ozone for solar radiation models. Sol. Energy 22, 63–68. https://doi.org/10.1016/0038-092X(79)90060-4

solar_radtran <- function(latitude, longitude, doy, hour, minute, alpha, surface_pressure, water_vapor, aerosol_scale_height = NA, visual_range = NA, air_mass, relative_humidity, wind_speed = 0, cloud_modification = TRUE, cloud_droplet_radius = 11.8, liquid_water_path = 125, water_refraction = 1.341, foam_spectrum = FALSE, ozone = "vanheuklon") {
  
  # Numbers in parentheses correspond with equations in Gregg and Carder (1990) 
  
  # Photon flux constant
  CONST <- (1/(6.626176e-34*299792458))*10^-10
  
  # Wavelength (um)
  w_v <- c(350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 
           368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 
           386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 
           404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 
           422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 
           440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 
           458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 
           476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 
           494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 
           512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 
           530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 
           548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 
           566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 
           584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 
           602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 
           620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 
           638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 
           656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 
           674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 
           692, 693, 694, 695, 696, 697, 698, 699, 700)/1000
  
  # Extraterrestrial solar radiation W m^-2 um^-1)
  e_v <- c(0.961, 0.953, 0.949, 1.056, 1.122, 1.078, 1.047, 0.879, 0.752, 0.919, 1.062, 1.054, 1.047, 
           1.024, 0.998, 1.108, 1.259, 1.221, 1.156, 1.184, 1.197, 1.162, 1.144, 1.027, 0.953, 1.004, 
           1.004, 1.317, 1.317, 1.141, 1.139, 1.115, 1.083, 0.821, 0.858, 1.029, 1.026, 0.995, 1.01, 
           1.145, 1.152, 1.263, 1.115, 0.733, 0.852, 1.25, 1.071, 0.853, 1.25, 1.515, 1.674, 1.721, 
           1.799, 1.719, 1.638, 1.651, 1.663, 1.681, 1.698, 1.65, 1.621, 1.74, 1.812, 1.755, 1.74, 
           1.781, 1.791, 1.715, 1.701, 1.663, 1.124, 1.823, 1.76, 1.657, 1.693, 1.748, 1.691, 1.673, 
           1.656, 1.65, 1.407, 1.351, 1.727, 1.805, 1.69, 1.767, 1.835, 1.845, 1.792, 1.673, 1.711, 
           1.796, 1.892, 1.957, 1.961, 1.963, 1.856, 1.874, 2.036, 2.054, 2.135, 2.111, 2.004, 2.007, 
           2.024, 2.03, 2.066, 2.06, 2.028, 2.028, 2.029, 2.039, 2.101, 2.086, 1.992, 1.987, 1.959, 
           1.966, 2.01, 2.001, 1.946, 1.957, 2.022, 2.025, 2.038, 2.029, 1.982, 1.996, 2.063, 2.064, 
           2.067, 2.065, 2.054, 2.047, 2.011, 1.95, 1.687, 1.723, 1.874, 1.949, 1.938, 1.9, 1.909, 
           1.941, 1.954, 2.003, 2.003, 2.003, 1.973, 1.933, 1.871, 1.832, 1.89, 1.928, 1.925, 1.924, 
           1.956, 1.977, 1.953, 1.941, 1.939, 1.939, 1.913, 1.9, 1.872, 1.858, 1.744, 1.688, 1.725, 
           1.743, 1.828, 1.862, 1.891, 1.908, 1.922, 1.873, 1.843, 1.834, 1.83, 1.921, 1.959, 1.952, 
           1.948, 1.872, 1.859, 1.951, 1.941, 1.861, 1.858, 1.851, 1.84, 1.819, 1.847, 1.875, 1.893, 
           1.911, 1.889, 1.867, 1.875, 1.883, 1.878, 1.874, 1.87, 1.869, 1.889, 1.896, 1.84, 1.826, 
           1.817, 1.815, 1.823, 1.824, 1.86, 1.868, 1.848, 1.844, 1.844, 1.844, 1.852, 1.853, 1.805, 
           1.799, 1.884, 1.894, 1.861, 1.857, 1.857, 1.857, 1.822, 1.823, 1.852, 1.853, 1.866, 1.866, 
           1.861, 1.859, 1.81, 1.808, 1.765, 1.765, 1.761, 1.766, 1.797, 1.797, 1.794, 1.798, 1.818, 
           1.81, 1.763, 1.761, 1.752, 1.747, 1.729, 1.737, 1.772, 1.768, 1.751, 1.749, 1.742, 1.734, 
           1.726, 1.735, 1.744, 1.726, 1.709, 1.693, 1.677, 1.705, 1.733, 1.732, 1.731, 1.717, 1.704, 
           1.684, 1.666, 1.669, 1.671, 1.689, 1.701, 1.674, 1.656, 1.655, 1.654, 1.654, 1.654, 1.658, 
           1.661, 1.662, 1.663, 1.643, 1.63, 1.622, 1.616, 1.624, 1.629, 1.619, 1.612, 1.611, 1.61, 
           1.582, 1.564, 1.585, 1.6, 1.599, 1.598, 1.462, 1.371, 1.377, 1.415, 1.46, 1.505, 1.548, 
           1.581, 1.584, 1.576, 1.566, 1.557, 1.55, 1.543, 1.537, 1.531, 1.525, 1.519, 1.512, 1.506, 
           1.5, 1.494, 1.488, 1.481, 1.476, 1.472, 1.469, 1.466, 1.463, 1.46, 1.457, 1.454, 1.45, 
           1.447, 1.444, 1.441, 1.438, 1.435, 1.432, 1.429, 1.426, 1.423, 1.42, 1.417, 1.414, 1.411)
  
  # Absorption: Water
  absorption_water <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.006, 0.014, 0.022, 0.03, 
                        0.034, 0.032, 0.027, 0.022, 0.016, 0.012, 0.01, 0.011, 0.01, 0.008, 0, 0, 0, 
                        0, 0, 0.003, 0.123, 0.284, 0.454, 0.605, 0.7, 0.697, 0.636, 0.549, 0.454, 
                        0.361, 0.278, 0.202, 0.132, 0.069, 0.023, 0.022, 0.017, 0.011, 0.001, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001, 0.002, 
                        0.003, 0.004, 0.005, 0.004, 0.003, 0.002, 0.001, 0, 0, 0.001, 0.001, 0.002, 
                        0.001, 0, 0, 0, 0, 0.011, 0.038, 0.074, 0.112, 0.147, 0.173, 0.181, 0.179, 
                        0.171, 0.16, 0.142, 0.117, 0.087, 0.057, 0.029, 0.008, 0, 0, 0, 0, 0.001, 
                        0.002, 0.002, 0.001, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001, 0.001, 
                        0.001, 0.001, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.084, 0.196, 0.317, 
                        0.434, 0.53, 0.588, 0.621, 0.637, 0.636, 0.602)
  
  # Absorption: Ozone
  absorption_ozone <- c(0.009, 0.012, 0.009, 0.006, 0.004, 0.002, 0.002, 0.001, 0.001, 0.001, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001, 0.001, 0.001, 0.001, 
                        0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.003, 
                        0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.004, 0.004, 
                        0.005, 0.005, 0.006, 0.007, 0.007, 0.008, 0.008, 0.009, 0.009, 0.009, 0.009, 
                        0.008, 0.008, 0.007, 0.001, 0.007, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 
                        0.013, 0.014, 0.015, 0.017, 0.018, 0.019, 0.02, 0.02, 0.021, 0.02, 0.02, 
                        0.019, 0.019, 0.018, 0.018, 0.019, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 
                        0.025, 0.026, 0.028, 0.03, 0.032, 0.034, 0.036, 0.038, 0.039, 0.041, 0.041, 
                        0.041, 0.04, 0.039, 0.038, 0.037, 0.037, 0.038, 0.039, 0.041, 0.042, 0.044, 
                        0.045, 0.047, 0.049, 0.051, 0.052, 0.054, 0.056, 0.058, 0.059, 0.061, 0.063, 
                        0.064, 0.066, 0.068, 0.069, 0.07, 0.071, 0.072, 0.072, 0.072, 0.072, 0.073, 
                        0.074, 0.075, 0.076, 0.077, 0.079, 0.08, 0.081, 0.083, 0.084, 0.085, 0.086, 
                        0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.094, 0.096, 0.098, 0.1, 0.103, 
                        0.105, 0.107, 0.109, 0.112, 0.113, 0.115, 0.116, 0.117, 0.118, 0.119, 0.119, 
                        0.119, 0.119, 0.119, 0.119, 0.118, 0.117, 0.116, 0.115, 0.114, 0.112, 0.111, 
                        0.11, 0.109, 0.109, 0.108, 0.108, 0.109, 0.11, 0.112, 0.113, 0.115, 0.117, 
                        0.119, 0.121, 0.122, 0.124, 0.125, 0.125, 0.125, 0.125, 0.124, 0.123, 0.122, 
                        0.12, 0.119, 0.118, 0.116, 0.115, 0.114, 0.112, 0.111, 0.11, 0.108, 0.107, 
                        0.105, 0.104, 0.103, 0.101, 0.1, 0.099, 0.097, 0.096, 0.094, 0.093, 0.092, 
                        0.09, 0.089, 0.088, 0.086, 0.085, 0.083, 0.082, 0.081, 0.079, 0.078, 0.077, 
                        0.075, 0.074, 0.073, 0.071, 0.07, 0.068, 0.067, 0.066, 0.065, 0.064, 0.063, 
                        0.062, 0.061, 0.06, 0.059, 0.058, 0.057, 0.056, 0.056, 0.055, 0.054, 0.053, 
                        0.052, 0.051, 0.051, 0.05, 0.049, 0.048, 0.047, 0.046, 0.045, 0.045, 0.044, 
                        0.043, 0.042, 0.041, 0.04, 0.04, 0.039, 0.038, 0.037, 0.036, 0.035, 0.034, 
                        0.034, 0.033, 0.032, 0.031, 0.03, 0.03, 0.029, 0.028, 0.027, 0.026, 0.025, 
                        0.024, 0.024, 0.023, 0.022, 0.022)
  
  # Absorption: Oxygen/Uniform mixed gasses
  absorption_umg <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0.002, 0.005, 0.008, 0.01, 0.011, 0.01, 0.008, 0.005, 0.002, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.067, 
                        0.81, 0.65, 0.505, 0.36, 0.325, 0.248, 0.157, 0.068, 0.001, 0, 0, 0, 0, 0)
  
  # Calculate Solar Zenith Angle
  day_angle <- 6.283185*(doy-1)/365
  declination <- (0.006918-0.399912*cos(day_angle)+0.070257*sin(day_angle)-0.006758*cos(2*day_angle) +0.000907*sin(2*day_angle)-0.002697*cos(3*day_angle)+0.00148*sin(3*day_angle))*(180/3.14159)
  eqntime <- (0.000075+0.001868*cos(day_angle)-0.032077*sin(day_angle)-0.014615*cos(2*day_angle) -0.040849*sin(2*day_angle))*(229.18)
  hour_angle <- 15*(hour+minute/60+eqntime/60+((as.integer(floor(longitude/15))*15-longitude)*4)/60+12)-360
  zenith_radians <- acos(cos(declination*pi/180)*cos(latitude*pi/180)*cos(hour_angle*pi/180)+sin(declination*pi/180)*sin(latitude*pi/180))
  zenith_degrees <- 180/pi*zenith_radians
  
  cos_z <- cos(zenith_radians)
  
  # (12) Extraterrestrial irradiance correction for earth-sun distance (Gordon 1983)
  etc_factor <- (1 + 0.0167 * cos((2*pi*(doy-3)/365)))^2
  
  # (13) Geometric airmass path length (Kasten and Young 1989)
  geometric_airmass <- 1/(cos_z+0.50572*(96.07995-zenith_degrees)^(-1.16364))
  
  # (14) Effective airmass of ozone (Paltridge and Platt 1976)
  amoz <- (1+22/6370)/(cos_z^2+2*(22/6370))^0.5
  
  # (16) Pressure-corrected airmass
  pressure_corr_airmass <- geometric_airmass*(surface_pressure/1013)
  
  # Total ozone scale height (Van Heuklon 1979)
  if(ozone == "vanheuklon") {
  if(latitude > 0) {
    oz_par <- c(0.15, 0.04, -30, 3, 1.28)
  } else {
    oz_par <- c(0.1, 0.03, 152.625, 2, 1.25)
  }
  
  if(longitude > 0 & latitude > 0) {
    oz_par <- c(oz_par, 20)
  } else if(longitude > 0 & latitude <= 0) {
    oz_par <- c(oz_par, -75)
  } else {
    oz_par <- c(oz_par, 0)
  }
  
  total_ozone <- 0.235+(oz_par[1] + oz_par[2]*sin(0.9865*(doy+oz_par[3])*0.017453)+0.02*sin(0.017453*(oz_par[4])*(oz_par[6])))*(sin(oz_par[5]*0.017453*latitude))^2
  
  } else {
    total_ozone <- ozone
  }

  
  # (17) Transmittance: ozone (absorption coefficients from Inn and Tanaka [1953])
  transmission_toz <- exp(-1*absorption_ozone*amoz*total_ozone)
  
  # (18) Transmittance: uniform mixed gasses (Bird and Riordan 1986)
  transmission_umg <- exp(-1.41*absorption_umg*pressure_corr_airmass/(1+118.3*absorption_umg*pressure_corr_airmass)^0.45) 
  
  # (19) Transmittance: Water (Bird and Riordan 1986)
  transmission_water <- exp(-0.2385*absorption_water*water_vapor*geometric_airmass/(1+20.07*absorption_water*absorption_water*water_vapor*geometric_airmass)^0.45)
  
  # (35) Asymmetry parameter for scattering phase function
  if(alpha <0) {
    asym <- 0.82
  } else if(alpha > 1.2) {
    asym <- 1.2
  } else{
    asym <- -0.1417*alpha+0.82
  }

  # Forward scattering probability
  b_3 <- log(1-asym)
  b_1 <- b_3*(1.459+b_3*(0.1595+b_3*0.4129))
  b_2 <- b_3*(0.0783+b_3*(-0.3824-b_3*0.5874))
  fsp <- 1-0.5*exp((b_1+b_2*cos_z)*cos_z)
  
  # (36) Single-scattering albedo
  omega_a <- (-0.0032*air_mass+0.972) * exp(3.06*relative_humidity*10^-4)
  
  # (39-43) Sea-foam reflectance as a function of wind stress (Trenberth et al. 1989)
  if(wind_speed <= 4) {
    rho_f <- 0
  } else if(wind_speed >4 & wind_speed <= 7) {
    drag <- (0.62+1.56*wind_speed^-1)*1e-3
    rho_f <- (2.2e-5*1.2e3*drag*wind_speed^2)-4e-4
  } else {
    drag <- (0.49+0.065*windspeed)*1e-3
    rho_f <- (4.5e-5*1.2e3*drag-4e-5)*wind_speed^2
  }
  
  # Foam spectral reflectance (Frouin et al. 1996; Gregg 2002) *** In development ***
  # if(foam_spectrum) {
  #   a_seawater
  #   b_seawater
  #   transmittance_seawater <- exp(a_seawater + 0.5*b_seawater) 
  #   foam_spectral_correction <- rep(1, length(w_v))
  #   foam_spectral_correction[w_v < 0.9] <- 0.9976 + 0.2194 * log(transmittance_seawater[w_v < 0.9]) + 0.0554 * log(transmittance_seawater[w_v < 0.9])^2 +
  #     0.0067 * log(transmittance_seawater[w_v < 0.9])^3
  #   foam_spectral_correction[w_v >= 0.9] <- 5.026 - 0.0114*(w_v[w_v >= 0.9]*1000) + 9.552e-6*(w_v[w_v >= 0.9]*1000)^2 - 2.698e-9*(w_v[w_v >= 0.9]*1000)^3  
  # } else {
  #   foam_spectral_correction <- 1
  # }
  
  foam_spectral_correction <- 1
  
  rho_f <- rho_f * foam_spectral_correction

  # (44-47) Direct specular reflectance
  rho_dsp <- 0
  
  if(zenith_degrees >= 40) {
    if(wind_speed > 2) {
      rho_dsp <- 0.0253*exp((-7.14e-4*wind_speed+0.0618)*(zenith_degrees-40))
    } else {
      refracted_angle <- asin(sin(zenith_degrees)/water_refraction)
      rho_dsp <- 0.5*(0.5*sin(zenith_radians-refracted_angle)^2/sin(zenith_radians+refracted_angle)^2+
                        tan(zenith_radians-refracted_angle)^2/tan(zenith_radians+refracted_angle)^2)
    }
  }
  
  # Diffuse specular reflectance (Burt 1954)
  if(wind_speed > 4) {
    rho_ssp <- 0.057
  } else {
    rho_ssp <- 0.066
  }
  
  # (37) Surface reflectance (direct)
  rho_direct <- rho_dsp + rho_f
  
  # (38) Surface reflectance (diffuse)
  rho_diffuse <- rho_ssp + rho_f
  
  # Cloudless sky = 100% transmisson
  transmission_cloud_direct <- 1
  transmission_cloud_diffuse <- 1
  
  # Cloudy sky (Slingo 1989, Gregg 2002)
  if(cloud_modification) {
    
    # Parameters from Slingo (1989) Table 1.
    cloud_a_i <- 1e-2*c(3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 3.308, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.801, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.668, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.698, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.672, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.838, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.831, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 2.895, 3.115)
    cloud_b_i <- c(1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.293, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.307, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.317, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.244)
    cloud_c_i <- c(-3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, -3e-07, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-06, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -1.2e-07, -2.7e-07)
    cloud_d_i <- c(2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 2.36e-07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 4.4e-07, 1.4e-06)
    cloud_e_i <- c(0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.839, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.836, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.825, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.828, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.818, 0.804)
    cloud_f_i <- 1e-3*c(1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 1.946, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 2.153, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 1.881, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 3.004, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.467, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.776, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.492, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 2.989, 3.52)
    
    # Cloud Optical Depth (Slingo 1989 - Eqn. 1) with a mean cloud droplet radius
    tau_clouds <- liquid_water_path * (cloud_a_i + cloud_b_i/cloud_droplet_radius)
    
    # Single Scatter Albedo (Slingo 1989 - Eqn. 2)
    omega_clouds <- 1 - cloud_c_i - cloud_d_i*cloud_droplet_radius
    
    # Asymmetry parameter (g) (Slingo 1989 - Eqn. 3)
    asym_clouds <- cloud_e_i + cloud_f_i*cloud_droplet_radius
    
    # Fraction of scattered diffuse radiation scattered backwards (Slingo 1989 - Eqn. 6)
    beta_diffuse <- 3/7*(1-asym_clouds)
    
    # Fraction of direct radiation scattered backwards (Slingo 1989 - Eqn. 7)
    beta_direct <- 0.5 - 3*cos_z*asym_clouds/(4*(1+asym_clouds))
    
    # Fraction of scattered direct flux emerging at zenith angles close to incident beam (Slingo 1989 - Eqn. 8)
    f_foreward <- asym_clouds^2
    
    # Reciprocal of cosine of diffuse upward (Slingo 1989 - Eqn. 9)
    u_1 <- 7/4
    
    # Reciprocal of cosine of diffuse downward (Slingo 1989 - Eqn. 10)
    u_2 <- 7/4*(1-(1-omega_clouds)/(7*omega_clouds*beta_diffuse))
  
    # (Slingo 1989 - Eqn. 11)
    a_1 <- u_1*(1-omega_clouds*(1-beta_diffuse))
    
    # (Slingo 1989 - Eqn. 12)
    a_2 <- u_2*omega_clouds*beta_diffuse
    
    # (Slingo 1989 - Eqn. 13)
    a_3 <- (1-f_foreward) * omega_clouds * beta_direct
    
    # (Slingo 1989 - Eqn. 14)
    a_4 <- (1-f_foreward) * omega_clouds * (1-beta_direct)
    
    # (Slingo 1989 - Eqn. 15)
    epsilon <- sqrt(abs(a_1^2 - a_2^2))
    
    # (Slingo 1989 - Eqn. 16)
    M <- a_2 / (a_1+epsilon)
    
    # (Slingo 1989 - Eqn. 17)
    E <- exp(-1*epsilon*tau_clouds)
    
    # (Slingo 1989 - Eqn. 18)
    gamma_1 <- ((1-omega_clouds*f_foreward)*a_3-cos_z*(a_1*a_3+a_2*a_4))/((1-omega_clouds*f_foreward)^2-epsilon^2*cos_z^2)
    
    # (Slingo 1989 - Eqn. 19)
    gamma_2 <- -1*((1-omega_clouds*f_foreward)*a_4-cos_z*(a_1*a_4+a_2*a_3))/((1-omega_clouds*f_foreward)^2-epsilon^2*cos_z^2)
    
    # Direct cloud transmission (Slingo 1989 - Eqn. 20)
    transmission_cloud_direct <- exp(-1*(1-omega_clouds*f_foreward)*tau_clouds/cos_z)
    
    # Diffuse reflectivity for diffuse incident radiation (Slingo 1989 - Eqn. 21)
    ref_diffuse <- M*(1-E^2)/(1-E^2*M^2)
    
    # Diffuse transmissiviity for diffuse incident radiation (Slingo 1989 - Eqn. 22)
    transmission_cloud_diffuse <- E*(1-M^2)/(1-E^2*M^2)
    
    # Diffuse reflectivity for direct incident radiation (Slingo 1989 - Eqn. 23)
    ref_direct <- -1*gamma_2*ref_diffuse-gamma_1*transmission_cloud_diffuse*transmission_cloud_direct + gamma_1
    
    # Diffuse transmissivity for direct incident radiation
    transmission_cloud_direct_incident <- -1*gamma_2*transmission_cloud_diffuse - gamma_1*transmission_cloud_diffuse*ref_diffuse + gamma_2*transmission_cloud_direct 
    
    # Total transmission to direct
    transmission_cloud_total_direct <- transmission_cloud_direct + transmission_cloud_direct_incident
    
    # Aerosol and Rayleign transmission ignored in the presence of clouds
    
    transmission_rayleigh <- 1
    transmission_aerosol <- 1
    transmission_as <- 0
    transmission_aa <- 1
    
  } else {
    
    # (28) Kochsmieder formula for aerosol beam attenuation at 550 (Fitzgerald 1989)
    ca_550 <- 3.91/visual_range
    
    # (29) Aerosol optical thickness
    taua_550 <- ca_550*aerosol_scale_height
    
    # (Rearranging 27) - Calculate turbidity coefficient (beta)
    beta <- taua_550 * 0.55^alpha
    
    # (27) Wavelength-dependent aerosol optical thickness
    tau_a <- beta*(w_v^-1) 
    
    # (15) Transmittance: Rayleigh (Bird and Riordan 1986)
    transmission_rayleigh <- exp(-1*pressure_corr_airmass/(w_v^4*(115.6406*-1.335/w_v^2)))
    
    # (30) Transmittance: Aerosol 
    transmission_aerosol <- exp(-1*geometric_airmass*tau_a)
    
    # (7) Transmittance due to aerosol scattering
    transmission_as <- exp(-1*omega_a*tau_a*geometric_airmass)
    
    # (5) Transmittance after aerosol absorption
    transmission_aa <- exp(-1*(1-omega_a)*tau_a*geometric_airmass)
  }
  
  # (9) Direct downwelling irradiance
  direct_surface_wm2 <- e_v*etc_factor*cos_z*transmission_aerosol*transmission_umg*transmission_water*transmission_toz*transmission_rayleigh*transmission_cloud_total_direct
  
  direct_subsurface_wm2 <- direct_surface_wm2*(1-rho_direct)
  
  # (4) Diffuse component from Rayleigh scattering
  diffuse_rayleigh <- e_v*etc_factor*cos_z*transmission_toz*transmission_umg*transmission_water*transmission_aa*(1-transmission_rayleigh^0.95)*0.5
  
  diffuse_scattering <- e_v*etc_factor*cos_z*transmission_toz*transmission_umg*transmission_water*transmission_rayleigh*(1-transmission_as)*fsp*transmission_cloud_diffuse
  
  # (10) Diffuse downwelling irradiance
  diffuse_surface_wm2 <- diffuse_rayleigh + diffuse_scattering
  diffuse_subsurface_wm2 <- diffuse_surface_wm2 * (1-rho_diffuse)
  
  diffuse_surface_pfd <- w_v * diffuse_surface_wm2 * CONST
  diffuse_subsurface_pfd <- w_v * diffuse_subsurface_wm2 * CONST
  
  direct_surface_pfd <- w_v * direct_surface_wm2 * CONST
  direct_subsurface_pfd <- w_v * direct_subsurface_wm2 * CONST
  
  return(data.frame(wavelength = w_v,
                    direct_surface_pfd = direct_surface_pfd, 
                    direct_subsurface_pfd = direct_subsurface_pfd,
                    diffuse_surface_pfd = diffuse_surface_pfd,
                    diffuse_subsurface_pfd = diffuse_subsurface_pfd,
                    direct_surface_wm2_nm = direct_surface_wm2,
                    direct_subsurface_wm2_nm = direct_subsurface_wm2,
                    diffuse_surface_wm2_nm = diffuse_surface_wm2,
                    diffuse_subsurface_wm2_nm = diffuse_subsurface_wm2))
}