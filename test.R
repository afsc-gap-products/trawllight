library(mk9process)
library(trawllight)

trawllight::tlu_setup_dir(survey = "BS", 
                          light_data_root = "G:/RACE_LIGHT/LightData/Data")


# Retrieve upcast and downcast data and run filter_stepwise, calculate_attenuation, assign_orientation
trawllight::tlu_get_casts(
  directory_structure = data.frame(path = read.csv(file = here::here("imports", "directories.csv"))[3:12,]),
  survey = "BS",
  cast.dir = "downcast",
  time.buffer = 20,
  silent = FALSE,
  bin.size = 2,
  bin.gap = 6,
  agg.fun = trawllight::geometric.mean)

trawllight::tlu_get_casts(
  directory_structure = data.frame(path = read.csv(file = here::here("imports", "directories.csv"))[3:12,]),
  survey = "BS",
  cast.dir = "upcast",
  time.buffer = 20,
  silent = FALSE,
  bin.size = 2,
  bin.gap = 6,
  agg.fun = trawllight::geometric.mean)

# Retrieve surface measurement data and calculate average near-surface irradiance for downcast start and upcast end.
trawllight::tlu_get_surface(
  data.frame(path = 
               read.csv(file = here::here("imports", "directories.csv"))[3:12,]),
                survey = "BS")

trawllight::tlu_run_trawllight()
