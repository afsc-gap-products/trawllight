#' Function to generate table of sunrise times and maps
#'
#' @param min_date Minimum date as a character vector in YYYY-MM-DD format ("2023-05-31").
#' @param max_date Maximum date as a character vector  in YYYY-MM-DD format.
#' @param dates Optional. Must provide dates or a min/max date. Dates as a POSIXct vector.
#' @param select_stations Optional character vector of stations.
#' @param region Region. Passed to akgfmaps::get_base_layers()
#' @param timezone Valid time zone for POSIXct objects as a character vector (default = "US/Alaska)
#' @export

make_sunrise_table <- function(min_date,
                               max_date,
                               dates = NULL,
                               select_stations = NULL,
                               region,
                               timezone = "US/Alaska") {
  
  extract_dt <- function(x, index) {

    out <-
      as.numeric(
        paste(
          unlist(
            strsplit(x = as.character(x), split = ""))[index],
          collapse = "")
      )

    return(out)

  }

  format_date <- function(x) {
    
    if(x > 24) {
      x <- x - 24
    }

    hour <- unlist(strsplit(as.character(floor(x)), split = ""))

    hour_vec <- c("0", "0")

    if(length(hour) == 2) {
      hour_vec <- hour
    } else {
      hour_vec[2] <- hour
    }

    min_vec <- c("0", "0")

    minutes <- unlist(strsplit(as.character(floor(x%%1*60)), split = ""))

    if(length(minutes) == 2) {
      min_vec <- minutes
    } else {
      min_vec[2] <- minutes
    }

    out <- paste0(paste(hour_vec, collapse = ""), ":", paste(min_vec, collapse = ""))

    return(out)

  }

  map_layers <- akgfmaps::get_base_layers(select.region = region, set.crs = "EPSG:4326")

  map_layers$survey.grid <- sf::st_make_valid(map_layers$survey.grid)

  stn_grid <- dplyr::bind_cols(map_layers$survey.grid,
                               as.data.frame(sf::st_coordinates(sf::st_centroid(map_layers$survey.grid))))

  if(!is.null(select_stations)) {
    stn_grid <- dplyr::filter(stn_grid, STATIONID %in% select_stations)
  }

  n_stn <- nrow(stn_grid)


  if(timezone == "US/Alaska") {
    sel_tz <- -8
  }

  if(is.null(dates)) {
    min_date <- as.POSIXct(x = as.character(min_date), tz = timezone)
    max_date <- as.POSIXct(x = as.character(max_date), tz = timezone)

    dates <- seq.POSIXt(min_date, max_date, by = "day")
  } else {
    if(!all(class(dates) == "POSIXct")) {
      dates <- as.POSIXct(x = as.character(dates), tz = timezone)
    } else {
      dates <- lubridate::force_tz(dates, tzone = timezone)
    }
  }


  time_df <- data.frame()

  for(ii in 1:length(dates)) {

    day <- extract_dt(x = dates[ii], index = 9:10)
    month <- extract_dt(x = dates[ii], index = 6:7)
    year <- extract_dt(x = dates[ii], index = 1:4)

    date_vars <- fishmethods::astrocalc4r(day = rep(day, n_stn),
                                          month = rep(month, n_stn),
                                          year = rep(year, n_stn),
                                          hour = rep(1, n_stn),
                                          timezone = rep(sel_tz, n_stn),
                                          lat = stn_grid$Y,
                                          lon = stn_grid$X)[, c('sunrise', 'sunset', 'daylight')]

    date_vars$STATIONID <- stn_grid$STATIONID

    date_vars$sunrise <- as.POSIXct(x =
                                      paste0(
                                        as.Date(dates[ii]), " ",
                                        unlist(
                                          lapply(X = date_vars$sunrise,
                                                 FUN = format_date)
                                        )
                                      ),
                                    tz = timezone
    )

    date_vars$sunset <- date_vars$sunrise + 3600 * date_vars$daylight


    time_df <- rbind(time_df, date_vars)

  }


  time_df <- time_df |>
    dplyr::mutate(sunrise_hr = as.numeric(lubridate::hour(sunrise) + lubridate::minute(sunrise)/60),
                  sunset_hr = as.numeric(lubridate::hour(sunset) + lubridate::minute(sunset)/60))
  
  time_df$sunset_hr[time_df$sunset_hr < 6] <- 23.99
  
  time_df <- time_df |>
    dplyr::group_by(STATIONID) |>
    dplyr::summarise(max_sunrise = max(sunrise_hr),
                     min_sunset = min(sunset_hr))


  sunrise_range <- dplyr::inner_join(stn_grid,
                                     time_df)

  sunrise_range_centroid <- sf::st_centroid(sunrise_range)
  sunrise_range_centroid$label_sunrise <- unlist(
    lapply(X = sunrise_range_centroid$max_sunrise,
           FUN = format_date)
  )
  
  sunrise_range_centroid$label_sunset <- unlist(
    lapply(X = sunrise_range_centroid$min_sunset,
           FUN = format_date)
  )

  png(here::here(paste0("latest_sunrise_", min(dates), "_", max(dates), ".png")), width = 10, height = 10, units = "in", res = 300)
  print(
    ggplot() +
      geom_sf(data = map_layers$akland |>
                dplyr::filter(is.na(DESC_))) +
      geom_sf(data = map_layers$survey.grid,
              color = "black",
              fill = NA) +
      geom_sf_text(data = sunrise_range_centroid,
                   mapping = aes(label = label_sunrise,
                                 color = max_sunrise),
                   size = rel(2.2)) +
      scale_color_gradient(name = "Sunrise (AKT)") +
      scale_x_continuous(limits = c(map_layers$plot.boundary$x)) +
      scale_y_continuous(limits = c(map_layers$plot.boundary$y)) +
      ggtitle(label = paste0("Latest sunrise (", min(dates), " to ", max(dates), ")")) +
      theme_bw() +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
  )
  dev.off()
  
  
  png(here::here(paste0("earliest_sunset_", min(dates), "_", max(dates), ".png")), width = 10, height = 10, units = "in", res = 300)
  print(
    ggplot() +
      geom_sf(data = map_layers$akland |>
                dplyr::filter(is.na(DESC_))) +
      geom_sf(data = map_layers$survey.grid,
              color = "black",
              fill = NA) +
      geom_sf_text(data = sunrise_range_centroid,
                   mapping = aes(label = label_sunset,
                                 color = min_sunset),
                   size = rel(2.2)) +
      scale_color_gradient(name = "Sunset (AKT)") +
      scale_x_continuous(limits = c(map_layers$plot.boundary$x)) +
      scale_y_continuous(limits = c(map_layers$plot.boundary$y)) +
      ggtitle(label = paste0("Earliest sunset (", min(dates), " to ", max(dates), ")")) +
      theme_bw() +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
  )
  dev.off()


  sunrise_range_centroid |>
    as.data.frame() |>
    dplyr::select(STATIONID, 
                  LONGITUDE = X, 
                  LATITUDE = Y, 
                  MAX_SUNRISE = max_sunrise, 
                  MAX_SUNRISE_TIME = label_sunrise,
                  MIN_SUNSET = min_sunset,
                  MIN_SUNRISE_TIME = label_sunset) |>
    dplyr::arrange(STATIONID) |>
    write.csv(file = here::here(paste0("latest_sunrise_", min(dates), "_", max(dates), ".csv")),
              row.names = FALSE)

}


### ----

library(akgfmaps)
library(fishmethods)


# Leg 1
make_sunrise_table(min_date = "2023-05-20",
                   max_date = "2023-06-09",
                   region = "ebs",
                   timezone = "US/Alaska")

# Leg 2
make_sunrise_table(min_date = "2023-06-11",
                   max_date = "2023-07-03",
                   region = "ebs",
                   timezone = "US/Alaska")

# Leg 3
make_sunrise_table(min_date = "2023-07-05",
                   max_date = "2023-07-28",
                   region = "ebs",
                   timezone = "US/Alaska")

# Leg 4
make_sunrise_table(min_date = "2023-07-30",
                   max_date = "2023-08-21",
                   region = "ebs",
                   timezone = "US/Alaska")