# ----- convert NetCDF time -----

nc_time <- function(time, units) {

  parts  <- strsplit(units, " since ", fixed = TRUE)[[1]]
  unit   <- parts[1]
  origin <- as.POSIXct(parts[2], tz = "UTC")

  scale <- switch(unit,
                  "seconds" = 1,
                  "minutes" = 60,
                  "hours"   = 3600,
                  "days"    = 86400,
                  stop(paste("Unsupported time unit:", unit)))

  origin + time * scale
}
