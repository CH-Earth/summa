#!/usr/bin/env Rscript

library(ncdf4)

# ==============================================================================
# User settings
# ==============================================================================

ncfile1 <- "~/data/great-slave-lake/Athabasca-river/results/SUMMA/run_1_stable_serial_timestep.nc"

#ncfile2 <- "~/data/great-slave-lake/Athabasca-river/results/SUMMA/run_1_dev_serial_timestep.nc"
ncfile2 <- "~/data/great-slave-lake/Athabasca-river/results/SUMMA/run_1_dev_np16_G01-05_timestep.nc"

label1 <- "stable"
label2 <- "development"

hru_id <- 16

variables <- c(
  "pptrate",
  "airtemp",
  "scalarCanopyTemp",
  "scalarSurfaceTemp",
  "scalarRootZoneTemp",
  "scalarCanopyWat",
  "scalarSWE",
  "scalarTotalSoilWat"
)

# ==============================================================================
# Helper: read one HRU variable
# ==============================================================================

get_hru_var <- function(nc, var, hru_index) {

  dims <- sapply(
    nc$var[[var]]$dim,
    function(x) x$name
  )

  start <- rep(1, length(dims))
  count <- sapply(
    nc$var[[var]]$dim,
    function(x) x$len
  )

  hru_dim <- which(dims == "hru")

  start[hru_dim] <- hru_index
  count[hru_dim] <- 1

  drop(
    ncvar_get(
      nc,
      var,
      start = start,
      count = count
    )
  )
}

# ==============================================================================
# Open NetCDF files
# ==============================================================================

nc1 <- nc_open(ncfile1)
nc2 <- nc_open(ncfile2)

# ==============================================================================
# Identify HRU in each file
# ==============================================================================

hru_ids1 <- ncvar_get(nc1, "hruId")
hru_ids2 <- ncvar_get(nc2, "hruId")

hru_index1 <- match(hru_id, hru_ids1)
hru_index2 <- match(hru_id, hru_ids2)

if (is.na(hru_index1)) {
  stop("hruId ", hru_id, " not found in file 1")
}

if (is.na(hru_index2)) {
  stop("hruId ", hru_id, " not found in file 2")
}

cat("hruId:        ", hru_id, "\n")
cat("File 1 index: ", hru_index1, "\n")
cat("File 2 index: ", hru_index2, "\n")

# ==============================================================================
# Read time
# ==============================================================================

time1 <- ncvar_get(nc1, "time")
time2 <- ncvar_get(nc2, "time")

time_units1 <- ncatt_get(nc1, "time", "units")$value
time_units2 <- ncatt_get(nc2, "time", "units")$value

origin1 <- sub("^seconds since ", "", time_units1)
origin1 <- sub(" -0:00$", "", origin1)

origin2 <- sub("^seconds since ", "", time_units2)
origin2 <- sub(" -0:00$", "", origin2)

datetime1 <- as.POSIXct(origin1, tz = "UTC") + time1
datetime2 <- as.POSIXct(origin2, tz = "UTC") + time2

# ==============================================================================
# Set up graphics window
# ==============================================================================

quartz(width = 12, height = 6)

par(
  mfrow = c(2, 4),
  mar   = c(4, 4, 3, 1),
  oma   = c(0, 0, 2, 0)
)

# ==============================================================================
# Plot variables
# ==============================================================================

for (var in variables) {

  if (!var %in% names(nc1$var)) {
    warning("Variable not found in file 1: ", var)
    next
  }

  if (!var %in% names(nc2$var)) {
    warning("Variable not found in file 2: ", var)
    next
  }

  x1 <- get_hru_var(nc1, var, hru_index1)
  x2 <- get_hru_var(nc2, var, hru_index2)

  units <- ncatt_get(nc1, var, "units")$value

  ylim <- range(c(x1, x2), na.rm = TRUE)

  plot(
    datetime1,
    x1,
    type = "l",
    xlab = "",
    ylab = units,
    main = var,
    ylim = ylim
  )

  lines(
    datetime2,
    x2,
    col = "red",
    lty = 2
  )

  legend(
    "topright",
    legend = c(label1, label2),
    col = c("black", "red"),
    lty = c(1, 2),
    bty = "n",
    cex = 0.7
  )
}

mtext(
  paste("SUMMA science check — hruId", hru_id),
  outer = TRUE,
  cex = 1.1
)

# ==============================================================================
# Close NetCDF files
# ==============================================================================

nc_close(nc1)
nc_close(nc2)
