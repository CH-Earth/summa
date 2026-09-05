library(ncdf4)

# ------------------------------------------------------------------------------
# Test case
# ------------------------------------------------------------------------------

test_root <- path.expand("~/models/summa_test_cases/test_cases/output")

experiment_dir  <- "reynolds_canopySrad_windPrfile_stomResist"
#experiment_name <- "reynoldsUEB2stream"
experiment_name <- "reynoldsExponential"

experiment_dir  <- "reynolds_groundwatr"
experiment_name <- "reynoldsLumpedQTopmodel"

exe0 <- "_exe_0"
exe1 <- "_exe_1"

# ------------------------------------------------------------------------------
# Find output files
# ------------------------------------------------------------------------------

output_dir <- file.path(test_root, experiment_dir)

file0 <- Sys.glob(
  file.path(output_dir, paste0(experiment_name, exe0, "*.nc"))
)

file1 <- Sys.glob(
  file.path(output_dir, paste0(experiment_name, exe1, "*.nc"))
)

if (length(file0) != 1) {
  stop("Expected exactly one exe_0 file; found ", length(file0))
}

if (length(file1) != 1) {
  stop("Expected exactly one exe_1 file; found ", length(file1))
}

cat("Executable 0:\n  ", file0, "\n")
cat("Executable 1:\n  ", file1, "\n")

# ------------------------------------------------------------------------------
# Read NetCDF
# ------------------------------------------------------------------------------

nc0 <- nc_open(file0)
nc1 <- nc_open(file1)

time0 <- ncvar_get(nc0, "time")
time1 <- ncvar_get(nc1, "time")

# Convert model time to days from start
time0 <- (time0 - time0[1]) / 86400
time1 <- (time1 - time1[1]) / 86400


# ------------------------------------------------------------------------------
# State variables
# ------------------------------------------------------------------------------

swe0 <- ncvar_get(nc0, "scalarSWE")
swe1 <- ncvar_get(nc1, "scalarSWE")

surface_temp0 <- ncvar_get(nc0, "scalarSurfaceTemp")
surface_temp1 <- ncvar_get(nc1, "scalarSurfaceTemp")

root_temp0 <- ncvar_get(nc0, "scalarRootZoneTemp")
root_temp1 <- ncvar_get(nc1, "scalarRootZoneTemp")

soil_water0 <- ncvar_get(nc0, "scalarTotalSoilWat")
soil_water1 <- ncvar_get(nc1, "scalarTotalSoilWat")


# ------------------------------------------------------------------------------
# Fluxes
# ------------------------------------------------------------------------------

snow_drain0 <- ncvar_get(nc0, "scalarSnowDrainage")
snow_drain1 <- ncvar_get(nc1, "scalarSnowDrainage")

transpire0 <- ncvar_get(nc0, "scalarCanopyTranspiration")
transpire1 <- ncvar_get(nc1, "scalarCanopyTranspiration")


# ------------------------------------------------------------------------------
# Cumulative fluxes
#
# scalarSnowDrainage:         m s-1
# scalarCanopyTranspiration:  kg m-2 s-1
#
# 1 m water       = 1000 mm
# 1 kg m-2 water  = 1 mm
# ------------------------------------------------------------------------------

dt0 <- c(0, diff(ncvar_get(nc0, "time")))
dt1 <- c(0, diff(ncvar_get(nc1, "time")))

cum_snow_drain0 <- cumsum(snow_drain0 * dt0 * 1000)
cum_snow_drain1 <- cumsum(snow_drain1 * dt1 * 1000)

cum_transpire0 <- cumsum(transpire0 * dt0)
cum_transpire1 <- cumsum(transpire1 * dt1)


# ------------------------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------------------------

plot_compare <- function(time0, x0, time1, x1, ylab, main) {

  ylim <- range(c(x0, x1), finite = TRUE)

  plot(time0, x0,
       type = "n",
       xlab = "Time (days)",
       ylab = ylab,
       main = main,
       ylim = ylim)

  lines(time1, x1, col="darkblue",  lty = 1)
  lines(time1, x1, col="lightblue", lty = 2)

  legend("topright",
         legend = c("exe_0", "exe_1"),
         col = c("darkblue","lightblue"),
         lty = c(1, 2),
         bty = "n")
}


# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

quartz(width = 9, height = 6)

par(mfrow = c(3, 2),
    mar = c(4, 4.5, 2.5, 1))

plot_compare(time0, swe0,
             time1, swe1,
             "SWE (kg m-2)",
             "Snow water equivalent")

plot_compare(time0, surface_temp0,
             time1, surface_temp1,
             "Temperature (K)",
             "Surface temperature")

plot_compare(time0, root_temp0,
             time1, root_temp1,
             "Temperature (K)",
             "Root-zone temperature")

plot_compare(time0, soil_water0,
             time1, soil_water1,
             "Water storage (kg m-2)",
             "Total soil water")

plot_compare(time0, cum_snow_drain0,
             time1, cum_snow_drain1,
             "Cumulative drainage (mm)",
             "Cumulative snow drainage")

plot_compare(time0, cum_transpire0,
             time1, cum_transpire1,
             "Cumulative transpiration (mm)",
             "Cumulative transpiration")

# ------------------------------------------------------------------------------
# Plot differences
# ------------------------------------------------------------------------------

plot_difference <- function(time, x0, x1, ylab, main) {

  diff <- x1 - x0

  ylim <- range(diff, finite = TRUE)

  plot(time, diff,
       type = "l",
       xlab = "Time (days)",
       ylab = ylab,
       main = main,
       ylim = ylim)

  abline(h = 0, lty = 2)
}


quartz(width = 9, height = 6)

par(mfrow = c(3, 2),
    mar = c(4, 4.5, 2.5, 1))

plot_difference(
  time0, swe0, swe1,
  "Difference (kg m-2)",
  "Snow water equivalent"
)

plot_difference(
  time0, surface_temp0, surface_temp1,
  "Difference (K)",
  "Surface temperature"
)

plot_difference(
  time0, root_temp0, root_temp1,
  "Difference (K)",
  "Root-zone temperature"
)

plot_difference(
  time0, soil_water0, soil_water1,
  "Difference (kg m-2)",
  "Total soil water"
)

plot_difference(
  time0, cum_snow_drain0, cum_snow_drain1,
  "Difference (mm)",
  "Cumulative snow drainage"
)

plot_difference(
  time0, cum_transpire0, cum_transpire1,
  "Difference (mm)",
  "Cumulative transpiration"
)

# ------------------------------------------------------------------------------
# Clean up
# ------------------------------------------------------------------------------

nc_close(nc0)
nc_close(nc1)




