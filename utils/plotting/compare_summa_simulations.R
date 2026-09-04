library(ncdf4)
library(hydroGOF)

source("utils/plotting/plot_utils.R")

quartz(width = 10, height = 6)

# ----- files -----

ref_path <- "~/data/great-slave-lake/Athabasca-river/orig"
new_path <- "~/data/great-slave-lake/Athabasca-river"

ref_file <- file.path(ref_path, "results/SUMMA/run_1_timestep.nc")
new_file <- file.path(new_path, "results/SUMMA/run_1_coupled_timestep.nc")


# ----- reference SUMMA -----

nc <- nc_open(ref_file)

  time_ref   <- ncvar_get(nc, "time")
  units_ref  <- ncatt_get(nc, "time", "units")$value
  gruid_ref  <- ncvar_get(nc, "gruId")
  runoff_ref <- ncvar_get(nc, "averageRoutedRunoff")

nc_close(nc)


# ----- coupled SUMMA -----

nc <- nc_open(new_file)

  time_new   <- ncvar_get(nc, "time")
  units_new  <- ncatt_get(nc, "time", "units")$value
  gruid_new  <- ncvar_get(nc, "gruId")
  runoff_new <- ncvar_get(nc, "averageRoutedRunoff")

nc_close(nc)


# ----- time -----

date_ref <- nc_time(time_ref, units_ref)
date_new <- nc_time(time_new, units_new)


# ----- plot -----

par(mfrow = c(4,3), mar = c(3,3,2,1), oma = c(0,3,0,0))

date_range <- as.POSIXct(c("2016-07-01", "2024-01-01"), tz = "UTC")
ylim       <- c(0, 10)

# indices within comparison period
ix_ref <- date_ref >= date_range[1] & date_ref <= date_range[2]
ix_new <- date_new >= date_range[1] & date_new <= date_range[2]

for(iGRU in 1:12) {

  gru_id <- gruid_ref[iGRU]
  iNew   <- which(gruid_new == gru_id)

  conv <- 86400 * 1000
  sim  <- runoff_new[iNew, ix_new]*conv
  obs  <- runoff_ref[iGRU, ix_ref]*conv

  kge <- KGE(sim = sim, obs = obs)
  nse <- NSE(sim = sim, obs = obs)

  cat(sprintf("GRU: %d  KGE: %.6f  NSE: %.6f\n",
               gru_id, kge, nse))

  plot(date_ref[ix_ref], obs,
       type = "n",
       xlim = date_range,
       ylim = ylim,
       xlab = "",
       ylab = "",
       main = sprintf("GRU %s: KGE=%.3f, NSE=%.3f",
                      gru_id, kge, nse))

  lines(date_ref[ix_ref], obs, col = "darkblue", lwd=2)
  lines(date_new[ix_new], sim, col = "lightblue", lty = 2)

  legend("topleft",
         legend = c("ref", "new"),
         col = c("darkblue", "lightblue"),
         lwd = c(2, 1),
         lty = c(1, 2))
}

mtext("Routed runoff (mm/day)", side = 2, outer = TRUE, line = 1)
