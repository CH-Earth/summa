library(ncdf4)
library(hydroGOF)

source("utils/plotting/plot_utils.R")

quartz(width = 10, height = 6)
#pdf("compare_coupled_sequential_mizuroute.pdf", width = 10, height = 6)

# ----- file paths/names -----

data_path <- "~/data/great-slave-lake/Athabasca-river"

mizu_path <- file.path(data_path, "results/mizuRoute")
mizu_file <- file.path(mizu_path, "run_1_standalone.h.2014-09-01-03600.nc")

summa_path <- file.path(data_path, "results/SUMMA")
summa_file <- file.path(summa_path, "run_1_coupled_timestep.nc")

cat("mizuRoute file:", mizu_file, "\n")
cat("SUMMA file:    ", summa_file, "\n")


# ----- sequential SUMMA-mizuRoute -----

nc <- nc_open(mizu_file)

  segid_mizu   <- ncvar_get(nc, "reachID")
  time_mizu    <- ncvar_get(nc, "time")
  units_mizu   <- ncatt_get(nc, "time", "units")$value
  qrouted_mizu <- ncvar_get(nc, "KWroutedRunoff")

nc_close(nc)


# ----- coupled SUMMA-mizuRoute -----

nc <- nc_open(summa_file)

  segid_summa   <- ncvar_get(nc, "seg")
  time_summa    <- ncvar_get(nc, "time")
  units_summa   <- ncatt_get(nc, "time", "units")$value
  uparea_summa  <- ncvar_get(nc, "upArea")
  qrouted_summa <- ncvar_get(nc, "q_reach")

nc_close(nc)


# ----- time -----

date_mizu  <- nc_time(time_mizu, units_mizu)
date_summa <- nc_time(time_summa, units_summa)


# ----- choose most downstream reach -----

iRch_summa <- which.max(uparea_summa)
reach_id   <- segid_summa[iRch_summa]

iRch_mizu <- which(segid_mizu == reach_id)

cat("Reach:", reach_id, "\n")
cat("Upstream area:", uparea_summa[iRch_summa] / 1e6, "km2\n")


# ----- extract routed streamflow -----

routed_mizu  <- qrouted_mizu[iRch_mizu, ]
routed_summa <- qrouted_summa[iRch_summa, ]


# ----- comparison period -----

date_range <- as.POSIXct(c("2016-07-01", "2024-01-01"), tz = "UTC")
ylim       <- c(0, 2.5)

tol <- 1

ix_mizu <- as.numeric(date_mizu) >= as.numeric(date_range[1]) - tol &
           as.numeric(date_mizu) <= as.numeric(date_range[2]) + tol

ix_summa <- as.numeric(date_summa) >= as.numeric(date_range[1]) - tol &
            as.numeric(date_summa) <= as.numeric(date_range[2]) + tol

date_mizu_sub  <- date_mizu[ix_mizu]
date_summa_sub <- date_summa[ix_summa]

dt <- abs(as.numeric(date_summa_sub) - as.numeric(date_mizu_sub))
stopifnot(max(dt) < tol)

# ----- convert streamflow to runoff depth -----

area <- uparea_summa[iRch_summa]

conv <- 86400 * 1000 / area

sim <- routed_summa[ix_summa] * conv
ref <- routed_mizu[ix_mizu]   * conv


# ----- performance metrics -----

kge <- KGE(sim = sim, obs = ref)
nse <- NSE(sim = sim, obs = ref)

cat(sprintf("Reach: %d  KGE: %.6f  NSE: %.6f\n",
            reach_id, kge, nse))


# ----- plot -----

plot(date_mizu_sub, ref,
     type = "n",
     xlim = date_range,
     ylim = ylim,
     xlab = "Date",
     ylab = "Runoff (mm/day)",
     main = sprintf("Reach %s: KGE=%.3f, NSE=%.3f",
                    reach_id, kge, nse))

lines(date_mizu_sub,  ref, col = "darkblue",  lwd = 2)
lines(date_summa_sub, sim, col = "lightblue", lwd = 1, lty = 2)

legend("topleft",
       legend = c("Sequential", "Coupled"),
       col = c("darkblue", "lightblue"),
       lwd = c(2, 1),
       lty = c(1, 2))

#dev.off()

