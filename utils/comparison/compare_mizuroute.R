library(ncdf4)
library(hydroGOF)

source("utils/plotting/plot_utils.R")

quartz(width = 10, height = 6)
#pdf("compare_coupled_sequential_mizuroute.pdf", width = 10, height = 6)

# ----- file paths/names -----

data_path <- "~/data/century/test"

mizu_path <- file.path(data_path, "mizuroute_results")

gru_file <- file.path(
  mizu_path,
  "run1_gru_remap.h.1981-10-01-00000.nc"
)

hru_file <- file.path(
  mizu_path,
  "run1_hru_direct.h.1981-10-01-00000.nc"
)

summa_file <- file.path(
  data_path,
  "summa_results/run1_coupled_timestep.nc"
)

cat("GRU remap file:", gru_file, "\n")
cat("HRU direct file:", hru_file, "\n")
cat("Coupled file:   ", summa_file, "\n")


# ----- standalone mizuRoute: GRU runoff + remapping -----

nc <- nc_open(gru_file)

  segid_gru   <- ncvar_get(nc, "reachID")
  time_gru    <- ncvar_get(nc, "time")
  units_gru   <- ncatt_get(nc, "time", "units")$value
  qrouted_gru <- ncvar_get(nc, "KWroutedRunoff")

nc_close(nc)


# ----- standalone mizuRoute: HRU runoff directly -----

nc <- nc_open(hru_file)

  segid_hru   <- ncvar_get(nc, "reachID")
  time_hru    <- ncvar_get(nc, "time")
  units_hru   <- ncatt_get(nc, "time", "units")$value
  qrouted_hru <- ncvar_get(nc, "KWroutedRunoff")

nc_close(nc)


# ----- coupled SUMMA-mizuRoute -----

nc <- nc_open(summa_file)

  segid_summa   <- ncvar_get(nc, "seg")
  time_summa    <- ncvar_get(nc, "time")
  units_summa   <- ncatt_get(nc, "time", "units")$value
  uparea_summa  <- ncvar_get(nc, "upArea")
  qrouted_summa <- ncvar_get(nc, "Q_reach")

nc_close(nc)


# ----- time -----

date_gru   <- nc_time(time_gru,   units_gru)
date_hru   <- nc_time(time_hru,   units_hru)
date_summa <- nc_time(time_summa, units_summa)


# ----- choose most downstream reach -----

iRch_summa <- which.max(uparea_summa)
reach_id   <- segid_summa[iRch_summa]

iRch_gru <- which(segid_gru == reach_id)
iRch_hru <- which(segid_hru == reach_id)

stopifnot(length(iRch_gru) == 1)
stopifnot(length(iRch_hru) == 1)

cat("Reach:", reach_id, "\n")
cat("Upstream area:", uparea_summa[iRch_summa] / 1e6, "km2\n")


# ----- extract routed streamflow -----

routed_gru   <- qrouted_gru[iRch_gru, ]
routed_hru   <- qrouted_hru[iRch_hru, ]
routed_summa <- qrouted_summa[iRch_summa, ]


# ----- comparison period -----

# use standalone mizuRoute period, skipping first year as spinup
date_range <- range(date_gru)
date_range[1] <- seq(date_range[1], by = "1 year", length.out = 2)[2]

tol  <- 1

ix_gru <- as.numeric(date_gru) >= as.numeric(date_range[1]) - tol &
          as.numeric(date_gru) <= as.numeric(date_range[2]) + tol

ix_hru <- as.numeric(date_hru) >= as.numeric(date_range[1]) - tol &
          as.numeric(date_hru) <= as.numeric(date_range[2]) + tol

ix_summa <- as.numeric(date_summa) >= as.numeric(date_range[1]) - tol &
            as.numeric(date_summa) <= as.numeric(date_range[2]) + tol

date_gru_sub   <- date_gru[ix_gru]
date_hru_sub   <- date_hru[ix_hru]
date_summa_sub <- date_summa[ix_summa]

stopifnot(length(date_gru_sub) == length(date_hru_sub))
stopifnot(length(date_gru_sub) == length(date_summa_sub))

dt_gru_hru <- abs(as.numeric(date_gru_sub) - as.numeric(date_hru_sub))
dt_gru_sum <- abs(as.numeric(date_gru_sub) - as.numeric(date_summa_sub))

stopifnot(max(dt_gru_hru) <= tol)
stopifnot(max(dt_gru_sum) <= tol)


# ----- convert streamflow to runoff depth -----

area <- uparea_summa[iRch_summa]
conv <- 86400 * 1000 / area

ref_gru <- routed_gru[ix_gru]     * conv
sim_hru <- routed_hru[ix_hru]     * conv
sim_sum <- routed_summa[ix_summa] * conv


# ----- performance metrics -----

kge_hru <- KGE(sim = sim_hru, obs = ref_gru)
nse_hru <- NSE(sim = sim_hru, obs = ref_gru)

kge_sum <- KGE(sim = sim_sum, obs = ref_gru)
nse_sum <- NSE(sim = sim_sum, obs = ref_gru)

cat(sprintf(
  "HRU direct vs GRU remap: KGE = %.8f  NSE = %.8f\n",
  kge_hru, nse_hru
))

cat(sprintf(
  "Coupled vs GRU remap:    KGE = %.8f  NSE = %.8f\n",
  kge_sum, nse_sum
))


# ----- plot -----

ylim <- c(0,10)

plot(date_gru_sub, ref_gru,
     type = "n",
     xlim = date_range,
     ylim = ylim,
     xlab = "Date",
     ylab = "Runoff (mm/day)",
     main = sprintf("Reach %s", reach_id))

lines(date_gru_sub, ref_gru,
      col = "seagreen", lwd = 2)

lines(date_hru_sub, sim_hru,
      col = "darkblue", lwd = 1, lty = 2)

lines(date_summa_sub, sim_sum,
      col = "lightblue", lwd = 1, lty = 3)

legend("topleft",
       legend = c("GRU remap", "HRU direct", "Coupled"),
       col = c("seagreen", "darkblue", "lightblue"),
       lwd = c(2, 1, 1),
       lty = c(1, 2, 3))

# ----- performance metrics -----

compare_stats <- function(x, y) {
  c(
    KGE     = KGE(sim = x, obs = y),
    NSE     = NSE(sim = x, obs = y),
    MAE     = mean(abs(x - y), na.rm = TRUE),
    MaxAE   = max(abs(x - y), na.rm = TRUE)
  )
}

stats <- rbind(
  "HRU direct vs GRU remap" = compare_stats(sim_hru, ref_gru),
  "Coupled vs GRU remap"    = compare_stats(sim_sum, ref_gru),
  "Coupled vs HRU direct"   = compare_stats(sim_sum, sim_hru)
)

print(stats, digits = 8)

#dev.off()
