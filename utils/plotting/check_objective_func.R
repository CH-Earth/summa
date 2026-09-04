library(ncdf4)
library(hydroGOF)

# files
sim_file <- "~/data/century/test/summa_results/run1_test_timestep.nc"
obs_file <- "~/data/century/test/mizuroute_input/CAN_05BB001_daily_flow_observations.nc"

# evaluation period
start_date <- as.POSIXct("1982-10-01", tz="UTC")
end_date   <- as.POSIXct("1983-09-30 23:59:59", tz="UTC")

# ------------------------------------------------------------
# observations
# ------------------------------------------------------------

nc_obs <- nc_open(obs_file)

time_obs <- ncvar_get(nc_obs, "time")
q_obs    <- ncvar_get(nc_obs, "q_obs")

nc_close(nc_obs)

date_obs <- as.POSIXct("1950-01-01", tz="UTC") + time_obs * 60

# ------------------------------------------------------------
# simulation
# ------------------------------------------------------------

nc_sim <- nc_open(sim_file)

time_sim <- ncvar_get(nc_sim, "time")
q_reach  <- ncvar_get(nc_sim, "Q_reach")
up_area  <- ncvar_get(nc_sim, "upArea")
seg      <- ncvar_get(nc_sim, "seg")

time_units <- ncatt_get(nc_sim, "time", "units")$value

nc_close(nc_sim)

# outlet = segment with largest upstream drainage area
i_seg <- which.max(up_area)

q_sim <- q_reach[i_seg, ]

# SUMMA time is seconds since 1990-01-01
date_sim <- as.POSIXct("1990-01-01", tz="UTC") + time_sim

# ------------------------------------------------------------
# aggregate hourly simulation to daily period-ending means
# exactly as done in Fortran:
#
#     time_obs - 1 day < time_sim <= time_obs
# ------------------------------------------------------------

sim_daily <- data.frame(
  time  = date_obs,
  q_sim = NA_real_
)

for(i in seq_along(date_obs)) {

  t_end   <- date_obs[i]
  t_start <- t_end - 86400

  ix <- date_sim > t_start & date_sim <= t_end

  if(any(ix)) {
    sim_daily$q_sim[i] <- mean(q_sim[ix], na.rm=TRUE)
  }
}

names(sim_daily)[2] <- "q_sim"

# ------------------------------------------------------------
# merge model simulations and observations 
# ------------------------------------------------------------

obs <- data.frame(
  time=date_obs,
  q_obs=q_obs
)

dat <- merge(obs, sim_daily, by="time")

# evaluation period
dat <- subset(
  dat,
  time >= start_date &
  time <= end_date
)

# remove missing values
dat <- dat[
  is.finite(dat$q_obs) &
  is.finite(dat$q_sim),
]

# ------------------------------------------------------------
# calculate performance metrics 
# ------------------------------------------------------------

kge <- KGE(dat$q_sim, dat$q_obs)
nse <- NSE(dat$q_sim, dat$q_obs)
rmse <- sqrt(mean((dat$q_obs - dat$q_sim)^2))
mae  <- mean(abs(dat$q_obs - dat$q_sim))

cat("\nR\n")
cat("KGE  =", kge, "\n")
cat("NSE  =", nse, "\n")
cat("RMSE =", rmse, "\n")
cat("MAE  =", mae, "\n")

# ------------------------------------------------------------
# plot 
# ------------------------------------------------------------

# read flows aligned by Fortran
x <- read.csv("~/data/century/test/aligned_flow.txt")
x$date <- as.POSIXct("1950-01-01", tz="UTC") + x$time * 60

kge <- KGE(x$flowSim,  x$flowObs)
nse <- NSE(x$flowSim,  x$flowObs)
rmse <- sqrt(mean((x$flowSim - x$flowObs)^2))
mae  <- mean(abs(x$flowSim - x$flowObs))

cat("\nFortran\n")
cat("KGE  =", kge, "\n")
cat("NSE  =", nse, "\n")
cat("RMSE =", rmse, "\n")
cat("MAE  =", mae, "\n")



plot(
  dat$time,
  dat$q_obs,
  type="n",
  xlab="Date",
  ylab=expression(Streamflow~(m^3/s))
)

# R alignment
lines(dat$time, dat$q_obs, col="darkblue")
lines(dat$time, dat$q_sim, col="lightblue")

# Fortran alignment
lines(x$date, x$flowObs, col="red")
lines(x$date, x$flowSim, col="orange")

legend(
  "topright",
  c("Observed (R)", "Simulated (R)",
    "Observed (Fortran)", "Simulated (Fortran)"),
  col=c("darkblue", "lightblue", "red", "orange"),
  lty=1
)
