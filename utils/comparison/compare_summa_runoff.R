library(ncdf4)
library(hydroGOF)

source("utils/plotting/plot_utils.R")

# ----- git commit date -----

git_commit_date <- function(hash) {
  system2(
    "git",
    c("show", "-s", "--format=%ad", "--date=short", hash),
    stdout = TRUE
  )
}

# -----------------------------------------------------------------------------
# ----- COMPARE FILES ---------------------------------------------------------
# -----------------------------------------------------------------------------

compare_runs <- function(ref_file, new_file) {

  # ----- reference SUMMA -----
  
  nc <- nc_open(ref_file)
  
  time_ref   <- ncvar_get(nc, "time")
  units_ref  <- ncatt_get(nc, "time", "units")$value
  runoff_ref <- ncvar_get(nc, "averageRoutedRunoff")
  
  hash_ref    <- ncatt_get(nc, 0, "gitHash")$value
  branch_ref  <- ncatt_get(nc, 0, "gitBranch")$value
  
  nc_close(nc)
  
  date_ref_commit <- git_commit_date(hash_ref)
  
  cat("\nReference:\n")
  cat("  file:      ", ref_file, "\n")
  cat("  branch:    ", branch_ref, "\n")
  cat("  hash:      ", hash_ref, "\n")
  cat("  date commit", date_ref_commit, "\n")
  
  # ----- new SUMMA -----
  
  nc <- nc_open(new_file)
  
  time_new   <- ncvar_get(nc, "time")
  units_new  <- ncatt_get(nc, "time", "units")$value
  runoff_new <- ncvar_get(nc, "averageRoutedRunoff")
  
  hash_new    <- ncatt_get(nc, 0, "gitHash")$value
  branch_new  <- ncatt_get(nc, 0, "gitBranch")$value
  
  nc_close(nc)
  
  date_new_commit <- git_commit_date(hash_new)
  
  cat("\nNew:\n")
  cat("  file:      ", new_file, "\n")
  cat("  branch:    ", branch_new, "\n")
  cat("  hash:      ", hash_new, "\n")
  cat("  date commit", date_new_commit, "\n")
  
  
  # ----- time -----
  
  date_ref <- nc_time(time_ref, units_ref)
  date_new <- nc_time(time_new, units_new)
  
  start_date <- max(min(date_ref), min(date_new))
  end_date   <- min(max(date_ref), max(date_new))
  
  # (skip one year spinup)
  start_date <- seq(start_date, by = "1 year", length.out = 2)[2]
  
  ix_ref <- date_ref >= start_date & date_ref <= end_date
  ix_new <- date_new >= start_date & date_new <= end_date
  
  date_ref_sub <- date_ref[ix_ref]
  date_new_sub <- date_new[ix_new]
  
  stopifnot(length(date_ref_sub) == length(date_new_sub))
  
  dt <- abs(as.numeric(date_ref_sub) - as.numeric(date_new_sub))
  stopifnot(max(dt) < 1)
  
  # ----- runoff -----
  
  conv <- 86400 * 1000       # m/s -> mm/day
  
  ref <- as.numeric(runoff_ref[ix_ref]) * conv
  sim <- as.numeric(runoff_new[ix_new]) * conv
  
  # ----- diagnostics -----
  
  diff <- sim - ref
  
  kge <- KGE(sim = sim, obs = ref)
  nse <- NSE(sim = sim, obs = ref)
  
  cat("\nStatistics:\n")
  cat(sprintf("KGE:               %.8f\n", kge))
  cat(sprintf("NSE:               %.8f\n", nse))
  cat(sprintf("Mean ref:          %.8f mm/day\n", mean(ref)))
  cat(sprintf("Mean new:          %.8f mm/day\n", mean(sim)))
  cat(sprintf("Mean difference:   %.8f mm/day\n", mean(diff)))
  cat(sprintf("Mean abs diff:     %.8f mm/day\n", mean(abs(diff))))
  cat(sprintf("Max abs diff:      %.8f mm/day\n", max(abs(diff))))
  cat(sprintf("Correlation:       %.8f\n", cor(ref, sim)))

} # end function compare_files

# -----------------------------------------------------------------------------
# ----- LOOP THROUGH DIFFERENT OUTPUT FILES -----------------------------------
# -----------------------------------------------------------------------------

data_path <- '~/data/century'

ref_path <- file.path(data_path, 'old/CAN_05BB001/7-FA_mod_IC_newSUMMA_inf_GA_bsflwParams_MP/best_run/summa/summa_results')
ref_file <- file.path(ref_path, 'run1_timestep.nc')

new_path <- file.path(data_path, 'test/summa_results')

files    <- file.path(
  new_path,
  c(
    'run1_10Jul2025_timestep.nc', 
    'run1_26Aug2025_timestep.nc', 
    'run1_13Dec2025_timestep.nc', 
    'run1_09Jan2026_timestep.nc', 
    'run1_03Feb2026_timestep.nc', 
    'run1_stable_timestep.nc'
  )
)

for (new_file in files) {

 cat("\n")
 cat(strrep("-", 50), "\n")
 cat(strrep("-", 50), "\n")

 compare_runs(ref_file, new_file)

 ref_file <- file.path(new_path, 'run1_10Jul2025_timestep.nc')

}

ref_path <- file.path(data_path, 'test/summa_results')
ref_file <- file.path(ref_path,  'run1_stable_timestep.nc')

for (new_file in files) {

 cat("\n")
 cat(strrep("-", 50), "\n")
 cat(strrep("-", 50), "\n")

 compare_runs(ref_file, new_file)

}











