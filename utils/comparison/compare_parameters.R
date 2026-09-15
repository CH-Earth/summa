library(ncdf4)

options(scipen = 999)
options(width = 200)

data_path <- '~/data/century'

exp1_path <- file.path(data_path, 'test/summa_results')
exp2_path <- file.path(data_path, 'old/CAN_05BB001/7-FA_mod_IC_newSUMMA_inf_GA_bsflwParams_MP/best_run/summa/summa_results')

exp1_file <- file.path(exp1_path, 'run1_coupled_timestep.nc')
exp2_file <- file.path(exp2_path, 'run1_timestep.nc')

nc1 <- nc_open(exp1_file)
nc2 <- nc_open(exp2_file)

# Variables present in either file
vars <- union(names(nc1$var), names(nc2$var))

# Exclude dimensions time and seg
has_excluded_dim <- function(nc, var) {
  if (!var %in% names(nc$var)) return(FALSE)

  dims <- sapply(nc$var[[var]]$dim, function(x) x$name)

  any(dims %in% c("time", "seg"))
}

# Keep variables that do not have a time or segment dimension in either file
vars_keep <- vars[
  !sapply(vars, function(v)
    has_excluded_dim(nc1, v) || has_excluded_dim(nc2, v)
  )
]

# Compare
out <- lapply(vars_keep, function(v) {

  in1 <- v %in% names(nc1$var)
  in2 <- v %in% names(nc2$var)

  x1 <- if (in1) as.vector(ncvar_get(nc1, v)) else NA
  x2 <- if (in2) as.vector(ncvar_get(nc2, v)) else NA

  n <- max(length(x1), length(x2))

  if (length(x1) == 1 && n > 1) x1 <- rep(x1, n)
  if (length(x2) == 1 && n > 1) x2 <- rep(x2, n)

  # If variable is absent from one file
  if (!in1) x1 <- rep(NA, n)
  if (!in2) x2 <- rep(NA, n)

  data.frame(
    parameter = v,
    index     = seq_len(n),
    file1     = x1,
    file2     = x2,
    difference = x2 - x1,
    identical = x1 == x2
  )
})

param_compare <- do.call(rbind, out)

nc_close(nc1)
nc_close(nc2)

# Show everything
print(param_compare, row.names = FALSE)

# Or just parameters that differ
param_compare[
  is.na(param_compare$identical) | !param_compare$identical,
]
