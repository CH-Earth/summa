# Lookup table provenance
This information cannot easily be encoded in the files themselves, due to the way files are being read.

#### GENPARM
- 2021-08-10: Added NOAH-MP default parameters (NCAR Research Applications Laboratory, RAL, 2021)

#### MPTABLE
- 2021-08-10: Added NOAH-MP default parameters (NCAR RAL, 2021)

#### SOILPARM
- 2021-08-10: Added NOAH-MP default parameters (NCAR RAL, 2021)
- 2021-08-10: Added ROSETTA table values (see `ROSETTA soil parameter values`)

#### VEGPARM
- 2021-08-10: Added NOAH-MP default parameters (NCAR RAL, 2021)


## ROSETTA soil parameter values 
ROSETTA soil parameters are a combination of the NOAH-MP STAS and STAS-RUC tables (NCAR RAL, 2021), updated with five additional parameters derived by the ROSETTA model (U.S. Department of Agriculture: Agricultural Research Service (USDA ARS), 2021).

Columns `BB`, `DRYSMC`, `HC`, `MAXSMC`, `SATPIS`, `SATDK`, `SATDW`, `WLTSMC` and `QTZ` are duplicated from the STAS-RUC table.

Column `REFSMC` is duplicated from the STAS table.

Columns `theta_res`, `theta_sat`, `vGn_alpha`, `vGn_n` and `k_soil` are duplicated/derived from table values provided by the U.S. Department of Agriculture, Agricultural Research Service,  (USDA ARS, 2021) as follows:

- `theta_res` = `x`, where x is taken from column `qr` [-];
- `theta_sat` = `x`, where x is taken from column `qs` [-];
- `vGn_alpha` = `-1 * 10^x * 100`, where x is taken from column `log(a)` [log(cm-1)] > [m-1] while also changing the sign of the alpha parameter so that matric head calculations use SUMMA's convention of negative values for matric head;
- `vGn_n` = `10^x`, where x is taken from column `log(n)` [-];
- `k_soil` = `10^x / 100 / (60 * 60 * 24)`, where x is taken from column `Ks` [log(cm d-1)] > [m s-1].


## References
NCAR Research Applications Laboratory, RAL, 2021. Noah-Multiparameterization Land Surface Model (Noah-MP LSM) [WWW Document]. URL https://ral.ucar.edu/solutions/products/noah-multiparameterization-land-surface-model-noah-mp-lsm (accessed 11.15.21).

U.S. Department of Agriculture: Agricultural Research Service (USDA ARS), 2021. ROSETTA Class Average Hydraulic Parameters [WWW Document]. URL https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-research-unit/docs/model/rosetta-class-average-hydraulic-parameters/ (accessed 11.15.21).
