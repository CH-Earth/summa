#!/bin/bash
#SBATCH --job-name=atha_fdds
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=96:00:00
#SBATCH --output=athabasca_fortran_dds_%j.out
#SBATCH --error=athabasca_fortran_dds_%j.err
# SUMMA develop's own async DDS (Fortran) on Athabasca, built with the seed + default-start patch
# (scripts/fortran_dds/patches/fortran_dds_seed_defstart.patch) as summa_mizuroute_opt_defstart.exe.
# 1 coordinator + 49 workers, the same layout as the Python async DDS runs.
#   SUMMA_DDS_SEED=<n>        random seed (default 42, as in the unpatched code)
#   SUMMA_DDS_START=default   start from the SUMMA default values (otherwise random, as unpatched)
# Usage (from summa-spotpy/), e.g. seed 42 from the defaults:
#   sbatch --partition=wholenode --time=3-12:00:00 --job-name=fdds_d42 \
#     --export=ALL,SUMMA_DDS_SEED=42,SUMMA_DDS_START=default \
#     my_work_athabasca/scripts/fortran_dds/run_fortran_dds_defstart.slurm defstart_s42
# The run folder my_work_athabasca/fortran_dds/<tag>/ must first be made with setup_fortran_dds.py
# (copy coldState.nc into its state folder).
TAG=${1:?run tag, e.g. defstart_s42}
BASIN=<WORK_ROOT>/summa-spotpy/my_work_athabasca
DEV=<WORK_ROOT>/summa-develop
RUN=$BASIN/fortran_dds/$TAG
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3
export OMP_NUM_THREADS=1 OMPI_MCA_btl='^openib'
EXE=$DEV/bin/summa_mizuroute_opt_defstart.exe   # built with scripts/fortran_dds/patches/fortran_dds_seed_defstart.patch
echo "Job $SLURM_JOB_ID on $(hostname): Fortran async DDS, $EXE, run folder $RUN, $SLURM_NTASKS processes, started $(date)"
ls -l $EXE $RUN/fileManager.txt $RUN/config.toml || exit 1
grep -E "n_samples|theta_sat" $RUN/config.toml $BASIN/settings/SUMMA/localParamInfo_fortran.txt
cd $RUN
mpirun -np $SLURM_NTASKS $EXE -m $RUN/fileManager.txt -c $RUN/config.toml > $RUN/dds_run.log 2>&1
EXIT_CODE=$?
echo "new-best lines: $(grep -c 'new DDS best' $RUN/dds_run.log); last: $(grep 'new DDS best' $RUN/dds_run.log | tail -1)"
echo "output files:"; ls -l $RUN/output | head
echo "Finished $(date) with exit code $EXIT_CODE"
exit $EXIT_CODE
