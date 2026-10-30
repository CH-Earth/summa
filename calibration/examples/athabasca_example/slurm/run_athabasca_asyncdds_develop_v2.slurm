#!/bin/bash
#SBATCH --job-name=atha_dds_v2
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=96:00:00
#SBATCH --output=athabasca_dds_v2_%j.out
#SBATCH --error=athabasca_dds_v2_%j.err
# Async DDS v2 (log transforms + ordered constraints from athabasca_dds_v2_settings.toml).
# Usage (from summa-spotpy/); arguments go to the driver:
#   settings check: python my_work_athabasca/scripts/calibrate_athabasca_asyncdds_mpi_develop_v2.py --check
#   smoke test    : sbatch --time=02:00:00 my_work_athabasca/scripts/run_athabasca_asyncdds_develop_v2.slurm 147 0.2 smoke
#   full run      : sbatch my_work_athabasca/scripts/run_athabasca_asyncdds_develop_v2.slurm
#   whole node    : sbatch --partition=wholenode my_work_athabasca/scripts/run_athabasca_asyncdds_develop_v2.slurm
#   other options : ... run_athabasca_asyncdds_develop_v2.slurm 5000 0.2 <tag> --bounds reject --start random
echo "Job $SLURM_JOB_ID on $(hostname): Athabasca async DDS v2 (develop, noTS, period-end obs), $SLURM_NTASKS processes, args: $*, started $(date)"
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
REPO=<WORK_ROOT>/summa-spotpy
cd $REPO || exit 1
ls -l ../summa-develop/build/cmake_build/libsumma.so my_work_athabasca/settings/athabasca_config_develop_noTS.toml \
      my_work_athabasca/scripts/eval_one_develop.py my_work_athabasca/scripts/athabasca_dds_v2_settings.toml || exit 1
mpirun -np 1 python -u my_work_athabasca/scripts/calibrate_athabasca_asyncdds_mpi_develop_v2.py --check "$@" || exit 1
export OMPI_MCA_btl='^openib'
mpirun -np $SLURM_NTASKS python -u my_work_athabasca/scripts/calibrate_athabasca_asyncdds_mpi_develop_v2.py "$@"
EXIT_CODE=$?
echo "Finished $(date) with exit code $EXIT_CODE"
exit $EXIT_CODE
