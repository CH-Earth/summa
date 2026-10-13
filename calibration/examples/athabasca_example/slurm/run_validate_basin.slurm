#!/bin/bash
# Usage: sbatch run_validate_basin.slurm <ross|athabasca>
#SBATCH -J dev_check
#SBATCH --account=<ACCOUNT>
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=0-2
#SBATCH -o devcheck_%x_%A_%a.out
#SBATCH -e devcheck_%x_%A_%a.err
# One J() evaluation of the basin's best parameters with the old build and the develop port
# (output on and off); the three tasks run side by side.
BASIN=${1:?basin: ross or athabasca}
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
V=(old dev_ts dev_nots)
cd <WORK_ROOT>/summa-spotpy/my_work_$BASIN
python -u ../my_work_athabasca/scripts/develop_port/validate_basin.py $BASIN ${V[$SLURM_ARRAY_TASK_ID]} $SLURM_ARRAY_JOB_ID
