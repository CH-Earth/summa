#!/bin/bash
#SBATCH -J atha_hydro_dev
#SBATCH --account=<ACCOUNT>
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=0-2
#SBATCH -o athabasca_hydro_dev_%A_%a.out
#SBATCH -e athabasca_hydro_dev_%A_%a.err
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
export OMPI_MCA_btl='^openib'
POINTS=(develop_best spotpy_best symfluence_best)
cd <WORK_ROOT>/summa-spotpy/my_work_athabasca
python -u scripts/run_best_hydrograph_develop.py ${POINTS[$SLURM_ARRAY_TASK_ID]} ${SLURM_ARRAY_JOB_ID}
