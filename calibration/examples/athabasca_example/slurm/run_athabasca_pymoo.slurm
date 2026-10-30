#!/bin/bash
#SBATCH --job-name=atha_pymoo
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=48:00:00
#SBATCH --output=athabasca_pymoo_%j.out
#SBATCH --error=athabasca_pymoo_%j.err
# Usage (from summa-spotpy/); all arguments go to calibrate_athabasca_pymoo.py:
#   smoke test : sbatch --time=02:00:00 my_work_athabasca/scripts/pymoo/run_athabasca_pymoo.slurm --evals 147 --tag smoke
#                (98 starting sets + one generation of 49 = 3 waves of 49 workers)
#   full run   : sbatch my_work_athabasca/scripts/pymoo/run_athabasca_pymoo.slurm
#                (5,000 budget -> 98 + 100 x 49 = 4,998 evaluations; results in results/pymoo/ga_<jobid>/)
#   resume     : sbatch my_work_athabasca/scripts/pymoo/run_athabasca_pymoo.slurm --tag <jobid of the first run> --resume
#   other seed : sbatch my_work_athabasca/scripts/pymoo/run_athabasca_pymoo.slurm --seed 43
echo "Job $SLURM_JOB_ID on $(hostname): Athabasca Pymoo GA (develop, noTS, period-end obs), $SLURM_NTASKS processes, args: $*, started $(date)"
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
REPO=<WORK_ROOT>/summa-spotpy
cd $REPO || exit 1
ls -l ../summa-develop/build/cmake_build/libsumma.so my_work_athabasca/settings/athabasca_config_develop_noTS.toml my_work_athabasca/scripts/eval_one_develop.py my_work_athabasca/scripts/athabasca_params.csv || exit 1
python -c "import pymoo, dill, mpi4py; print('pymoo', pymoo.__version__, 'mpi4py', mpi4py.__version__)" || exit 1
export OMPI_MCA_btl='^openib'
mpirun -np $SLURM_NTASKS python -u my_work_athabasca/scripts/pymoo/calibrate_athabasca_pymoo.py "$@"
EXIT_CODE=$?
echo "Finished $(date) with exit code $EXIT_CODE"
exit $EXIT_CODE
