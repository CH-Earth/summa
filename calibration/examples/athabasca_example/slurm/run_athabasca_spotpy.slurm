#!/bin/bash
#SBATCH --job-name=atha_spotpy
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=96:00:00
#SBATCH --output=athabasca_spotpy_%j.out
#SBATCH --error=athabasca_spotpy_%j.err
# Usage (from summa-spotpy/):
#   SCE-UA smoke : sbatch --time=02:00:00 my_work_athabasca/scripts/spotpy/run_athabasca_spotpy.slurm sceua 62 2 1
#   SCE-UA full  : sbatch my_work_athabasca/scripts/spotpy/run_athabasca_spotpy.slurm sceua 5000
#   DDS smoke    : sbatch --time=02:00:00 my_work_athabasca/scripts/spotpy/run_athabasca_spotpy.slurm dds 99 0 1
#   DDS full     : sbatch --time=48:00:00 my_work_athabasca/scripts/spotpy/run_athabasca_spotpy.slurm dds 5000
#   LHS full     : sbatch my_work_athabasca/scripts/spotpy/run_athabasca_spotpy.slurm lhs 5000
# (4th argument 1 = smoke test: results go to <algo>_smoke_<jobid>; for SCE-UA it also stops after burn-in)
ALGO=${1:-sceua}; REPS=${2:-5000}; NGS=${3:-0}; MAXLOOP=${4:-0}
TAG=$SLURM_JOB_ID; [ "$MAXLOOP" = "1" ] && TAG=smoke_$SLURM_JOB_ID
echo "Job $SLURM_JOB_ID on $(hostname): Athabasca SPOTPY $ALGO (develop, noTS, period-end obs), $REPS repetitions, ngs=$NGS (0 = one per worker), max_loop_inc=$MAXLOOP, $SLURM_NTASKS processes, started $(date)"
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
REPO=<WORK_ROOT>/summa-spotpy
cd $REPO || exit 1
ls -l ../summa-develop/build/cmake_build/libsumma.so my_work_athabasca/settings/athabasca_config_develop_noTS.toml my_work_athabasca/scripts/eval_one_develop.py || exit 1
python -c "import spotpy, mpi4py; print('spotpy', spotpy.__version__, 'mpi4py', mpi4py.__version__)" || exit 1
export OMPI_MCA_btl='^openib'
mpirun -np $SLURM_NTASKS python -u my_work_athabasca/scripts/spotpy/calibrate_athabasca_spotpy.py $ALGO $REPS $NGS $MAXLOOP $TAG
EXIT_CODE=$?
echo "Results: my_work_athabasca/results/spotpy/${ALGO}_${TAG}/"
echo "Finished $(date) with exit code $EXIT_CODE"
exit $EXIT_CODE
