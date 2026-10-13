#!/bin/bash
#SBATCH --job-name=atha_asyncdds
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=96:00:00
#SBATCH --output=athabasca_asyncdds_%j.out
#SBATCH --error=athabasca_asyncdds_%j.err
EVALS=${1:-5000}
echo "Job $SLURM_JOB_ID on $(hostname): Athabasca async DDS, $EVALS evaluations, $SLURM_NTASKS processes, started $(date)"
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
REPO=<WORK_ROOT>/summa-spotpy
cd $REPO || exit 1
ls -l calibration/libsumma.dylib build/cmake_build/libsumma.so || exit 1
export OMPI_MCA_btl='^openib'
grep -n "N_PARAMS = " build/source/driver/summa_c_api.f90
mpirun -np $SLURM_NTASKS python -u my_work_athabasca/scripts/calibrate_athabasca_asyncdds_mpi.py $EVALS 0.2 $SLURM_JOB_ID
EXIT_CODE=$?
echo "Finished $(date) with exit code $EXIT_CODE"
exit $EXIT_CODE
