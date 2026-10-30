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
# SUMMA develop's own async DDS (Fortran, summa_mizuroute_opt.exe) on Athabasca:
# 1 coordinator + 49 workers, the same layout as the Python async DDS runs.
# Usage (from summa-spotpy/):
#   smoke test : sbatch --time=02:00:00 my_work_athabasca/scripts/fortran_dds/run_fortran_dds.slurm smoke
#   full run   : sbatch my_work_athabasca/scripts/fortran_dds/run_fortran_dds.slurm full
# The run folder my_work_athabasca/fortran_dds/<tag>/ must first be made with setup_fortran_dds.py.
TAG=${1:?run tag (smoke or full)}
BASIN=<WORK_ROOT>/summa-spotpy/my_work_athabasca
DEV=<WORK_ROOT>/summa-develop
RUN=$BASIN/fortran_dds/$TAG
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3
export OMP_NUM_THREADS=1 OMPI_MCA_btl='^openib'
EXE=$(ls $DEV/bin/summa*_mizuroute_opt.exe | head -1)
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
