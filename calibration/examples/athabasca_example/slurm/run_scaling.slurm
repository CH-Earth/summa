#!/bin/bash
# Strong-scaling run of the Athabasca async DDS calibration on whole Anvil nodes.
# Submitted by submit_scaling.sh, which creates the run folder and sets nodes, tasks, time and log paths.
# Arguments: $1 = number of evaluations, $2 = run tag (= folder name under scaling_tests/)
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=wholenode
#SBATCH --cpus-per-task=1
EVALS=${1:?number of evaluations}
TAG=${2:?run tag}
BASIN=<WORK_ROOT>/summa-spotpy/my_work_athabasca
export ASYNC_OUTDIR=$BASIN/scaling_tests/$TAG                       # everything for this run ends up here
export ASYNC_WORKROOT=<SCRATCH>/athabasca_scaling          # per-rank SUMMA output (large, temporary)
mkdir -p $ASYNC_OUTDIR $ASYNC_WORKROOT
{
echo "tag            $TAG"
echo "job_id         $SLURM_JOB_ID"
echo "nodes          $SLURM_JOB_NUM_NODES"
echo "node_list      $SLURM_JOB_NODELIST"
echo "processes      $SLURM_NTASKS"
echo "workers        $((SLURM_NTASKS-1))"
echo "evaluations    $EVALS"
echo "partition      $SLURM_JOB_PARTITION"
echo "time_limit     $(squeue -h -j $SLURM_JOB_ID -o %l)"
echo "start          $(date '+%F %T %Z')"
echo "work_dir       $ASYNC_WORKROOT/asyncdds_$TAG"
echo "setup          develop build, athabasca_config_develop_noTS.toml (period-end obs, KW); v2: only successful SUMMA runs count"
} > $ASYNC_OUTDIR/run_info.txt
cat $ASYNC_OUTDIR/run_info.txt
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
export OMP_NUM_THREADS=1
export OMPI_MCA_btl='^openib'
cd $BASIN/.. || exit 1
T0=$(date +%s)
mpirun -np $SLURM_NTASKS python -u my_work_athabasca/scripts/calibrate_athabasca_asyncdds_mpi_scaling_develop_v2.py $EVALS 0.2 $TAG
EXIT_CODE=$?
WALL=$(( $(date +%s) - T0 ))
{
echo "end            $(date '+%F %T %Z')"
echo "walltime_s     $WALL"
echo "walltime_h     $(awk -v s=$WALL 'BEGIN{printf "%.3f", s/3600}')"
echo "exit_code      $EXIT_CODE"
} >> $ASYNC_OUTDIR/run_info.txt
echo "WALLTIME_S $WALL  nodes $SLURM_JOB_NUM_NODES  workers $((SLURM_NTASKS-1))  evals $EVALS  tag $TAG"
exit $EXIT_CODE
