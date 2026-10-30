#!/bin/bash
#SBATCH --job-name=atha_fdds_chk
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --output=athabasca_fortran_dds_check_%j.out
# One model run with the 15 DEFAULT parameters through the Fortran build and the Fortran-DDS config files.
# Must print KGE 0.5528 (= J() with the defaults), which confirms that the config files, the
# seg_outlet change and the build reproduce the J() setup.
# Usage (from summa-spotpy/): sbatch my_work_athabasca/scripts/fortran_dds/run_fortran_dds_check.slurm <tag>
TAG=${1:?run tag, e.g. smoke}
BASIN=<WORK_ROOT>/summa-spotpy/my_work_athabasca
DEV=<WORK_ROOT>/summa-develop
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3
EXE=$(ls $DEV/bin/summa*_mizuroute.exe | grep -v _opt | grep -v _mpi | head -1)
RUN=$BASIN/fortran_dds/$TAG
echo "Check run: $EXE, $RUN, started $(date)"
PARAMS=$(tail -n +2 $BASIN/scripts/athabasca_params.csv | awk -F, '{printf " --param %s %s", $1, $2}')
echo "parameters:$PARAMS"
cd $RUN || exit 1
mpirun -np 1 $EXE -m $RUN/fileManager.txt -c $RUN/config.toml $PARAMS > check_run.log 2>&1
EXIT=$?
echo "exit code $EXIT"
grep -i "subset mode\|error\|stop" check_run.log | head -5
echo "objective (J() with the defaults gives 5.528168198395311E-01):"
grep -E "^ *-?[0-9]\.[0-9]{6,}E[+-][0-9]{2} *$" check_run.log | tail -1
echo "Finished $(date)"
