#!/bin/bash
#SBATCH --job-name=atha_xcheck
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mem=24GB
#SBATCH --time=02:00:00
#SBATCH --output=athabasca_crosscheck_%j.out
# Cross-check of the two spin-up methods with the best parameter sets (three single SUMMA runs at once):
#   A: async DDS v2 best  under the Fortran method (start 2015-09-01 from the default-parameter spinup restart)
#   B: Fortran DDS best   under the J() method     (cold start 2014-09-01, own parameters from the start)
#   C: Fortran DDS best   under the Fortran method (sanity check: must reproduce its calibration KGE)
# Usage (from summa-spotpy/): sbatch my_work_athabasca/scripts/fortran_dds/run_crosscheck.slurm
BASIN=<WORK_ROOT>/summa-spotpy/my_work_athabasca
DEV=<WORK_ROOT>/summa-develop
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3 conda/2024.09
source <WORK_ROOT>/venv-spotpy/bin/activate
EXE=$(ls $DEV/bin/summa*_mizuroute.exe | grep -v _opt | grep -v _mpi | head -1)
FULL=$BASIN/fortran_dds/full
CC=$BASIN/fortran_dds/crosscheck
V2=$BASIN/results/asyncdds_mpi_v2_21608093/summary.json
FD=$FULL/results/summary.json
RESTART=$(ls $FULL/state/*_restart_*_spinup_rank0000.nc | head -1)
echo "Cross-check on $(hostname), $EXE, spinup restart $RESTART, started $(date)"
rm -rf $CC; mkdir -p $CC

params() {   # --param name value ... from a summary.json
  python -c "import json,sys; p=json.load(open(sys.argv[1]))['best_params']; print(' '.join(f'--param {k} {v!r}' for k,v in p.items()))" "$1"
}
fortran_method() {   # <case>: fileManager as in the Fortran run, starting from its spinup restart
  mkdir -p $CC/$1/output $CC/$1/state
  ln -s $RESTART $CC/$1/state/
  sed -e "s|^statePath .*|statePath    '$CC/$1/state/'|" \
      -e "s|^outputPath .*|outputPath '$CC/$1/output/'|" \
      -e "s|^initConditionFile .*|initConditionFile    '$(basename $RESTART)'|" $FULL/fileManager.txt > $CC/$1/fileManager.txt
}
j_method() {         # <case>: fileManager as in the J() develop runs (cold start 2014-09-01)
  mkdir -p $CC/$1/output
  sed -e "s|^outputPath .*|outputPath '$CC/$1/output/'|" \
      -e "s|^outputControlFile .*|outputControlFile    'outputControl_develop.txt'|" $BASIN/settings/SUMMA/fileManager.txt > $CC/$1/fileManager.txt
}
fortran_method A; j_method B; fortran_method C
for c in A B C; do grep -n "simStartTime\|statePath\|initConditionFile\|globalHruParamFile" $CC/$c/fileManager.txt | sed "s|^|$c: |"; done

cd $CC
( cd A && mpirun -np 1 $EXE -m $CC/A/fileManager.txt -c $FULL/config.toml $(params $V2) > run.log 2>&1 ) &
( cd B && mpirun -np 1 $EXE -m $CC/B/fileManager.txt -c $FULL/config.toml $(params $FD) > run.log 2>&1 ) &
( cd C && mpirun -np 1 $EXE -m $CC/C/fileManager.txt -c $FULL/config.toml $(params $FD) > run.log 2>&1 ) &
wait

kge() { grep -E "^ *-?[0-9]\.[0-9]{6,}E[+-][0-9]{2} *$" $CC/$1/run.log | tail -1; }
echo
echo "A  v2 best      under Fortran method : $(kge A)    (under J(): 0.8925)"
echo "B  Fortran best under J() method     : $(kge B)    (under Fortran method: 0.8056)"
echo "C  Fortran best under Fortran method : $(kge C)    (sanity check, must be ~0.8056)"
for c in A B C; do grep -i -m1 "error" $CC/$c/run.log | sed "s|^|$c error: |"; done
echo "Finished $(date)"
