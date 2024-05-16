#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=0-05:00
#SBATCH --job-name=STATB
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out
echo "$SLURM_ARRAY_TASK_ID"

# ----------------------------------------------------------------------------------------------
# RUN WITH:
# sbatch --array1-[number of jobs] [script name]
# sbatch --array=1-200 process_statBatch.sh
# sbatch --array=201-201 process_statBatch.sh
# ----------------------------------------------------------------------------------------------

module load StdEnv/2020
module load gcc/9.3.0
module load geo-stack/2022c
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

python timeseries_to_statistics.py sundials_1en6 $SLURM_ARRAY_TASK_ID 200
python timeseries_to_statistics.py be1 $SLURM_ARRAY_TASK_ID 200
python timeseries_to_statistics.py be32 $SLURM_ARRAY_TASK_ID 200
python timeseries_to_statistics.py be16 $SLURM_ARRAY_TASK_ID 200
