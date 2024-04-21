#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=6GB
#SBATCH --time=1-00:00
#SBATCH --job-name=STAT
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load StdEnv/2020
module load gcc/9.3.0
module load geo-stack/2022c
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

python timeseries_to_statistics.py sundials_1en6 1 1
python timeseries_to_statistics.py sundials_1en6 2 1
python timeseries_to_statistics.py be1 1 1
python timeseries_to_statistics.py be1 2 1
python timeseries_to_statistics.py be32 1 1
python timeseries_to_statistics.py be32 2 1
python timeseries_to_statistics.py be16 1 1
python timeseries_to_statistics.py be16 2 1
