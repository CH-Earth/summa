#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50GB
#SBATCH --time=0-04:00
#SBATCH --job-name=GDB
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load StdEnv/2020
module load gcc/9.3.0
module load geo-stack/2022c
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

python plot_per_GRU.py sundials_1en6 rmse
python plot_per_GRU.py sundials_1en6 maxe
python plot_per_GRU.py sundials_1en6 kgem
