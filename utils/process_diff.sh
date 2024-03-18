#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00
#SBATCH --job-name=DIFF
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load gcc/9.3.0
module load geo-stack
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

python check_bit_4_bit_withTol.py /home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/summa-sundials_be_noJac/ /home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/summa-be/ 0.1 > out.txt
