#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=1-00:00
#SBATCH --job-name=CONCAT
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load gcc/9.3.0
module load geo-stack
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index h5netcdf
pip install --no-index h5py
pip install --no-index xarray
pip install --no-index netCDF4

python concat_groups_split_summa.py sundials_1en8
