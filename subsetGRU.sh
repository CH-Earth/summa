#!/bin/bash
# Subset out a NA HRU forcing, parameter, and attribute files where GRU matches HRU
#
# Inside forcingFile_out need to change forcingPath to = desForcingPath
# Inside forcingFile_out need to change initConditionFile, attributeFile, trialParamFile to = *out versions

module load StdEnv/2020
module load gcc/9.3.0
module load nco/5.0.6

# GRU want to subset, change to do another GRU
# first basin in the slurm batch run, 20th script
nbasin_slurm = 1
nscript_slurm = 20
GRU_id =$((154 + nbasin_slurm + 207*nscript_slurm))
HRU_id =$((154 + nbasin_slurm + 207*nscript_slurm))

# top paths, could change
homePath = /globalhome/gwu479/HPC/TestScripts/
summa_exe= /globalhome//gwu479/HPC/SummaSundials/summa/bin/summa_sundials.exe
forcingPath  = /project/gwf/gwf_cmt/rez299/summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input/

# in paths, probably won't change
fileManager_in  = ${homePath}fileManager.txt
settingsPath = ${homePath}settings/
initConditionFile_in = ${settingsPath}coldState.nc
attributeFile_in = ${settingsPath}attributes.nc
trialParamFile_in = ${settingsPath}trialParams.nc

# out paths, probably won't change
fileManager_out  = ${homePath}fileManager_${GRU_id}.txt
initConditionFile_out = ${settingsPath}coldState_${GRU_id}.nc
attributeFile_out = ${settingsPath}attributes_${GRU_id}.nc
trialParamFile_out = ${settingsPath}trialParams_${GRU_id}.nc
desForcingPath = ${homePath}forcing_${GRU_id})/

# set up directory and new file Manager (will have to change things in it manually as above)
mkdir $desForcingPath
cp $fileManager_in $fileManager_out

# do the subset
ncks -d gru $GRU_id,$GRU_id -d hru $HRU_id,$HRU_id $initConditionFile_in $initConditionFile_out
ncks -d gru $GRU_id,$GRU_id -d hru $HRU_id,$HRU_id $attributeFile_in $attributeFile_out
ncks -d gru $GRU_id,$GRU_id -d hru $HRU_id,$HRU_id $trialParamFile_in $trialParamFile_out

# forcing subset has multiple files
for fn in ${forcingPath}*.nc; do
    output_fn = ${desForcingPath}${fn}
    output = ncks -d gru $GRU_id,$GRU_id -d hru $HRU_id,$HRU_id $fn $output_fn

# write summa command call file
runFile = ${homePath}run_${GRU_id}.sh
summaCommand = $summa_exe -p never -s _testSumma -m $fileManger_out -r e
$summaCommand > runFile
