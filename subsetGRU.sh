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
# could also just set this as id from G* on output nc files (here 4295)
# However, this is not the actual GRUid, it is one off from the nc files, so subtract one before ncks
nbasin_slurm=1
nscript_slurm=20
GRU_file=$((154 + nbasin_slurm + 207*nscript_slurm))
GRU_id=$((GRU_file - 1))
HRU_id=$GRU_id
echo "file name id is 

# top paths, could change
homePath=/globalhome/gwu479/HPC/TestScripts/
summa_exe=/globalhome/gwu479/HPC/SummaSundials/summa/bin/summa_sundials.exe
forcingPath=/project/gwf/gwf_cmt/rez299/summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input/

# in paths, probably won't change
fileManager_in=${homePath}fileManager.txt
settingsPath=${homePath}settings/
initConditionFile_in=${settingsPath}coldState.nc
attributeFile_in=${settingsPath}attributes.nc
trialParamFile_in=${settingsPath}trialParams.nc

# out paths, probably won't change
fileManager_out=${homePath}fileManager_${GRU_file}.txt
initConditionFile_out=${settingsPath}coldState_${GRU_file}.nc
attributeFile_out=${settingsPath}attributes_${GRU_file}.nc
trialParamFile_out=${settingsPath}trialParams_${GRU_file}.nc
desForcingPath=/project/gwf/gwf_cmt/gwu479/summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input_${GRU_file}/

# set up directory and new file Manager (will have to change things in it manually as above)
mkdir -p "$desForcingPath"
cp $fileManager_in $fileManager_out

# do the subset
ncks -d hru,$HRU_id,$HRU_id $initConditionFile_in $initConditionFile_out
echo "coldState.nc HRU ${HRU_id} subsetted"
ncks -d gru,$GRU_id,$GRU_id -d hru,$HRU_id,$HRU_id $attributeFile_in $attributeFile_out
echo "attributes.nc GRU ${GRU_id} HRU ${HRU_id} subsetted"
ncks -d hru,$HRU_id,$HRU_id $trialParamFile_in $trialParamFile_out
echo "trialParams.nc HRU ${HRU_id} subsetted"

# forcing subset has multiple files
cd $forcingPath
for fn in *.nc; do
    output_fn=${desForcingPath}${fn}
    ncks -d hru,$HRU_id,$HRU_id $fn $output_fn
    echo "${fn} HRU ${HRU_id} subsetted"
done
cd $homePath

# write summa command call file
runFile=${homePath}run_${GRU_file}.sh
echo "${summa_exe} -p never -s _testSumma -m ${fileManager_out} -r e" > $runFile
