#!/bin/bash
# Subset out a NA HRU forcing, parameter, and attribute files where GRU matches HRU
#
# Inside forcingFile_out need to change forcingPath to = desForcingPath
# Inside forcingFile_out need to change initConditionFile, attributeFile, trialParamFile to = _${GRU_file} versions
# written originally by A. Van Beusekom 2023, hardwired paths run on Copernicus

module load StdEnv/2020
module load gcc/9.3.0
module load nco/5.0.6

# GRU want to subset, change to do another GRU
# 155th basin in the slurm batch run, offset-th script (offset is SLURM_ARRAY_TASK_ID
# could also just set this as GRU_file from G* on output nc files (here the equation output is 487981).

offset=942
gruCount=518

# If know basin number in batch, add here, or comment out and use next if know GRU_file from the failure message id
#nbasin_slurm=25
#GRU_file=$((nbasin_slurm + gruCount*offset))

GRU_file=487980
nbasin_slurm=$((GRU_file-gruCount*offset))

GRU_id=$GRU_file
HRU_id=$GRU_id
echo "file name id is ${GRU_file}, GRU id is ${GRU_id}, HRU id is ${HRU_id}, number in batch is ${nbasin_slurm}"

# top paths, change these to yours
homeDir=/globalhome/gwu479/HPC/
desforceDir=/project/gwf/gwf_cmt/gwu479/
inpforceDir=/project/gwf/gwf_cmt/rez299/
homePath=${homeDir}TestScripts/
summa_exe=${homeDir}SummaSundials/summa/bin/summa_sundials.exe
forcingPath=${inpforceDir}summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input/

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
desForcingPath=${desforceDir}summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input_${GRU_file}/

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
for fn in NorthAmerica_remapped_*00-00-00.nc; do
    output_fn=${desForcingPath}${fn}
    ncks -d hru,$HRU_id,$HRU_id $fn $output_fn
    echo "${fn} HRU ${HRU_id} subsetted"
done
cd $homePath

# write summa command call file
runFile=${homePath}run_${GRU_file}.sh
echo "${summa_exe} -p never -s _testSumma -m ${fileManager_out} -r e" > $runFile
