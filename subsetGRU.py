# Subset out a NA HRU forcing, parameter, and attribute files where GRU matches HRU
#
import os
os.system('ls')

nbasin_slurm = 1 # first basin in the slurm batch run
nscript_slurm = 20 # 20th script
GRU_id =154 + nbasin_slurm + 207*nscript_slurm
HRU_id =154 + nbasin_slurm + 207*nscript_slurm

settingsPath = '/globalhome/gwu479/HPC/TestScripts/settings/'
forcingPath  = '/project/gwf/gwf_cmt/rez299/summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input/'
attributeFile_in = settingsPath + 'attributes.nc' # Relative to settingsPath
trialParamFile_in = settingsPath + 'trialParams.nc' # Relative to settingsPath
forcingFile_in = forcingPath + 'forcing.nc' # Relative to forcingPath WHAT IS THIS

desPath = '/globalhome/gwu479/HPC/TestScripts/subset_'+str(GRU_id)+'/' # Will need to change this in fileManger for settingsPath and forcingPath
attributeFile_out = desPath + 'attributes.nc' # Relative to settingsPath
trialParamFile_out = desPath + 'trialParams.nc' # Relative to settingsPath
forcingFile_out = desPath + 'forcing.nc' # Relative to forcingPath WHAT IS THIS

os.system('ncks -d gru GRU_id,GRU_id -d hru HRU_id,HRU_id forcingFile_in forcingFile_out')
os.system('ncks -d gru GRU_id,GRU_id -d hru HRU_id,HRU_id attributeFile_in attributeFile_out')
os.system('ncks -d gru GRU_id,GRU_id -d hru HRU_id,HRU_id trialParamFile_in trialParamFile_out')
