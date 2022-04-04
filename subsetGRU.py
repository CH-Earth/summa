# Subset out a NA HRU forcing, parameter, and attribute files where GRU matches HRU
#
# Inside forcingFile_out need to change forcingPath to = desForcingPath
# Inside forcingFile_out need to change initConditionFile, attributeFile, trialParamFile to = *out versions

import glob
import subprocess

# GRU want to subset, change to do another GRU
nbasin_slurm = 1 # first basin in the slurm batch run
nscript_slurm = 20 # 20th script
GRU_id =154 + nbasin_slurm + 207*nscript_slurm
GRU_bounds = str(GRU_id) + ',' + str(GRU_id)
HRU_id =154 + nbasin_slurm + 207*nscript_slurm
HRU_bounds = str(HRU_id) + ',' + str(HRU_id)

# top paths, could change
homePath = '/globalhome/gwu479/HPC/TestScripts/'
summa_exe= '/globalhome//gwu479/HPC/SummaSundials/summa/bin/summa_sundials.exe'
forcingPath  = '/project/gwf/gwf_cmt/rez299/summaWorkflow_data/domain_NorthAmerica/forcing/4_SUMMA_input/'

# in paths, probably won't change
fileManager_in  = homePath + 'fileManager.txt'
settingsPath = homePath + 'settings/'
initConditionFile_in = settingsPath + 'coldState.nc' # Relative to settingsPath
attributeFile_in = settingsPath + 'attributes.nc' # Relative to settingsPath
trialParamFile_in = settingsPath + 'trialParams.nc' # Relative to settingsPath

# out paths, probably won't change
fileManager_out  = homePath + 'fileManager_' + str(GRU_id) + '.txt'
initConditionFile_out = settingsPath + 'coldState_' + str(GRU_id) + '.nc' # Relative to settingsPath
attributeFile_out = settingsPath + 'attributes_' + str(GRU_id) + '.nc' # Relative to settingsPath
trialParamFile_out = settingsPath + 'trialParams_' + str(GRU_id) + '.nc' # Relative to settingsPath
desForcingPath = homePath + 'forcing_' + str(GRU_id)+'/'

# set up directory and new file Manager (will have to change things in it manually as above)
subprocess.call(['mkdir', desForcingPath])
subprocess.call(['cp', fileManager_in, fileManager_out])

# do the subset
subprocess.call(['ncks', '-d', 'gru', GRU_bounds, '-d', 'hru', HRU_bounds, initConditionFile_in initConditionFile_out])
subprocess.call(['ncks', '-d', 'gru', GRU_bounds, '-d', 'hru', HRU_bounds, attributeFile_in attributeFile_out])
subprocess.call(['ncks', '-d', 'gru', GRU_bounds, '-d', 'hru', HRU_bounds, trialParamFile_in trialParamFile_out])

# forcing subset has multiple files
for fn in glob.glob(forcingPath +'*.nc'):
    output_fn = desForcingPath + fn
    output = subprocess.call(['ncks', '-d', 'gru', GRU_bounds, '-d', 'hru', HRU_bounds, fn, output_fn])

# write summa command call file
runFile = homePath + 'run_' + str(GRU_id) + '.sh'
summaCommand = summa_exe + '-p never -s _testSumma -m' + fileManger_out + '-r e'
runFile = open(homePath + 'run_' + str(GRU_id) + '.sh', "w")
a = runFile.write(summaCommand)
runFile.close()
