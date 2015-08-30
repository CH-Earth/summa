#!/opt/local/bin/bash
#
# used to run SUMMA for all the PLUMBER sites

# define the SUMMA directory
summaDir=/Users/mclark/summa

# Define the SUMMA executable
summaEXE=${summaDir}/bin/summa.exe

# Define the path for the log files
logPath=${summaDir}/output/plumber/log

# Define the file path for the SUMMA settings
settingsPath=${summaDir}/settings/plumber

# Define the manager path
managerPath=${settingsPath}/fileManager

# define the experiment name
expName=_initialPlumberTest

# run summa for each site
for managerFile in $( ls  ${managerPath}/* ); do

 # get the logfile name
 IFS='_' read -a strarr <<< "${managerFile}"
 logName=${logPath}/summaLog_${strarr[-1]}

 # run the model for a given site
 echo $summaEXE $expName $managerFile $logName
 #$summaEXE $expName $managerFile > ${logName}

done
