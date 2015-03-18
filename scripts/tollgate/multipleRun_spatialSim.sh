#!/bin/bash
#
# used to run SUMMA multiple times
#
#
# define the number of processors
nProcessors=24
#
# define the executable
exRun='/home/mclark/summa/bin/summa.exe '
#
# define the path for the settings file
settingsPath='/home/mclark/summa/settings/tollgate/fileManager/'
#
# define the path for the log files
logPath='/home/mclark/summa/logFilez/tollgate/'

#
# loop through basins
for fileManager in $( ls ${settingsPath} )
do

 # print progress
 echo $fileManager

 # split the string using underscores
 IFS='_' read -a strarr <<< "${fileManager}"

 # split last string element using periods (to remove ".txt")
 IFS='.' read -a suffix <<< "${strarr[-1]}"

 # define the basin
 iBasin=${suffix[0]}

 # define the log file
 exLog=${logPath}ex${iBasin}.log

 # define the name of the experiment
 exName=exBasin${iBasin}

 # run the experiment
 ./runExperiment.sh $exRun $exName ${settingsPath}${fileManager} $exLog &
 sleep 1

 # sleep if using the desired number of processors
 while [ `ls -1 summaCopy.*.control | wc -l` -ge $nProcessors ]
 do
  sleep 30
 done

done

exit 0
