#!/bin/bash
#
# used to spinup SUMMA fora given PLUMBER site

# define the site name
siteName=$1
#echo $siteName

# define the SUMMA directory
summaDir=/home/mclark/summa

# Define the SUMMA executable
summaEXE=${summaDir}/bin/summa.exe

# Define the path for the log files
logPath=${summaDir}/output/plumber/log

# Define the path for the output files
outputPath=${summaDir}/output/plumber/orig

# Define the file path for the SUMMA settings
settingsPath=${summaDir}/settings/plumber

# Define the manager path
managerPath=${settingsPath}/spinupFileManager

# define the manager file
managerFile=${managerPath}/summa_fileManager_${siteName}.txt

# get the initial conditions
cp ${settingsPath}/summa_zInitialCond.txt ${settingsPath}/initialCond/summa_zInitialCond_${siteName}.txt

# loop through trials
for iTrial in `seq -w 1 10`; do

 # define the experiment name
 expName=_spinupPlumber${iTrial}
 #echo $siteName $expName

 # define the log name
 logName=${logPath}/summaLog${expName}_${siteName}.txt

 # run the model for a given site
 $summaEXE $expName $managerFile > ${logName}

 # define the initial conditions file
 initCondFile=`ls -1rt ${outputPath}/${siteName}*${expName}.txt | tail -n1`

 # copy the initial conditions file to the initial conditions directory
 cp ${initCondFile} ${settingsPath}/initialCond/summa_zInitialCond_${siteName}.txt

done
