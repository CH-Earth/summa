#!/bin/bash
#
# used to spinup SUMMA for all the PLUMBER sites

# define the SUMMA directory
summaDir=/home/mclark/summa

# Define the file path for the SUMMA settings
settingsPath=${summaDir}/settings/plumber

# Define the manager path
managerPath=${settingsPath}/fileManager

# run summa for each site
for managerFile in $( ls  ${managerPath}/* ); do

 # get the site name
 IFS='_' read -a strarr0 <<< "${managerFile}"
 IFS='.' read -a strarr1 <<< "${strarr0[-1]}"
 siteName=${strarr1[0]}

 # spinup the model for a given site
 ./spinupPlumber.sh $siteName &
 echo spinning up SUMMA for $siteName

done
