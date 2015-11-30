#!/bin/bash
#
# used to make a list of forcing files for each GRU
# (based on the localAttributes files)

# define the decisions template
template=summa_fileManager_template.txt

# define file path for the vegetation data
dataPath=/home/mclark/summa/input/plumber

# define the path to the file manager
managerPath=/home/mclark/summa/settings/plumber/fileManager

# define pattern to replace with experiment name
cPatternExp='XX_expName_XX'

# define pattern to replace with site name
cPatternSite='XX_siteName_XX'

# loop through experiments
for ix in $( seq -f "%03g" 01 09 ); do

 # define experiment name
 expName=exp.03.${ix}
 echo $expName

 # loop through files in the local attributes folder
 for dataFile in $( ls  ${dataPath}/* ); do

  # remove the directory
  IFS='/' read -a strarr <<< "${dataFile}"
  fileName=${strarr[-1]}

  # split the string using underscores
  IFS='_' read -a strTemp0 <<< "${fileName}"

  # split the string using periods
  IFS='.' read -a strTemp1 <<< "${strTemp0[-1]}"

  # define file suffix
  siteName=${strTemp1[0]}

  # define the manager file
  managerFile=${managerPath}/${expName}/summa_fileManager_${siteName}.txt

  # make the file manager file
  cp templates/${template} ${managerFile}

  # replace the search string with the experiment name
  sed -i 's/'${cPatternExp}'/'${expName}'/g' ${managerFile}

  # replace the search string with the site name
  sed -i 's/'${cPatternSite}'/'${siteName}'/g' ${managerFile}

  # print progress
  echo $managerFile

 done
done
