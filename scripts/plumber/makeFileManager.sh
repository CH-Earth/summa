#!/opt/local/bin/bash
#
# used to make a list of forcing files for each GRU
# (based on the localAttributes files)

# define file path for the vegetation data
dataPath=/Users/mclark/summa/input/plumber/

# define the path to the file manager
managerPath=/Users/mclark/summa/settings/plumber/fileManager/

# define pattern to replace with actual data file
cPattern='XX_siteName_XX'

# loop through files in the local attributes folder
for dataFile in $( ls  ${dataPath} ); do

 # split the string using underscores
 IFS='_' read -a strTemp0 <<< "${dataFile}"

 # split the string using periods
 IFS='.' read -a strTemp1 <<< "${strTemp0[-1]}"

 # define file suffix
 siteName=${strTemp1[0]}

 # define the manager file
 managerFile=${managerPath}summa_fileManager_${siteName}.txt

 # make the file manager file
 cp templates/summa_fileManager_template.txt ${managerFile}

 # replace the search string with the file suffix
 gsed -i 's/'${cPattern}'/'${siteName}'/g' ${managerFile}

 # print progress
 echo $managerFile
 
done
