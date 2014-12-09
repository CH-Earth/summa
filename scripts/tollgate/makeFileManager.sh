#!/bin/bash
#
# used to create multiple duplicates of the fileManager files
#
# define file path for the vegetation data
vegPath=/home/mclark/summa/settings/tollgate/localAttributes/

# define file path for the forcing file list
managerPath=/home/mclark/summa/settings/tollgate/fileManager/

# define pattern to replace with actual data file
cPattern='XXXX'

# loop through files in the local attributes folder
for vegFile in $( ls  ${vegPath} ); do

 # split the string using underscores
 IFS='_' read -a strarr <<< "${vegFile}"

 # split last string element using periods (to remove ".txt")
 IFS='.' read -a suffix <<< "${strarr[-1]}"

 # define the basin
 iBasin=${suffix[0]}

 # define the manager file
 managerFile=snow_fileManager_basin_$iBasin.txt

 # make the filemanager file
 cp snow_fileManager__template.txt $managerPath$managerFile
 
 # replace the search string with the basin ID
 sed -i 's/'${cPattern}'/'${iBasin}'/g' $managerPath$managerFile
 
 # print progress
 echo filenm: $managerFile - $iBasin

done
