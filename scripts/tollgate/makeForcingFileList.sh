#!/bin/bash
#
# used to make a list of forcing files for each GRU
# (based on the localAttributes files)

# define file path for the vegetation data
vegPath=/home/mclark/summa/settings/tollgate/localAttributes/

# define file path for the forcing file list
listPath=/home/mclark/summa/settings/tollgate/forcingFileList/

# define the file path for the forcing info files
infoPath=tollgate/forcingFileInfo/

for vegFile in $( ls  ${vegPath} ); do

 # split the string using underscores
 IFS='_' read -a strarr <<< "${vegFile}"

 # define filename (add the last element from the files in ${file_path})
 listFilenm=snow_zForcingFileList_basin_${strarr[-1]}
 
 # split last string element using periods (to remove ".txt")
 IFS='.' read -a suffix <<< "${strarr[-1]}"

 # define the basin
 iBasin=${suffix[0]}

 # define the forcing info file
 infoFilenm=snow_zForcingFileInfo_basin_$iBasin.txt

 # loop through the input filename
 while read iLine ; do
  
  # split the string using spaces
  IFS=' ' read -a strarr <<< "${iLine}"

  # define the file for a given HRU (if the first element of the line is a number)
  if [[ ${strarr[0]} =~ ^[0-9]+$ ]]; then
   echo ${strarr[0]} "'"${infoPath}${infoFilenm}"'" > tempFile.$iBasin.${strarr[0]}.txt
  fi

 done < ${vegPath}${vegFile}

 # concatenate temporary files
 cat tempFile.$iBasin.*.txt > $listPath$listFilenm

 # remove temporary files
 rm tempFile.$iBasin.*.txt

 # print results
 echo filenm: $listFilenm - $iBasin - $infoFilenm
 
done
