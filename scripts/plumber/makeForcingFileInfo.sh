#!/bin/bash
#
# used to make a list of forcing files for each GRU
# (based on the localAttributes files)

# define file path for the vegetation data
vegPath=/home/mclark/summa/settings/plumber/attributes/

# define file path for the forcing file info
infoPath=/home/mclark/summa/settings/plumber/forcingInfo/

# define file path for the forcing file list
listPath=/home/mclark/summa/settings/plumber/forcingList/

# define pattern to replace with actual data file
cPatternName='XX_forcingFileName_XX'
cPatternInfo='XX_forcingFileInfo_XX'

# loop through files in the local attributes folder
for vegFile in $( ls  ${vegPath} ); do

 # split the string using underscores
 IFS='_' read -a strarr <<< "${vegFile}"

 # define list filename (add the last element from the files in ${file_path})
 listFilenm=summa_zForcingFileList_${strarr[-1]}

 # define info filename (add the last element from the files in ${file_path})
 infoFilenm=summa_zForcingFileInfo_${strarr[-1]}

 # define data filename 
 dataFilenm=inputData_${strarr[-1]}

 # make the forcing list file
 cp templates/summa_zForcingList_template.txt $listPath$listFilenm

 # make the forcing info file
 cp templates/summa_zForcingInfo_template.txt $infoPath$infoFilenm

 # replace the search string with the info name 
 sed -i 's/'${cPatternInfo}'/'${infoFilenm}'/g' $listPath$listFilenm

 # replace the search string with the actual filename
 sed -i 's/'${cPatternName}'/'${dataFilenm}'/g' $infoPath$infoFilenm

 # print results
 echo filenm: $infoFilenm - $dataFilenm
 
done
