#!/bin/bash
#
# used to make a list of forcing files for each GRU
# (based on the localAttributes files)

# define file path for the input data
datPath=tollgate/hourlyAscii/

# define file path for the vegetation data
vegPath=/home/mclark/summa/settings/tollgate/localAttributes/

# define file path for the forcing file list
infoPath=/home/mclark/summa/settings/tollgate/forcingFileInfo/

# define pattern to replace with actual data file
cPattern='XX_forcingFileName_XX'

# loop through files in the local attributes folder
for vegFile in $( ls  ${vegPath} ); do

 # split the string using underscores
 IFS='_' read -a strarr <<< "${vegFile}"

 # define filename (add the last element from the files in ${file_path})
 infoFilenm=snow_zForcingFileInfo_basin_${strarr[-1]}

 # split last string element using periods (to remove ".txt")
 IFS='.' read -a suffix <<< "${strarr[-1]}"

 # define the basin
 iBasin=${suffix[0]}
 
 # define the forcing data file
 datFilenm=basinForcing.$iBasin.txt

 # make the forcing info file
 cp snow_zForcingInfo_template.txt $infoPath$infoFilenm

 # replace the search string with the actual filename
 sed -i 's/'${cPattern}'/'${datFilenm}'/g' $infoPath$infoFilenm

 # print results
 echo filenm: $infoFilenm - $iBasin - $datFilenm
 
done
