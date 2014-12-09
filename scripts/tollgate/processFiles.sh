#!/bin/bash
#
# used to run SUMMA multiple times
#
# define experiments
ixStart=21
ixEnd=80
#
# define start and end years
iyStart=2002
iyEnd=2007

# define file path
filePath=/home/mclark/summa/output/distributedRME/batchRuns/

# loop through years
for iYear in `seq -w $iyStart $iyEnd`
do

 # define string for the years
 cYear1=$((iYear))
 cYear2=$((iYear + 1))
 prefix=basinRunoff_${cYear1}-${cYear2}
 echo $prefix

 # loop through experiments
 for iExp in `seq -w $ixStart $ixEnd`
 do

  # define name of the file
  fileName=${filePath}${prefix}_trial0${iExp}.nc
  echo ${fileName}

  # define a temporary file name
  fileTemp=${filePath}${prefix}_fTemp0${iExp}.nc

  # fix the time dimension (need to ccreate space for a new record dimension)
  ncks -a --fix_rec_dmn time -O ${fileName} ${fileTemp}

  # add a new record dimension
  ncecat -O ${fileTemp} ${fileTemp}

 done

 # define the new filename
 fileNew=${filePath}${prefix}_allRuns.nc

 # concatenate files along the new record dimension
 ncrcat -O ${filePath}${prefix}_fTemp*.nc ${fileNew}

 # fix the record dimension
 ncks -a --fix_rec_dmn record -O ${fileNew} ${fileNew}
 
 # ensure the record dimension comes first
 ncpdq -O -a time,record ${fileNew} ${fileNew}

 # make the time dimension the record dimension
 ncks -O --mk_rec_dmn time ${fileNew} ${fileNew}

 # remove the temporary files
 rm ${filePath}${prefix}_fTemp*.nc

done
