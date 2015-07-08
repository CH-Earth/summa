#!/bin/bash
#
# used to create multiple duplicates of the fileManager files
#
# define parameters
nTrial=50
#
# Define file prefix
#fPrefix='summa_fileManager_lumpedTopmodel__'
fPrefix='summa_fileManager_1dRichards__'
#
# Define the pattern to match for the directory
cPatternPath='XXXXX'

# 
# loop through realizations
for exNum in `seq -w $nTrial`
do
 # define experiment
 exName=ex${exNum}
 # get new files
 oldfile=${fPrefix}template.txt
 newfile=${fPrefix}${exName}.txt
 # copy files
 cp ${oldfile} ${newfile}
 echo ${newfile}
 # define parameter trial
 sed -i 's/'${cPatternPath}'/'${exName}'/g' ${newfile}
done
