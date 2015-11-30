#!/bin/bash
#
# used to make a list of forcing files for each GRU
# (based on the localAttributes files)

# define the summa root
summaRoot=/home/mclark/summa

# define file path for the vegetation data
dataPath=${summaRoot}/input/plumber/

# define the path to the decisions file
decisionsPath=${summaRoot}/settings/plumber/spinupDecisions/

# define pattern to replace with actual data file
cPatternStart='XX_startYear_XX'
cPatternEnd='XX_endYear_XX'

# loop through files in the local attributes folder
for dataFile in $( ls  ${dataPath} ); do

 # split the string using underscores
 IFS='_' read -a strarr <<< "${dataFile}"

 # define decisions filename (add the last element from the files in ${file_path})
 decisionsFilenm=summa_zDecisions_${strarr[-1]}

 # get the first line of a data file
 cHead=$(head -n1 $dataPath$dataFile)

 # extract the first "word" (year)
 IFS=' ' read -a dataLine <<< "${cHead}"
 startYear=${dataLine[0]}

 # get the end year
 endYear=$(($startYear+1))

 # make the decisions file
 cp templates/summa_zDecisions_template.txt $decisionsPath$decisionsFilenm

 # replace the search string with the start and end year
 sed -i 's/'${cPatternStart}'/'${startYear}'/g' $decisionsPath$decisionsFilenm
 sed -i 's/'${cPatternEnd}'/'${endYear}'/g' $decisionsPath$decisionsFilenm

 # print progress
 echo $dataFile $startYear $endYear
 
done
