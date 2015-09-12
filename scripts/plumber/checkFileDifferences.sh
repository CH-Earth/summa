#!/bin/bash
#
# used to convert the summa output files to the plumber format
#
# define the name of the difference file
differenceFilename=test.nc

# define the variable to compare
varDesire=stomatalConductance

# define the file path
fPath=/d1/mclark/PLUMBER_data/model_output

# define the model instantiations
model_one=SUMMA.1.0.exp.01.009
model_two=SUMMA.1.0.exp.01.test

# loop through files in the master path
for outputFile in $( ls  ${fPath}/${model_one}/*.nc ); do

 # strip out the path
 IFS='/' read -a strTemp <<< "${outputFile}"
 outputName=${strTemp[-1]}

 # strip out the model
 IFS='_' read -a strTemp <<< "${outputName}"
 fileSuffix=${strTemp[-1]}

 # define the filenames
 file_one=${fPath}/${model_one}/${model_one}_${fileSuffix}
 file_two=${fPath}/${model_two}/${model_two}_${fileSuffix}
 echo $file_one

 # create the difference file
 ncbo -O --op_typ=sub -v $varDesire $file_one $file_two $differenceFilename

 # create the stats 
 ncap2 -v -A -s 'min=gsl_stats_min('$varDesire');print(min)' $differenceFilename
 ncap2 -v -A -s 'max=gsl_stats_max('$varDesire');print(max)' $differenceFilename

 # remove the difference file
 rm $differenceFilename

done
