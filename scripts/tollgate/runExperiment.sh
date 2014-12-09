#!/bin/bash

# get arguments
exRun=$1
exName=$2
fileManager=$3
exLog=$4

# make a control file
ctlFile=summaCopy.${exName}.control
touch $ctlFile

# run the model
${exRun} _${exName} ${fileManager} > ${exLog}

# remove the control file
rm $ctlFile

exit 0
