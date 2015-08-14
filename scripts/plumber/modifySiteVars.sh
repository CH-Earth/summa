#!/bin/bash
#
# used to extract site variables from the PLUMBER NetCDF output files
#  and place them in the local model parameter files

# Define the path to the local parameter files
localParamPath=/home/mclark/summa/settings/plumber/localParamTrial/

# Define path to the PLUMBER model output
plumberPath=/d1/mclark/PLUMBER_data/model_output/

# Define desired model
plumberModel=CABLE.2.0

# Define the desired site variable
plumberSiteVar=hc

# define the search string
searchString_hvt=XX_hvt_XX
searchString_hvb=XX_hvb_XX

# define the undesired suffix
undesiredSuffix=Fluxnet.1.4.nc

# loop through the sites
for plumberSite in $( ls ${plumberPath}${plumberModel}/* ); do

 # get the site name
 IFS='_' read -a strarr <<< "${plumberSite}"
 siteName=${strarr[-1]}
 siteName=${siteName%$undesiredSuffix}

 # get the name of the local parameter file
 paramFilenm=${localParamPath}summa_zLocalParamTrial_${siteName}.txt

 # create the local parameter file
 cp templates/summa_zLocalParamTrial_template.txt $paramFilenm

 # get the string that defines the site variable
 varString=`ncks -C -u -v $plumberSiteVar $plumberSite | grep $plumberSiteVar | tail -1`

 # extract the site variable (always last)
 IFS='=' read -a strarr <<< "${varString}"
 siteVar=${strarr[-1]}

 # convert the site variable to a number
 siteVar=$(printf "%09.4f" $siteVar)

 # compute a modified site var
 siteVarDiv=10.0
 siteVarMod=$(echo "$siteVar/$siteVarDiv" | bc -l)
 siteVarMod=$(printf "%09.4f" $siteVarMod)

 # replace the search string with the info name 
 sed -i 's/'${searchString_hvt}'/'${siteVar}'/g' $paramFilenm
 sed -i 's/'${searchString_hvb}'/'${siteVarMod}'/g' $paramFilenm

 echo $plumberSite $siteName $siteVar $siteVarMod

done

