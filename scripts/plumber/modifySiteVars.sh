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
 echo '*****' $siteName

 # get the name of the local parameter file
 paramFilenm=${localParamPath}summa_zLocalParamTrial_${siteName}.txt

 # create the local parameter file
 cp templates/summa_zLocalParamTrial_template.txt $paramFilenm

 # loop through site variables
 for plumberSiteVar in hc ssat sfc swilt hyds; do

  # get the value for the given site variable
  varStringValue=`ncks -C -u -v $plumberSiteVar $plumberSite | grep $plumberSiteVar | tail -1`
  IFS='=' read -a strarr <<< "${varStringValue}"
  siteVarValue=${strarr[-1]}

  # convert the string to a number
  if [ ${plumberSiteVar} == 'hyds' ]; then
   siteVarValue=$(printf "%10.4E" $siteVarValue)
  else
   siteVarValue=$(printf "%010.4f" $siteVarValue)
  fi

  # define the match string
  searchString=XX_${plumberSiteVar}_XX

  # replace the search string with the info name 
  sed -i 's/'${searchString}'/'${siteVarValue}'/g' $paramFilenm

  # save the canopy height
  if [ ${plumberSiteVar} == 'hc' ]; then
   canopyHeight=$siteVarValue
  fi

  # save the field capacity
  if [ ${plumberSiteVar} == 'sfc' ]; then
   fieldCapacity=$siteVarValue
  fi

  # save the residual volumetric water content
  if [ ${plumberSiteVar} == 'swilt' ]; then
   thetaRes=$siteVarValue
  fi

  # print progress
  echo $plumberSite $siteName $siteVarValue

 done

 # compute the height at the bottom of the canopy
 canopyDivisor=10.0
 canopyBottom=$(echo "$canopyHeight/$canopyDivisor" | bc -l)
 canopyBottom=$(printf "%010.4f" $canopyBottom)
 sed -i 's/XX_hvb_XX/'${canopyBottom}'/g' $paramFilenm

 # compute the critical point of plant wilting and plant transpiration
 fracWilting=0.05
 fracTranspr=0.50
 soilStorage=$(echo "$fieldCapacity - $thetaRes" | bc -l)
 critWilting=$(echo "$thetaRes + $fracWilting*$soilStorage" | bc -l)
 critTranspr=$(echo "$thetaRes + $fracTranspr*$soilStorage" | bc -l)
 critWilting=$(printf "%010.4f" $critWilting)
 critTranspr=$(printf "%010.4f" $critTranspr)
 sed -i 's/XX_critWilt_XX/'${critWilting}'/g' $paramFilenm
 sed -i 's/XX_critTrans_XX/'${critTranspr}'/g' $paramFilenm

 # check
 echo $fieldCapacity $thetaRes $soilStorage $critWilting $critTranspr

done
