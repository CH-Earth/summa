! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

program multi_driver
! used to evaluate different methods for simulating snow processes
! *****************************************************************************
! use desired modules
! *****************************************************************************
USE nrtype                                                  ! variable types, etc.
! provide access to subroutines and functions
USE summaFileManager,only:summa_SetDirsUndPhiles            ! sets directories and filenames
USE module_sf_noahmplsm,only:read_mp_veg_parameters         ! module to read NOAH vegetation tables
USE module_sf_noahmplsm,only:redprm                         ! module to assign more Noah-Mp parameters
USE allocspace_module,only:init_metad                       ! module to allocate space for metadata structures
USE allocspace_module,only:alloc_stim                       ! module to allocate space for scalar time structures
USE allocspace_module,only:alloc_time                       ! module to allocate space for model time structures
USE allocspace_module,only:alloc_forc                       ! module to allocate space for model forcing data strictures
USE allocspace_module,only:alloc_mpar                       ! module to allocate space for local column model parameter structures
USE allocspace_module,only:alloc_mvar                       ! module to allocate space for local column model variable structures
USE allocspace_module,only:alloc_indx                       ! module to allocate space for local column model indices
USE allocspace_module,only:alloc_bpar                       ! module to allocate space for basin-average model parameter structures
USE allocspace_module,only:alloc_bvar                       ! module to allocate space for basin-average model variable structures
USE mDecisions_module,only:mDecisions                       ! module to read model decisions
USE read_metad_module,only:read_metad                       ! module to populate metadata structures
USE def_output_module,only:def_output                       ! module to define model output
USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
USE read_attrb_module,only:read_attrb                       ! module to read local attributes
USE read_pinit_module,only:read_pinit                       ! module to read initial model parameter values
USE paramCheck_module,only:paramCheck                       ! module to check consistency of model parameters
USE read_icond_module,only:read_icond                       ! module to read initial conditions
USE read_param_module,only:read_param                       ! module to read model parameter sets
USE ConvE2Temp_module,only:E2T_lookup                       ! module to calculate a look-up table for the temperature-enthalpy conversion
USE var_derive_module,only:calcHeight                       ! module to calculate height at layer interfaces and layer mid-point
USE var_derive_module,only:v_shortcut                       ! module to calculate "short-cut" variables
USE var_derive_module,only:rootDensty                       ! module to calculate the vertical distribution of roots
USE var_derive_module,only:satHydCond                       ! module to calculate the saturated hydraulic conductivity in each soil layer
USE var_derive_module,only:fracFuture                       ! module to calculate the fraction of runoff in future time steps (time delay histogram)
USE read_force_module,only:read_force                       ! module to read model forcing data
USE derivforce_module,only:derivforce                       ! module to compute derived forcing data
USE modelwrite_module,only:writeAttrb,writeParam            ! module to write model attributes and parameters
USE modelwrite_module,only:writeForce                       ! module to write model forcing data
USE modelwrite_module,only:writeModel,writeBasin            ! module to write model output
USE coupled_em_module,only:coupled_em                       ! module to run the coupled energy and mass model
USE groundwatr_module,only:groundwatr                       ! module to simulate regional groundwater balance
USE qTimeDelay_module,only:qOverland                        ! module to route water through an "unresolved" river network
! provide access to data
USE summaFileManager,only:SETNGS_PATH                       ! define path to settings files (e.g., Noah vegetation tables)
USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
USE summaFileManager,only:LOCALPARAM_INFO,BASINPARAM_INFO   ! files defining the default values and constraints for model parameters
USE data_struc,only:doJacobian                              ! flag to compute the Jacobian
USE data_struc,only:localParFallback                        ! local column default parameters
USE data_struc,only:basinParFallback                        ! basin-average default parameters
USE data_struc,only:mpar_meta,bpar_meta                     ! metadata for local column and basin-average model parameters
USE data_struc,only:numtim                                  ! number of time steps
USE data_struc,only:time_data,time_hru,refTime              ! time and reference time
USE data_struc,only:forc_data,forc_hru                      ! model forcing data
USE data_struc,only:type_data,type_hru                      ! classification of veg, soils etc.
USE data_struc,only:attr_data,attr_hru                      ! local attributes (lat, lon, elev, etc.)
USE data_struc,only:mpar_data,mpar_hru                      ! local column model parameters
USE data_struc,only:mvar_data,mvar_hru                      ! local column model variables
USE data_struc,only:indx_data,indx_hru                      ! local column model indices
USE data_struc,only:bpar_data                               ! basin-average model parameters
USE data_struc,only:bvar_data                               ! basin-average model variables
USE data_struc,only:model_decisions                         ! model decisions
USE data_struc,only:urbanVegCategory                        ! vegetation category for urban areas
USE data_struc,only:globalPrintFlag                         ! global print flag
USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM                    ! 2-d tables for stem area index and leaf area index (vegType,month)
USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB                      ! height at the top and bottom of vegetation (vegType)
! named variables for elements of model structures
USE var_lookup,only:iLookTIME,iLookFORCE                    ! look-up values for time and forcing data structures
USE var_lookup,only:iLookTYPE                               ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookATTR                               ! look-up values for local attributes
USE var_lookup,only:iLookMVAR                               ! look-up values for local column model variables
USE var_lookup,only:iLookPARAM                              ! look-up values for local column model parameters
USE var_lookup,only:iLookINDEX                              ! look-up values for local column index variables
USE var_lookup,only:iLookBVAR                               ! look-up values for basin-average model variables
USE var_lookup,only:iLookBPAR                               ! look-up values for basin-average model parameters
USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions
! named variables for model decisions
USE mDecisions_module,only:  &                              ! look-up values for method used to compute derivative
 numerical,   & ! numerical solution
 analytical     ! analytical solution
USE mDecisions_module,only:&                                ! look-up values for LAI decisions
 monthlyTable,& ! LAI/SAI taken directly from a monthly table for different vegetation classes
 specified      ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
USE mDecisions_module,only:&                                ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, & ! separate groundwater representation in each local soil column
 singleBasin    ! single groundwater store over the entire basin
implicit none

! *****************************************************************************
! (0) variable definitions
! *****************************************************************************
! define counters
integer(i4b)              :: iHRU,jHRU,kHRU                 ! index of the hydrologic response unit
integer(i4b)              :: nHRU                           ! number of hydrologic response units
integer(i4b)              :: iStep=0                        ! index of model time step
integer(i4b)              :: jStep=0                        ! index of model output
! define the re-start file
logical(lgt)              :: printRestart                   ! flag to print a re-start file
integer(i4b),parameter    :: ixRestart_im=1001              ! named variable to print a re-start file once per month
integer(i4b),parameter    :: ixRestart_id=1002              ! named variable to print a re-start file once per day
integer(i4b),parameter    :: ixRestart_never=1003           ! named variable to print a re-start file never
integer(i4b)              :: ixRestart=ixRestart_never      ! define frequency to write restart files
! define output file
character(len=8)          :: cdate1=''                      ! initial date
character(len=10)         :: ctime1=''                      ! initial time
character(len=64)         :: output_fileSuffix=''           ! suffix for the output file
character(len=256)        :: summaFileManagerFile=''        ! path/name of file defining directories and files
character(len=256)        :: fileout=''                     ! output filename
! define pointers for model indices
integer(i4b),pointer      :: nSnow=>null()                  ! number of snow layers
integer(i4b),pointer      :: nSoil=>null()                  ! number of soil layers
integer(i4b),pointer      :: nLayers=>null()                ! total number of layers
integer(i4b),pointer      :: midSnowStartIndex=>null()      ! start index of the midSnow vector for a given timestep
integer(i4b),pointer      :: midSoilStartIndex=>null()      ! start index of the midSoil vector for a given timestep
integer(i4b),pointer      :: midTotoStartIndex=>null()      ! start index of the midToto vector for a given timestep
integer(i4b),pointer      :: ifcSnowStartIndex=>null()      ! start index of the ifcSnow vector for a given timestep
integer(i4b),pointer      :: ifcSoilStartIndex=>null()      ! start index of the ifcSoil vector for a given timestep
integer(i4b),pointer      :: ifcTotoStartIndex=>null()      ! start index of the ifcToto vector for a given timestep
real(dp),allocatable      :: dt_init(:)                     ! used to initialize the length of the sub-step for each HRU
real(dp),pointer          :: totalArea=>null()              ! total basin area (m2)
! exfiltration
real(dp),parameter        :: supersatScale=0.001_dp         ! scaling factor for the logistic function (-)
real(dp),parameter        :: xMatch = 0.99999_dp            ! point where x-value and function value match (-)
real(dp),parameter        :: safety = 0.01_dp               ! safety factor to ensure logistic function is less than 1
real(dp),parameter        :: fSmall = epsilon(xMatch)       ! smallest possible value to test
real(dp),allocatable      :: upArea(:)                      ! area upslope of each HRU
! general local variables
real(dp)                  :: fracHRU                        ! fractional area of a given HRU (-)
real(dp),allocatable      :: zSoilReverseSign(:)            ! height at bottom of each soil layer, negative downwards (m)
real(dp),dimension(12)    :: greenVegFrac_monthly           ! fraction of green vegetation in each month (0-1)
real(dp),parameter        :: doubleMissing=-9999._dp        ! missing value
! error control
integer(i4b)              :: err=0                          ! error code
character(len=1024)       :: message=''                     ! error message

! *****************************************************************************
! (1) inital priming -- get command line arguments, identify files, etc.
! *****************************************************************************
print*, 'start'
! get the initial time
call date_and_time(cdate1,ctime1)
print*,ctime1
! get command-line arguments for the output file suffix
call getarg(1,output_fileSuffix)
if (len_trim(output_fileSuffix) == 0) then
 print*,'1st command-line argument missing, expect text string defining the output file suffix'; stop
endif
! get command-line argument for the muster file
call getarg(2,summaFileManagerFile) ! path/name of file defining directories and files
if (len_trim(summaFileManagerFile) == 0) then
 print*,'2nd command-line argument missing, expect path/name of muster file'; stop
endif
! set directories and files -- summaFileManager used as command-line argument
call summa_SetDirsUndPhiles(summaFileManagerFile,err,message); call handle_err(err,message)
! initialize the Jacobian flag
doJacobian=.false.

! *****************************************************************************
! (2) read model metadata
! *****************************************************************************
! initialize model metadata structures
call init_metad(err,message); call handle_err(err,message)
! read metadata on all model variables
call read_metad(err,message); call handle_err(err,message)
! read default values and constraints for model parameters (local column, and basin-average)
call read_pinit(LOCALPARAM_INFO,.TRUE., mpar_meta,localParFallback,err,message); call handle_err(err,message)
call read_pinit(BASINPARAM_INFO,.FALSE.,bpar_meta,basinParFallback,err,message); call handle_err(err,message)

! *****************************************************************************
! (3) read information for each HRU and allocate space for data structures
! *****************************************************************************
! read local attributes for each HRU
call read_attrb(nHRU,err,message); call handle_err(err,message)
! allocate space for HRU data structures
! NOTE: attr_hru and type_hru are defined in read_attrb
call alloc_mpar(nHRU,err,message); call handle_err(err,message)
call alloc_mvar(nHRU,err,message); call handle_err(err,message)
call alloc_indx(nHRU,err,message); call handle_err(err,message)
! allocate space for basin data structures
call alloc_bpar(err,message); call handle_err(err,message)
call alloc_bvar(err,message); call handle_err(err,message)
! allocate space for the forcing and time structures
call alloc_forc(nHRU,err,message); call handle_err(err,message)
call alloc_time(nHRU,err,message); call handle_err(err,message)
call alloc_stim(refTime,err,message); call handle_err(err,message)
! allocate space for the time step (recycled for each HRU for subsequent calls to coupled_em)
allocate(dt_init(nHRU),stat=err); call handle_err(err,'problem allocating space for dt_init')

! *****************************************************************************
! (4a) read description of model forcing datafile used in each HRU
! *****************************************************************************
call ffile_info(nHRU,err,message); call handle_err(err,message)

! *****************************************************************************
! (4b) read model decisions
! *****************************************************************************
call mDecisions(err,message); call handle_err(err,message)

! *****************************************************************************
! (5a) read Noah vegetation and soil tables
! *****************************************************************************
! define monthly fraction of green vegetation
!                           J        F        M        A        M        J        J        A        S        O        N        D
greenVegFrac_monthly = (/0.01_dp, 0.02_dp, 0.03_dp, 0.07_dp, 0.50_dp, 0.90_dp, 0.95_dp, 0.96_dp, 0.65_dp, 0.24_dp, 0.11_dp, 0.02_dp/)
! read Noah soil and vegetation tables
call soil_veg_gen_parm(trim(SETNGS_PATH)//'VEGPARM.TBL',                              & ! filename for vegetation table
                       trim(SETNGS_PATH)//'SOILPARM.TBL',                             & ! filename for soils table
                       trim(SETNGS_PATH)//'GENPARM.TBL',                              & ! filename for general table
                       trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision),    & ! classification system used for vegetation
                       trim(model_decisions(iLookDECISIONS%soilCatTbl)%cDecision))      ! classification system used for soils
! read Noah-MP vegetation tables
call read_mp_veg_parameters(trim(SETNGS_PATH)//'MPTABLE.TBL',                         & ! filename for Noah-MP table
                            trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision)) ! classification system used for vegetation
! define urban vegetation category
select case(trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision))
 case('USGS');                     urbanVegCategory=1
 case('MODIFIED_IGBP_MODIS_NOAH'); urbanVegCategory=13
 case default; call handle_err(30,'unable to identify vegetation category')
end select

! *****************************************************************************
! (5b) read trial model parameter values for each HRU, and populate initial data structures
! *****************************************************************************
call read_param(nHRU,err,message); call handle_err(err,message)
bpar_data%var(:) = basinParFallback(:)%default_val

! *****************************************************************************
! (5c) compute derived model variables that are pretty much constant for the basin as a whole
! *****************************************************************************
call fracFuture(err,message); call handle_err(err,message) ! calculate the fraction of runoff in future time steps

! loop through HRUs
do iHRU=1,nHRU

 ! assign the structures to the appropriate HRUs
 attr_data => attr_hru(iHRU)
 type_data => type_hru(iHRU)
 mpar_data => mpar_hru(iHRU)
 mvar_data => mvar_hru(iHRU)
 indx_data => indx_hru(iHRU)

 ! check that the parameters are consistent
 call paramCheck(err,message); call handle_err(err,message)
 ! read description of model initial conditions -- also initializes model structure components
 ! NOTE: at this stage the same initial conditions are used for all HRUs -- need to modify
 call read_icond(err,message); call handle_err(err,message)
 print*, 'aquifer storage = ', mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1)
 ! assign pointers to model layers
 ! NOTE: layer structure is different for each HRU
 nSnow   => indx_data%var(iLookINDEX%nSnow)%dat(1)
 nSoil   => indx_data%var(iLookINDEX%nSoil)%dat(1)
 nLayers => indx_data%var(iLookINDEX%nLayers)%dat(1)
 ! re-calculate height of each layer
 call calcHeight(&
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 mvar_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,message); call handle_err(err,message)
 ! compute derived model variables that are pretty much constant over each HRU
 call E2T_lookup(err,message); call handle_err(err,message) ! calculate a look-up table for the temperature-enthalpy conversion
 call rootDensty(err,message); call handle_err(err,message) ! calculate vertical distribution of root density
 call satHydCond(err,message); call handle_err(err,message) ! calculate saturated hydraulic conductivity in each soil layer
 call v_shortcut(err,message); call handle_err(err,message) ! calculate "short-cut" variables such as volumetric heat capacity
 ! overwrite the vegetation height
 HVT(type_data%var(iLookTYPE%vegTypeIndex)) = mpar_data%var(iLookPARAM%heightCanopyTop)
 HVB(type_data%var(iLookTYPE%vegTypeIndex)) = mpar_data%var(iLookPARAM%heightCanopyBottom)
 ! overwrite the tables for LAI and SAI
 if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
  SAIM(type_data%var(iLookTYPE%vegTypeIndex),:) = mpar_data%var(iLookPARAM%winterSAI)
  LAIM(type_data%var(iLookTYPE%vegTypeIndex),:) = mpar_data%var(iLookPARAM%summerLAI)*greenVegFrac_monthly
 endif
 ! initialize canopy drip
 ! NOTE: canopy drip from the previous time step is used to compute throughfall for the current time step
 mvar_hru(iHRU)%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1) = 0._dp  ! not used
 ! define the filename for model spinup
 write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_spinup'//trim(output_fileSuffix)//'.nc'
 ! define the file if the first parameter set
 if(iHRU==1) then
  call def_output(nHRU,fileout,err,message); call handle_err(err,message)
 endif
 ! write local model attributes and parameters to the model output file
 call writeAttrb(fileout,iHRU,err,message); call handle_err(err,message)
 call writeParam(fileout,iHRU,err,message); call handle_err(err,message)
 ! initialize indices
 indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1

end do  ! (looping through HRUs)

! allocate space for the upslope area
allocate(upArea(nHRU),stat=err); call handle_err(err,'problem allocating space for upArea')

! identify the total basin area (m2)
totalArea => bvar_data%var(iLookBVAR%basin__totalArea)%dat(1)
totalArea = 0._dp
do iHRU=1,nHRU
 totalArea = totalArea + attr_hru(iHRU)%var(iLookATTR%HRUarea)
end do

! compute total area of the upstream HRUS that flow into each HRU
do iHRU=1,nHRU
 upArea(iHRU) = 0._dp
 do jHRU=1,nHRU
  ! check if jHRU flows into iHRU
  if(type_hru(jHRU)%var(iLookTYPE%downHRUindex) ==  type_hru(iHRU)%var(iLookTYPE%hruIndex))then
   upArea(iHRU) = upArea(iHRU) + attr_hru(jHRU)%var(iLookATTR%HRUarea)
  endif   ! (if jHRU is an upstream HRU)
 end do  ! jHRU
end do  ! iHRU

! initialize aquifer storage
! NOTE: this is ugly: need to add capabilities to initialize basin-wide state variables
select case(model_decisions(iLookDECISIONS%spatial_gw)%iDecision)
 case(localColumn)
  bvar_data%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 0._dp  ! not used
 case(singleBasin)
  bvar_data%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 1._dp
  do iHRU=1,nHRU
   mvar_hru(iHRU)%var(iLookMVAR%scalarAquiferStorage)%dat(1) = 0._dp  ! not used
  end do
 case default; call handle_err(20,'unable to identify decision for regional representation of groundwater')
endselect

! initialize time step length for each HRU
do iHRU=1,nHRU
 dt_init(iHRU) = mvar_hru(iHRU)%var(iLookMVAR%dt_init)%dat(1) ! seconds
end do

! initialize time step index
jstep=1

! ****************************************************************************
! (6) loop through time
! ****************************************************************************
do istep=1,numtim

 ! set print flag
 globalPrintFlag=.false.

 ! read a line of forcing data (if not already opened, open file, and get to the correct place)
 ! NOTE: only read data once: if same data used for multiple HRUs, data is copied across
 do iHRU=1,nHRU  ! loop through HRUs
  ! assign pointers to HRUs
  time_data => time_hru(iHRU)
  forc_data => forc_hru(iHRU)
  ! read forcing data
  call read_force(istep,iHRU,err,message); call handle_err(err,message)
 end do  ! (end looping through HRUs)
 print*, time_data%var


 ! *****************************************************************************
 ! (7) create a new NetCDF output file, and write parameters and forcing data
 ! *****************************************************************************
 ! check the start of a new water year
 if(time_data%var(iLookTIME%im)  ==10 .and. &   ! month = October
    time_data%var(iLookTIME%id)  ==1  .and. &   ! day = 1
    time_data%var(iLookTIME%ih)  ==1  .and. &   ! hour = 1
    time_data%var(iLookTIME%imin)==0)then       ! minute = 0
  ! define the filename
  write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_',&
                                 time_data%var(iLookTIME%iyyy),'-',time_data%var(iLookTIME%iyyy)+1,&
                                 trim(output_fileSuffix)//'.nc'
  ! define the file
  call def_output(nHRU,fileout,err,message); call handle_err(err,message)
  ! write parameters for each HRU, and re-set indices
  do iHRU=1,nHRU
   attr_data => attr_hru(iHRU)
   type_data => type_hru(iHRU)
   mpar_data => mpar_hru(iHRU)
   indx_data => indx_hru(iHRU)
   ! write model parameters to the model output file
   call writeAttrb(fileout,iHRU,err,message); call handle_err(err,message)
   call writeParam(fileout,iHRU,err,message); call handle_err(err,message)
   ! re-initalize the indices for midSnow, midSoil, midToto, and ifcToto
   jStep=1
   indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
   indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
   indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
   indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
   indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
   indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1
  end do  ! (looping through HRUs)
 endif  ! if start of a new water year, and defining a new file

 ! initialize runoff variables
 bvar_data%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._dp  ! surface runoff (m s-1)
 bvar_data%var(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._dp  ! outflow from all "outlet" HRUs (those with no downstream HRU)

 ! initialize baseflow variables
 bvar_data%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._dp ! recharge to the aquifer (m s-1)
 bvar_data%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._dp ! baseflow from the aquifer (m s-1)
 bvar_data%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._dp ! transpiration loss from the aquifer (m s-1)

 ! initialize total inflow for each layer in a soil column
 do iHRU=1,nHRU
  mvar_hru(iHRU)%var(iLookMVAR%mLayerColumnInflow)%dat(:) = 0._dp
 end do


 ! ****************************************************************************
 ! (8) loop through HRUs
 ! ****************************************************************************
 do iHRU=1,nHRU

  ! print progress
  !print*, 'iHRU = ', iHRU

  ! assign pointers to HRUs
  time_data => time_hru(iHRU)
  forc_data => forc_hru(iHRU)
  attr_data => attr_hru(iHRU)
  type_data => type_hru(iHRU)
  mpar_data => mpar_hru(iHRU)
  mvar_data => mvar_hru(iHRU)
  indx_data => indx_hru(iHRU)

  ! identify the area covered by the current HRU
  fracHRU =  attr_data%var(iLookATTR%HRUarea) / bvar_data%var(iLookBVAR%basin__totalArea)%dat(1)

  ! assign pointers to model layers
  ! NOTE: layer structure is different for each HRU
  nSnow   => indx_data%var(iLookINDEX%nSnow)%dat(1)
  nSoil   => indx_data%var(iLookINDEX%nSoil)%dat(1)
  nLayers => indx_data%var(iLookINDEX%nLayers)%dat(1)

  ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
  allocate(zSoilReverseSign(nSoil),stat=err); call handle_err(err,'problem allocating space for zSoilReverseSign')
  zSoilReverseSign(1:nSoil) = -mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow+1:nSnow+nSoil)

  ! assign pointers to model indices
  midSnowStartIndex => indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1)
  midSoilStartIndex => indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1)
  midTotoStartIndex => indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1)
  ifcSnowStartIndex => indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1)
  ifcSoilStartIndex => indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1)
  ifcTotoStartIndex => indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1)

  ! get NOAH-MP parameters
  call REDPRM(type_data%var(iLookTYPE%vegTypeIndex),                           & ! vegetation type index
              type_data%var(iLookTYPE%soilTypeIndex),                          & ! soil type
              type_data%var(iLookTYPE%slopeTypeIndex),                         & ! slope type index
              zSoilReverseSign,                                                & ! * not used: height at bottom of each layer [NOTE: negative] (m)
              nSoil,                                                           & ! number of soil layers
              urbanVegCategory)                                                  ! vegetation category for urban areas

  ! overwrite the vegetation height
  HVT(type_data%var(iLookTYPE%vegTypeIndex)) = mpar_data%var(iLookPARAM%heightCanopyTop)
  HVB(type_data%var(iLookTYPE%vegTypeIndex)) = mpar_data%var(iLookPARAM%heightCanopyBottom)

  ! overwrite the tables for LAI and SAI
  if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
   SAIM(type_data%var(iLookTYPE%vegTypeIndex),:) = mpar_data%var(iLookPARAM%winterSAI)
   LAIM(type_data%var(iLookTYPE%vegTypeIndex),:) = mpar_data%var(iLookPARAM%summerLAI)*greenVegFrac_monthly
  endif

  ! define the green vegetation fraction of the grid box (used to compute LAI)
  mvar_data%var(iLookMVAR%scalarGreenVegFraction)%dat(1) = greenVegFrac_monthly(time_data%var(iLookTIME%im))

  ! compute derived forcing variables
  call derivforce(err,message); call handle_err(err,message)

  ! ****************************************************************************
  ! (9) run the model
  ! ****************************************************************************
  ! define the need to calculate the re-start file
  select case(ixRestart)
   case(ixRestart_im);    printRestart = (time_data%var(iLookTIME%id) == 1 .and. time_data%var(iLookTIME%ih) == 1  .and. time_data%var(iLookTIME%imin) == 0)
   case(ixRestart_id);    printRestart = (time_data%var(iLookTIME%ih) == 1 .and. time_data%var(iLookTIME%imin) == 0)
   case(ixRestart_never); printRestart = .false.
   case default; call handle_err(20,'unable to identify option for the restart file')
  end select
  !printRestart = .true.

  ! run the model for a single parameter set and time step
  call coupled_em(printRestart,                    & ! flag to print a re-start file
                  output_fileSuffix,               & ! name of the experiment used in the restart file
                  dt_init(iHRU),                   & ! initial time step
                  err,message)                       ! error control
  call handle_err(err,message)

  kHRU = 0
  ! identify the downslope HRU
  do jHRU=1,nHRU
   if(type_hru(iHRU)%var(iLookTYPE%downHRUindex) ==  type_hru(jHRU)%var(iLookTYPE%hruIndex))then
    if(kHRU==0)then  ! check there is a unique match
     kHRU=jHRU
    else
     call handle_err(20,'multi_driver: only expect there to be one downslope HRU')
    endif  ! (check there is a unique match)
   endif  ! (if identified a downslope HRU)
  end do

  !write(*,'(a,1x,i4,1x,10(f20.10,1x))') 'iHRU, averageColumnOutflow = ', iHRU, mvar_data%var(iLookMVAR%averageColumnOutflow)%dat(:)

  ! add inflow to the downslope HRU
  if(kHRU > 0)then  ! if there is a downslope HRU
   mvar_hru(kHRU)%var(iLookMVAR%mLayerColumnInflow)%dat(:) = mvar_hru(kHRU)%var(iLookMVAR%mLayerColumnInflow)%dat(:) &
                                                              + mvar_data%var(iLookMVAR%averageColumnOutflow)%dat(:)

  ! increment basin column outflow (m3 s-1)
  else
   bvar_data%var(iLookBVAR%basin__ColumnOutflow)%dat(1) = bvar_data%var(iLookBVAR%basin__ColumnOutflow)%dat(1) + &
                                                          sum(mvar_data%var(iLookMVAR%averageColumnOutflow)%dat(:))
  endif

  ! increment basin surface runoff (m s-1)
  bvar_data%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)   = bvar_data%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    + &
                                                           mvar_data%var(iLookMVAR%averageSurfaceRunoff)%dat(1)    * fracHRU

  ! increment basin-average baseflow input variables (m s-1)
  bvar_data%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = bvar_data%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  + &
                                                            mvar_data%var(iLookMVAR%averageSoilDrainage)%dat(1)     * fracHRU
  bvar_data%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvar_data%var(iLookBVAR%basin__AquiferTranspire)%dat(1) + &
                                                            mvar_data%var(iLookMVAR%averageAquiferTranspire)%dat(1) * fracHRU

  ! increment aquifer baseflow -- ONLY if baseflow is computed individually for each HRU
  ! NOTE: groundwater computed later for singleBasin
  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn)then
   bvar_data%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = bvar_data%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  + &
                                                             mvar_data%var(iLookMVAR%averageAquiferBaseflow)%dat(1) * fracHRU
  endif

  ! write the forcing data to the model output file
  call writeForce(fileout,iHRU,jstep,err,message); call handle_err(err,message)

  ! write the model output to the NetCDF file
  call writeModel(fileout,iHRU,jstep,err,message); call handle_err(err,message)
  !if(istep>6) call handle_err(20,'stopping on a specified step: after call to writeModel')

  ! increment the model indices
  midSnowStartIndex = midSnowStartIndex + nSnow
  midSoilStartIndex = midSoilStartIndex + nSoil
  midTotoStartIndex = midTotoStartIndex + nLayers
  ifcSnowStartIndex = ifcSnowStartIndex + nSnow+1
  ifcSoilStartIndex = ifcSoilStartIndex + nSoil+1
  ifcTotoStartIndex = ifcTotoStartIndex + nLayers+1

  ! deallocate height at bottom of each soil layer(used in Noah MP)
  deallocate(zSoilReverseSign,stat=err); call handle_err(err,'problem deallocating space for zSoilReverseSign')

 end do  ! (looping through HRUs)

 ! compute water balance for the basin aquifer
 if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
  call handle_err(20,'multi_driver/bigBucket groundwater code not transferred from old code base yet')
 endif

 ! perform the routing
 call qOverland(&
                ! input
                model_decisions(iLookDECISIONS%subRouting)%iDecision,           &  ! intent(in): index for routing method
                bvar_data%var(iLookBVAR%basin__SurfaceRunoff)%dat(1),           &  ! intent(in): surface runoff (m s-1)
                bvar_data%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea, &  ! intent(in): outflow from all "outlet" HRUs (those with no downstream HRU)
                bvar_data%var(iLookBVAR%basin__AquiferBaseflow)%dat(1),         &  ! intent(in): baseflow from the aquifer (m s-1)
                bvar_data%var(iLookBVAR%routingFractionFuture)%dat,             &  ! intent(in): fraction of runoff in future time steps (m s-1)
                bvar_data%var(iLookBVAR%routingRunoffFuture)%dat,               &  ! intent(in): runoff in future time steps (m s-1)
                ! output
                bvar_data%var(iLookBVAR%averageInstantRunoff)%dat(1),           &  ! intent(out): instantaneous runoff (m s-1)
                bvar_data%var(iLookBVAR%averageRoutedRunoff)%dat(1),            &  ! intent(out): routed runoff (m s-1)
                err,message)                                                       ! intent(out): error control
 call handle_err(err,message)

 ! write basin-average variables
 call writeBasin(fileout,jstep,err,message); call handle_err(err,message)

 ! increment the time index
 jstep = jstep+1

 !stop 'end of time step'

end do  ! (looping through time)

! deallocate space for dt_init and upArea
deallocate(dt_init,upArea,stat=err); call handle_err(err,'unable to deallocate space for dt_init and upArea')

call stop_program('finished simulation')

contains

 ! **************************************************************************************************
 ! private subroutine handle_err: error handler
 ! **************************************************************************************************
 subroutine handle_err(err,message)
 ! used to handle error codes
 USE data_struc,only:mvar_data,mpar_data,indx_data     ! variable data structure
 USE var_lookup,only:iLookMVAR,iLookPARAM,iLookINDEX    ! named variables defining elements in data structure
 implicit none
 ! define dummy variables
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 ! return if A-OK
 if(err==0) return
 ! process error messages
 if (err>0) then
  write(*,'(a)') 'FATAL ERROR: '//trim(message)
 else
  write(*,'(a)') 'WARNING: '//trim(message); print*,'(can keep going, but stopping anyway)'
 endif
 ! dump variables
 print*, 'error, variable dump:'
 if(allocated(dt_init)) print*, 'dt = ', dt_init
 print*, 'istep = ', istep
 if(associated(type_data))then
  print*, 'HRU index = ', type_data%var(iLookTYPE%hruIndex)
 endif
 if(associated(forc_data))then
  print*, 'pptrate            = ', forc_data%var(iLookFORCE%pptrate)
  print*, 'airtemp            = ', forc_data%var(iLookFORCE%airtemp)
 endif
 if(associated(mpar_data))then
  print*, 'theta_res         = ', mpar_data%var(iLookPARAM%theta_res)            ! soil residual volumetric water content (-)
  print*, 'theta_sat         = ', mpar_data%var(iLookPARAM%theta_sat)            ! soil porosity (-)
  print*, 'plantWiltPsi      = ', mpar_data%var(iLookPARAM%plantWiltPsi)         ! matric head at wilting point (m)
  print*, 'soilStressParam   = ', mpar_data%var(iLookPARAM%soilStressParam)      ! parameter in the exponential soil stress function (-)
  print*, 'critSoilWilting   = ', mpar_data%var(iLookPARAM%critSoilWilting)      ! critical vol. liq. water content when plants are wilting (-)
  print*, 'critSoilTranspire = ', mpar_data%var(iLookPARAM%critSoilTranspire)    ! critical vol. liq. water content when transpiration is limited (-)
 endif
 if(associated(mvar_data))then
  if(associated(mvar_data%var(iLookMVAR%scalarSWE)%dat))then
   print*, 'scalarSWE = ', mvar_data%var(iLookMVAR%scalarSWE)%dat(1)
   print*, 'scalarSnowDepth = ', mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)
   print*, 'scalarCanopyTemp = ', mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1)
   print*, 'scalarRainPlusMelt = ', mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)
   write(*,'(a,100(i4,1x))'   ) 'layerType          = ', indx_data%var(iLookINDEX%layerType)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerDepth        = ', mvar_data%var(iLookMVAR%mLayerDepth)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerTemp         = ', mvar_data%var(iLookMVAR%mLayerTemp)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerVolFracIce   = ', mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerVolFracLiq   = ', mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat
   print*, 'mLayerMatricHead   = ', mvar_data%var(iLookMVAR%mLayerMatricHead)%dat
   print*, 'column inflow = ', mvar_data%var(iLookMVAR%mLayerColumnInflow)%dat

  endif
 endif
 print*,'error code = ', err
 if(associated(time_data)) print*, time_data%var
 write(*,'(a)') trim(message)
 stop
 end subroutine handle_err

 ! **************************************************************************************************
 ! private subroutine stop_program: stop program execution
 ! **************************************************************************************************
 subroutine stop_program(message)
 ! used to stop program execution
 implicit none
 ! define dummy variables
 character(*),intent(in)::message
 ! define the local variables
 integer(i4b),parameter :: outunit=6               ! write to screen
 character(len=8)       :: cdate2                  ! final date
 character(len=10)      :: ctime2                  ! final time
 ! get the final date and time
 call date_and_time(cdate2,ctime2)
 ! print initial and final date and time
 write(outunit,*) 'initial date/time = '//'ccyy='//cdate1(1:4)//' - mm='//cdate1(5:6)//' - dd='//cdate1(7:8), &
                                         ' - hh='//ctime1(1:2)//' - mi='//ctime1(3:4)//' - ss='//ctime1(5:10)
 write(outunit,*) 'final date/time   = '//'ccyy='//cdate2(1:4)//' - mm='//cdate2(5:6)//' - dd='//cdate2(7:8), &
                                         ' - hh='//ctime2(1:2)//' - mi='//ctime2(3:4)//' - ss='//ctime2(5:10)
 ! stop with message
 print*,'FORTRAN STOP: '//trim(message)
 stop
 end subroutine

end program multi_driver


 ! **************************************************************************************************
 ! private subroutine SOIL_VEG_GEN_PARM: Read soil, vegetation and other model parameters (from NOAH)
 ! **************************************************************************************************
!-----------------------------------------------------------------
SUBROUTINE SOIL_VEG_GEN_PARM(FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL, MMINLU, MMINSL)
!-----------------------------------------------------------------
  use module_sf_noahlsm, only : shdtbl, nrotbl, rstbl, rgltbl, &
       &                        hstbl, snuptbl, maxalb, laimintbl, &
       &                        bb, drysmc, f11, maxsmc, laimaxtbl, &
       &                        emissmintbl, emissmaxtbl, albedomintbl, &
       &                        albedomaxtbl, wltsmc, qtz, refsmc, &
       &                        z0mintbl, z0maxtbl, &
       &                        satpsi, satdk, satdw, &
       &                        theta_res, theta_sat, vGn_alpha, vGn_n, k_soil, &  ! MPC add van Genutchen parameters
       &                        fxexp_data, lvcoef_data, &
       &                        lutype, maxalb, &
       &                        slope_data, frzk_data, bare, cmcmax_data, &
       &                        cfactr_data, csoil_data, czil_data, &
       &                        refkdt_data, natural, refdk_data, &
       &                        rsmax_data, salp_data, sbeta_data, &
       &                        zbot_data, smhigh_data, smlow_data, &
       &                        lucats, topt_data, slcats, slpcats, sltype

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL
  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message

!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

  OPEN(19, FILE=trim(FILENAME_VEGTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF


  LUMATCH=0

  FIND_LUTYPE : DO WHILE (LUMATCH == 0)
     READ (19,*,END=2002)
     READ (19,*,END=2002)LUTYPE
     READ (19,*)LUCATS,IINDEX

     IF(LUTYPE.EQ.MMINLU)THEN
        WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
        ! CALL wrf_message( mess )
        LUMATCH=1
     ELSE
        call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
        DO LC = 1, LUCATS+12
           read(19,*)
        ENDDO
     ENDIF
  ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(SHDTBL)       < LUCATS .OR. &
       SIZE(NROTBL)       < LUCATS .OR. &
       SIZE(RSTBL)        < LUCATS .OR. &
       SIZE(RGLTBL)       < LUCATS .OR. &
       SIZE(HSTBL)        < LUCATS .OR. &
       SIZE(SNUPTBL)      < LUCATS .OR. &
       SIZE(MAXALB)       < LUCATS .OR. &
       SIZE(LAIMINTBL)    < LUCATS .OR. &
       SIZE(LAIMAXTBL)    < LUCATS .OR. &
       SIZE(Z0MINTBL)     < LUCATS .OR. &
       SIZE(Z0MAXTBL)     < LUCATS .OR. &
       SIZE(ALBEDOMINTBL) < LUCATS .OR. &
       SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
       SIZE(EMISSMINTBL ) < LUCATS .OR. &
       SIZE(EMISSMAXTBL ) < LUCATS ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F')
  ENDIF

  IF(LUTYPE.EQ.MMINLU)THEN
     DO LC=1,LUCATS
        READ (19,*)IINDEX,SHDTBL(LC),                        &
             NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
             SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
             LAIMAXTBL(LC),EMISSMINTBL(LC),             &
             EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
             ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC)
     ENDDO
!
     READ (19,*)
     READ (19,*)TOPT_DATA
     READ (19,*)
     READ (19,*)CMCMAX_DATA
     READ (19,*)
     READ (19,*)CFACTR_DATA
     READ (19,*)
     READ (19,*)RSMAX_DATA
     READ (19,*)
     READ (19,*)BARE
     READ (19,*)
     READ (19,*)NATURAL
  ENDIF
!
2002 CONTINUE

  CLOSE (19)
  IF (LUMATCH == 0) then
     CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
  ENDIF

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_SOILTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICATION = ', TRIM ( MMINSL )
  ! CALL wrf_message( mess )

  LUMATCH=0



  ! MPC add a new soil table
  FIND_soilTYPE : DO WHILE (LUMATCH == 0)
   READ (19,*)
   READ (19,*,END=2003)SLTYPE
   READ (19,*)SLCATS,IINDEX
   IF(SLTYPE.EQ.MMINSL)THEN
     WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
          SLCATS,' CATEGORIES'
     ! CALL wrf_message ( mess )
     LUMATCH=1
   ELSE
    call wrf_message ( "Skipping over SLTYPE = " // TRIM ( SLTYPE ) )
    DO LC = 1, SLCATS
     read(19,*)
    ENDDO
   ENDIF
  ENDDO FIND_soilTYPE
  ! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(BB    ) < SLCATS .OR. &
       SIZE(DRYSMC) < SLCATS .OR. &
       SIZE(F11   ) < SLCATS .OR. &
       SIZE(MAXSMC) < SLCATS .OR. &
       SIZE(REFSMC) < SLCATS .OR. &
       SIZE(SATPSI) < SLCATS .OR. &
       SIZE(SATDK ) < SLCATS .OR. &
       SIZE(SATDW ) < SLCATS .OR. &
       SIZE(WLTSMC) < SLCATS .OR. &
       SIZE(QTZ   ) < SLCATS  ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
  ENDIF

  ! MPC add new soil table
  select case(trim(SLTYPE))
   case('STAS','STAS-RUC')  ! original soil tables
     DO LC=1,SLCATS
        READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case('ROSETTA')          ! new soil table
     DO LC=1,SLCATS
        READ (19,*) IINDEX,&
             ! new soil parameters (from Rosetta)
             theta_res(LC), theta_sat(LC),        &
             vGn_alpha(LC), vGn_n(LC), k_soil(LC), &
             ! original soil parameters
             BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case default
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  end select

2003 CONTINUE

  CLOSE (19)

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_GENERAL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  READ (19,*)
  READ (19,*)
  READ (19,*) NUM_SLOPE

  SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
     CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
  ENDIF

  DO LC=1,SLPCATS
     READ (19,*)SLOPE_DATA(LC)
  ENDDO

  READ (19,*)
  READ (19,*)SBETA_DATA
  READ (19,*)
  READ (19,*)FXEXP_DATA
  READ (19,*)
  READ (19,*)CSOIL_DATA
  READ (19,*)
  READ (19,*)SALP_DATA
  READ (19,*)
  READ (19,*)REFDK_DATA
  READ (19,*)
  READ (19,*)REFKDT_DATA
  READ (19,*)
  READ (19,*)FRZK_DATA
  READ (19,*)
  READ (19,*)ZBOT_DATA
  READ (19,*)
  READ (19,*)CZIL_DATA
  READ (19,*)
  READ (19,*)SMLOW_DATA
  READ (19,*)
  READ (19,*)SMHIGH_DATA
  READ (19,*)
  READ (19,*)LVCOEF_DATA
  CLOSE (19)

!-----------------------------------------------------------------
END SUBROUTINE SOIL_VEG_GEN_PARM
!-----------------------------------------------------------------
