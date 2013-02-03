program multi_driver
! used to evaluate different methods for simulating snow processes
! *****************************************************************************
! use desired modules
! *****************************************************************************
USE nrtype                                                  ! variable types, etc.
USE snow_fileManager,only:fuse_SetDirsUndPhiles             ! sets directories and filenames
USE snow_fileManager,only:SETNGS_PATH                       ! define path to settings files (e.g., Noah vegetation tables)
USE snow_fileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
USE module_sf_noahmplsm,only:read_mp_veg_parameters         ! module to read NOAH vegetation tables
use module_sf_noahmplsm,only:redprm                         ! module to assign more Noah-Mp parameters
USE allocspace_module,only:init_metad                       ! module to allocate space for metadata structures
USE mDecisions_module,only:mDecisions                       ! module to read model decisions
USE read_metad_module,only:read_metad                       ! module to populate metadata structures
USE def_output_module,only:def_output                       ! module to define model output
USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
USE read_attrb_module,only:read_attrb                       ! module to read local attributes
USE read_pinit_module,only:read_pinit                       ! module to read initial model parameter values
USE pOverwrite_module,only:pOverwrite                       ! module to overwrite default parameter values with info from the Noah tables
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
USE modelwrite_module,only:writeAttrb                       ! module to write model attributes
USE modelwrite_module,only:writeParam,writeForce,writeModel ! module to write model output
USE coupled_em_module,only:coupled_em                       ! module to run the coupled energy and mass model
USE data_struc,only:forcFileInfo                            ! information on forcing data file
USE data_struc,only:time_data,forc_data                     ! time and forcing data structures
USE data_struc,only:type_data                               ! classification of veg, soils etc.
USE data_struc,only:mpar_data,mpar_sets                     ! model parameter information
USE data_struc,only:mvar_data                               ! model variable data
USE data_struc,only:indx_data,indx_meta                     ! index data structures
USE data_struc,only:model_decisions                         ! model decisions
USE data_struc,only:urbanVegCategory                        ! vegetation category for urban areas
USE var_lookup,only:iLookTIME,iLookFORCE                    ! look-up values for time and forcing data structures
USE var_lookup,only:iLookTYPE                               ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookMVAR                               ! look-up values for model variables
USE var_lookup,only:iLookINDEX                              ! look-up values for index variables
USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions
implicit none

! *****************************************************************************
! (0) variable definitions
! *****************************************************************************
! define counters
integer(i4b)              :: iParSet=0                      ! loop through parameter sets
integer(i4b)              :: nParSets=0                     ! number of parameter sets
integer(i4b)              :: iStep=0                        ! index of model time step
integer(i4b)              :: jStep=0                        ! index of model output
integer(i4b)              :: iMonth                         ! index of the current month
! define output file
character(len=8)          :: cdate1=''                      ! initial date
character(len=10)         :: ctime1=''                      ! initial time
character(len=32)         :: output_fileSuffix=''           ! suffix for the output file 
character(len=256)        :: fuseFileManager=''             ! path/name of file defining directories and files
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
real(dp)                  :: dt_init=0._dp                  ! used to initialize the length of the sub-step
! general local variables
real(dp),allocatable      :: zSoilReverseSign(:)            ! height at bottom of each soil layer, negative downwards (m)
! error control
integer(i4b)              :: err=0                          ! error code
character(len=512)        :: message=''                     ! error message

! *****************************************************************************
! (1) inital priming -- get command line arguments, identify files, etc.
! *****************************************************************************
! get the initial time
call date_and_time(cdate1,ctime1)
print*,ctime1
! get command-line arguments for the output file suffix
call getarg(1,output_fileSuffix)
if (len_trim(output_fileSuffix) == 0) then
 print*,'1st command-line argument missing, expect text string defining the output file suffix'; stop
endif
! get command-line argument for the muster file
call getarg(2,fuseFileManager) ! path/name of file defining directories and files
if (len_trim(fuseFileManager) == 0) then
 print*,'2nd command-line argument missing, expect path/name of muster file'; stop
endif
! set directories and files -- fuseFileManager used as command-line argument
call fuse_SetDirsUndPhiles(fuseFileManager,err,message); call handle_err(err,message)
! read model decisions
call mDecisions(err,message); call handle_err(err,message)

! *****************************************************************************
! (2) read model metadata
! *****************************************************************************
! initialize model metadata structures
call init_metad(err,message); call handle_err(err,message)
! read metadata on all model variables
call read_metad(err,message); call handle_err(err,message)
! read description of model forcing datafile
call ffile_info(err,message); call handle_err(err,message)
! read local attributes
call read_attrb(err,message); call handle_err(err,message)
! read default values and constraints for model parameters
call read_pinit(err,message); call handle_err(err,message) 

! *****************************************************************************
! (3) read Noah vegetation tables, and overwrite default values
! *****************************************************************************
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
! use values from Noah tables to overwrite default parameter values
call pOverwrite(err,message); call handle_err(err,message)

! *****************************************************************************
! (4) read trial model parameter values -- and allocate space for parameter structures
! *****************************************************************************
! read trial model parameter values -- and allocate space for parameter structures
call read_param(nParSets,err,message); call handle_err(err,message)

! *****************************************************************************
! (5) loop through the model parameter sets
! *****************************************************************************
do iParSet=1,nParSets

 ! assign the parameter structure to the appropriate parameter set
 mpar_data => mpar_sets(iParSet)
 ! read description of model initial conditions -- also initializes model structure components
 call read_icond(err,message); call handle_err(err,message)
 ! compute derived model variables that are pretty much constant
 call E2T_lookup(err,message); call handle_err(err,message) ! calculate a look-up table for the temperature-enthalpy conversion
 call rootDensty(err,message); call handle_err(err,message) ! calculate vertical distribution of root density
 call calcHeight(err,message); call handle_err(err,message) ! calculate height at layer interfaces and layer mid-point
 call satHydCond(err,message); call handle_err(err,message) ! calculate saturated hydraulic conductivity in each soil layer
 call fracFuture(err,message); call handle_err(err,message) ! calculate the fraction of runoff in future time steps
 call v_shortcut(err,message); call handle_err(err,message) ! calculate "short-cut" variables such as volumetric heat capacity
 ! define the filename for model spinup
 write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_spinup'//trim(output_fileSuffix)//'.nc'
 ! define the file if the first parameter set
 if(iParSet==1) then
  call def_output(fileout,err,message); call handle_err(err,message)
  call writeAttrb(fileout,err,message); call handle_err(err,message)
 endif
 ! write model parameters to the model output file
 call writeParam(fileout,iParSet,err,message); call handle_err(err,message)

 ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
 nSnow   => indx_data%var(iLookINDEX%nSnow)%dat(1)
 nSoil   => indx_data%var(iLookINDEX%nSoil)%dat(1)
 allocate(zSoilReverseSign(nSoil),stat=err); call handle_err(err,'problemAllocate')
 zSoilReverseSign(1:nSoil) = -mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow+1:nSnow+nSoil)

 ! initialize time step length
 dt_init = 900._dp ! seconds

 ! initialize time step index
 jstep=1

 ! ****************************************************************************
 ! (6) loop through time
 ! ****************************************************************************
 do istep=1,forcFileInfo%numtim

  ! assign pointers to model layers
  nSnow   => indx_data%var(iLookINDEX%nSnow)%dat(1)
  nSoil   => indx_data%var(iLookINDEX%nSoil)%dat(1)
  nLayers => indx_data%var(iLookINDEX%nLayers)%dat(1)

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

  ! ***************************************************************************
  ! (7) read forcing data
  ! ***************************************************************************
  ! read a line of forcing data (if not already opened, open file, and get to the correct place)
  call read_force(istep,err,message); call handle_err(err,message)
  ! compute derived forcing variables
  call derivforce(err,message); call handle_err(err,message)

  ! *****************************************************************************
  ! (8) create a new NetCDF output file, and write parameters and forcing data
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
   ! define the file if the first parameter set
   if(iParSet==1) then
    call def_output(fileout,err,message); call handle_err(err,message)
    call writeAttrb(fileout,err,message); call handle_err(err,message)
   endif
   ! write model parameters to the model output file
   call writeParam(fileout,iParSet,err,message); call handle_err(err,message)
    ! re-initalize the indices for midSnow, midSoil, midToto, and ifcToto
    jStep=1
    indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
    indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1
  endif
  ! write the forcing data to the model output file
  if(iParSet==1) call writeForce(fileout,jstep,err,message); call handle_err(err,message)
  !stop 'FORTRAN STOP: after call to writeForce'

  ! ****************************************************************************
  ! (9) run the model
  ! ****************************************************************************
  print*, time_data%var, nSnow
  ! run the model for a single parameter set and time step
  call coupled_em(dt_init,err,message); call handle_err(err,message) 
  if(istep>0) stop 'FORTRAN STOP: after call to coupled_em'

  ! write the model output to the NetCDF file
  call writeModel(fileout,iParSet,jstep,err,message); call handle_err(err,message)
  !if(istep>6) call handle_err(20,'stopping on a specified step: after call to writeModel')
  

  ! increment the model indices
  midSnowStartIndex = midSnowStartIndex + nSnow
  midSoilStartIndex = midSoilStartIndex + nSoil
  midTotoStartIndex = midTotoStartIndex + nLayers
  ifcSnowStartIndex = ifcSnowStartIndex + nSnow+1 
  ifcSoilStartIndex = ifcSoilStartIndex + nSoil+1 
  ifcTotoStartIndex = ifcTotoStartIndex + nLayers+1 

  ! increment the time index
  jstep = jstep+1

 end do  ! (looping through time)

 ! deallocate height at bottom of each soil layer(used in Noah MP)
 deallocate(zSoilReverseSign,stat=err); call handle_err(err,'problemDeallocate')
 call stop_program('end of first parameter set')

end do  ! (looping through model parameter sets)
call stop_program('finished looping through parameter sets')

contains

 subroutine handle_err(err,message)
 ! used to handle error codes
 USE data_struc,only:mvar_data            ! variable data structure
 USE var_lookup,only:iLookMVAR            ! named variables defining elements in data structure
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
 print*, 'dt = ', dt_init
 print*, 'istep = ', istep
 if(associated(forc_data))then
  print*, 'pptrate            = ', forc_data%var(iLookFORCE%pptrate)
  print*, 'airtemp            = ', forc_data%var(iLookFORCE%airtemp)
 endif
 if(associated(mvar_data))then
  print*, 'scalarRainPlusMelt = ', mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)
  write(*,'(a,100(f11.5,1x))') 'mLayerDepth        = ', mvar_data%var(iLookMVAR%mLayerDepth)%dat
  write(*,'(a,100(f11.5,1x))') 'mLayerTemp         = ', mvar_data%var(iLookMVAR%mLayerTemp)%dat
  write(*,'(a,100(f11.5,1x))') 'mLayerVolFracIce   = ', mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat
  write(*,'(a,100(f11.5,1x))') 'mLayerVolFracLiq   = ', mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat
  print*, 'mLayerMatricHead   = ', mvar_data%var(iLookMVAR%mLayerMatricHead)%dat
 endif
 print*,'error code = ', err
 if(associated(time_data)) print*, time_data%var, nSnow
 write(*,'(a)') trim(message)
 stop
 end subroutine handle_err

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
