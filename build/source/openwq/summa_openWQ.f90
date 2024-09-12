module summa_openwq
  USE nrtype
  USE openWQ,only:CLASSWQ_openwq
  USE data_types,only:gru_hru_doubleVec
  implicit none
  private
  ! Subroutines
  public :: openwq_init
  public :: openwq_run_time_start
  public :: openwq_run_space_step
  public :: openwq_run_time_end
  private:: openWQ_run_time_start_inner ! inner call at the HRU level

  ! Global Data for prognostic Variables of HRUs
  type(gru_hru_doubleVec),save,public   :: progStruct_timestep_start ! copy of progStruct at the start of timestep for passing fluxes
  type(CLASSWQ_openwq),save,public      :: openwq_obj
  

  contains

! Initialize the openWQ object
subroutine openwq_init(err)
  USE globalData,only:gru_struc                               ! gru-hru mapping structures
  USE globalData,only:prog_meta
  USE globalData,only:maxLayers,maxSnowLayers
  USE allocspace_progStuct_module,only:allocGlobal_porgStruct ! module to allocate space for global data structures
  implicit none

  ! Dummy Varialbes
  integer(i4b), intent(out)                       :: err

  ! local variables
  integer(i4b)                                    :: hruCount
  integer(i4b)                                    :: nSoil
  ! OpenWQ dimensions
  integer(i4b)                                    :: nCanopy_2openwq =  1    ! Canopy has only 1 layer
  integer(i4b)                                    :: nRunoff_2openwq  = 1   ! Runoff has only 1 layer (not a summa variable - openWQ keeps track of this)
  integer(i4b)                                    :: nAquifer_2openwq = 1   ! GW has only 1 layer
  integer(i4b)                                    :: nYdirec_2openwq  = 1   ! number of layers in the y-dir (not used in summa)
  ! error handling
  character(len=256)                              :: message

  ! nx -> num of HRUs)
  ! ny -> 1
  ! nz -> num of layers (snow + soil)
  openwq_obj = CLASSWQ_openwq() ! initalize openWQ object

  hruCount = sum( gru_struc(:)%hruCount )

  nSoil = maxLayers - maxSnowLayers

  ! intialize openWQ
  err=openwq_obj%decl(    &
    hruCount,             & ! num HRU
    nCanopy_2openwq,      & ! num layers of canopy (fixed to 1)
    maxSnowLayers,        & ! num layers of snow (fixed to max of 5 because it varies)
    nSoil,                & ! num layers of snoil (variable)
    nRunoff_2openwq,      & ! num layers of runoff (fixed to 1)
    nAquifer_2openwq,     & ! num layers of aquifer (fixed to 1)
    nYdirec_2openwq)             ! num of layers in y-dir (set to 1 because not used in summa)

  
  ! Create copy of state information, needed for passing to openWQ with fluxes that require
  ! the previous time_steps volume
  call allocGlobal_porgStruct(prog_meta,progStruct_timestep_start,maxSnowLayers,err,message) 

end subroutine openwq_init
  
! Pass Summa State to openWQ
subroutine openwq_run_time_start(summa1_struc)
  USE summa_type, only: summa1_type_dec            ! master summa data type
  USE var_lookup,only: iLookINDEX 
  implicit none

  ! Dummy Varialbes
  type(summa1_type_dec), intent(in)  :: summa1_struc
  ! local variables
  integer(i4b)                       :: openWQArrayIndex !index into OpenWQ's state structure
  integer(i4b)                       :: iGRU
  integer(i4b)                       :: iHRU
  integer(i4b)                       :: nHRU        ! number of HRUs in the GRU (used in looping)
  integer(i4b)                       :: nSoil
  integer(i4b)                       :: nSnow
  logical(1)                         :: lastHRUFlag
  summaVars: associate(&
      progStruct     => summa1_struc%progStruct             , &
      timeStruct     => summa1_struc%timeStruct             , &
      attrStruct     => summa1_struc%attrStruct             , &
      indxStruct     => summa1_struc%indxStruct             , & 
      nGRU           => summa1_struc%nGRU                     &
  )
  ! ############################

  openWQArrayIndex = 0
  lastHRUFlag = .false.

  do iGRU=1,nGRU
    nHRU = size(progStruct%gru(iGRU)%hru(:))
    do iHRU=1,nHRU
      if (iGRU == nGRU .and. iHRU == nHRU )then
        lastHRUFlag = .true.
      end if

      nSnow = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
      nSoil = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)
      
      call openwq_run_time_start_inner(openWQArrayIndex, iGRU, iHRU, &
                                       summa1_struc,&
                                       nSnow, nSoil, lastHRUFlag)
      
      openWQArrayIndex = openWQArrayIndex + 1 

    end do ! end HRU
  end do ! end GRU

  end associate summaVars
end subroutine


subroutine openWQ_run_time_start_inner(openWQArrayIndex, iGRU, iHRU, &
                                      summa1_struc, nSnow, nSoil, last_hru_flag)
  USE summa_type,only: summa1_type_dec            ! master summa data type
  USE var_lookup,only: iLookPROG  ! named variables for state variables
  USE var_lookup,only: iLookATTR  ! named variables for real valued attribute data structure
  USE var_lookup,only: iLookINDEX 
  USE var_lookup,only: iLookVarType  ! named variables for real valued attribute data structure
  USE var_lookup,only: iLookTIME  ! named variables for time data structure
  USE globalData,only:prog_meta
  USE globalData,only:realMissing
  USE multiconst,only:iden_water        ! intrinsic density of liquid water    (kg m-3)
  implicit none
  integer(i4b), intent(in)           :: openWQArrayIndex !index into OpenWQ's state structure
  integer(i4b), intent(in)           :: iGRU
  integer(i4b), intent(in)           :: iHRU
  type(summa1_type_dec), intent(in)  :: summa1_struc
  integer(i4b), intent(in)           :: nSnow
  integer(i4b), intent(in)           :: nSoil
  logical(1),intent(in)              :: last_hru_flag

  ! local variables
  integer(i4b)                       :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  real(rkind)                        :: canopyWatVol_stateVar_summa_m3     ! OpenWQ State Var
  real(rkind)                        :: sweWatVol_stateVar_summa_m3(nSnow) ! OpenWQ State Var
  real(rkind)                        :: soilTemp_depVar_summa_K(nSoil)     ! OpenWQ State Var
  real(rkind)                        :: soilWatVol_stateVar_summa_m3(nSoil)! OpenWQ State Var
  real(rkind)                        :: soilMoist_depVar_summa_frac(nSoil) ! OpenWQ State Var
  real(rkind)                        :: aquiferWatVol_stateVar_summa_m3    ! OpenWQ State Var
  ! counter variables
  integer(i4b)                       :: ilay
  integer(i4b)                       :: iVar
  integer(i4b)                       :: iDat
  integer(i4b)                       :: index
  integer(i4b)                       :: offset
  ! error handling
  integer(i4b)                       :: err

  summaVars: associate(&
    progStruct                  => summa1_struc%progStruct             , &
    timeStruct                  => summa1_struc%timeStruct             , &
    hru_area_m2                 => summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)                     ,&
    Tair_summa_K                => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanairTemp)%dat(1)     ,& ! air temperature (K)
    scalarCanopyWat_summa_kg_m2 => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)      ,& ! canopy water (kg m-2)
    mLayerDepth_summa_m         => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(:)          ,& ! depth of each layer (m)
    mLayerVolFracWat_summa_frac => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(:)     ,& ! volumetric fraction of total water in each layer  (-)
    Tsoil_summa_K               => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat(:)           ,& ! soil temperature (K) for each layer
    AquiferStorWat_summa_m      => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) & ! aquifer storage (m)
  )

  ! ############################
  ! Update unlayered variables and dependencies (1 layer only)
  ! ############################

  if(Tair_summa_K == realMissing) then
    stop 'Error: OpenWQ requires air temperature (K)'
  endif
  
  ! Vegetation
  ! unit for volume = m3 (summa-to-openwq unit conversions needed)
  ! scalarCanopyWat [kg m-2], so needs to  to multiply by hru area [m2] and divide by water density
  if(scalarCanopyWat_summa_kg_m2 == realMissing) then
    canopyWatVol_stateVar_summa_m3 = 0._rkind
  else
    canopyWatVol_stateVar_summa_m3 = scalarCanopyWat_summa_kg_m2 * hru_area_m2 / iden_water
  endif

  ! Aquifer
  ! unit for volume = m3 (summa-to-openwq unit conversions needed)
  ! scalarAquiferStorage [m], so needs to  to multiply by hru area [m2] only
  if(AquiferStorWat_summa_m == realMissing) then
    stop 'Error: OpenWQ requires aquifer storage (m3)'
  endif 
  aquiferWatVol_stateVar_summa_m3 = AquiferStorWat_summa_m * hru_area_m2
        
  ! ############################
  ! Update layered variables and dependenecies
  ! ############################

  if (nSnow .gt. 0)then
    do ilay = 1, nSnow
      ! Snow
      ! unit for volume = m3 (summa-to-openwq unit conversions needed)
      ! mLayerVolFracIce and mLayerVolFracLiq [-], so needs to  to multiply by hru area [m2] and divide by water density
      ! But needs to account for both ice and liquid, and convert to liquid volumes
      if(mLayerVolFracWat_summa_frac(ilay) /= realMissing) then
        sweWatVol_stateVar_summa_m3(ilay) =                             &
          mLayerVolFracWat_summa_frac(ilay)  * mLayerDepth_summa_m(ilay) * hru_area_m2
      else
        sweWatVol_stateVar_summa_m3(ilay) = 0._rkind
      endif
    enddo ! end snow layers
  endif ! end snow


  do ilay = 1, nSoil
    ! Soil
    ! Tsoil
    ! (Summa in K)
    if(Tsoil_summa_K(nSnow+ilay) == realMissing) then
      stop 'Error: OpenWQ requires soil temperature (K)'
    endif
    soilTemp_depVar_summa_K(ilay) = Tsoil_summa_K(nSnow+ilay)
    
    soilMoist_depVar_summa_frac(ilay) = 0     ! TODO: Find the value for this varaibles
    ! Soil
    ! unit for volume = m3 (summa-to-openwq unit conversions needed)
    ! mLayerMatricHead [m], so needs to  to multiply by hru area [m2]
    if(mLayerVolFracWat_summa_frac(nSnow+ilay) == realMissing) then
      stop 'Error: OpenWQ requires soil water (m3)'
    endif
    soilWatVol_stateVar_summa_m3(ilay) =                & 
        mLayerVolFracWat_summa_frac(nSnow+ilay) * hru_area_m2 * mLayerDepth_summa_m(nSnow+ilay)
  enddo


  ! Copy the prog structure
  do iVar = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var)
    do iDat = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat)
      select case(prog_meta(iVar)%vartype)
        case(iLookVarType%ifcSoil);
          offset = 0
        case(iLookVarType%ifcToto);
          offset = 0
        case default
          offset = 1
      end select         
      do index = offset , size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat) - 1 + offset
        progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iVar)%dat(index) = progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat(index)
      enddo
    end do
  end do

        
  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute
        
  err=openwq_obj%openwq_run_time_start(&
                                       last_hru_flag,                   & 
                                       openWQArrayIndex,                & ! total HRUs
                                       nSnow,                           &
                                       nSoil,                           &
                                       simtime,                         &
                                       soilMoist_depVar_summa_frac,     &                    
                                       soilTemp_depVar_summa_K,         &
                                       Tair_summa_K,                    & ! air temperature (K)
                                       sweWatVol_stateVar_summa_m3,     &
                                       canopyWatVol_stateVar_summa_m3,  &
                                       soilWatVol_stateVar_summa_m3,    &
                                       aquiferWatVol_stateVar_summa_m3)
  end associate summaVars

end subroutine openWQ_run_time_start_inner


subroutine openwq_run_space_step(summa1_struc)
  USE var_lookup,only: iLookPROG  ! named variables for state variables
  USE var_lookup,only: iLookTIME  ! named variables for time data structure
  USE var_lookup,only: iLookFLUX  ! named varaibles for flux data
  USE var_lookup,only: iLookATTR  ! named variables for real valued attribute data structure
  USE var_lookup,only: iLookINDEX 
  USE var_lookup,   only: iLookTYPE          ! look-up values for classification of veg, soils etc. 
  USE summa_type,only: summa1_type_dec            ! master summa data type
  USE data_types,only: var_dlength,var_i
  USE globalData,only: gru_struc
  USE globalData,only: data_step   ! time step of forcing data (s)
  USE globalData,only: realMissing
  USE multiconst,only:&
                        iden_ice,       & ! intrinsic density of ice             (kg m-3)
                        iden_water        ! intrinsic density of liquid water    (kg m-3)
  USE module_sf_noahmplsm,only:isWater       ! Identifier for water land cover type

  implicit none

  type(summa1_type_dec),   intent(in)    :: summa1_struc
  

  integer(i4b)                           :: hru_index ! needed because openWQ saves hrus as a single array
  integer(i4b)                           :: iHRU      ! variable needed for looping
  integer(i4b)                           :: iGRU      ! variable needed for looping
  integer(i4b)                           :: iLayer    ! varaible needed for looping

  integer(i4b)                           :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  integer(i4b)                           :: err

  ! compartment indexes in OpenWQ (defined in the hydrolink)
  integer(i4b)                           :: canopy_index_openwq    = 0
  integer(i4b)                           :: snow_index_openwq      = 1
  integer(i4b)                           :: runoff_index_openwq    = 2
  integer(i4b)                           :: soil_index_openwq      = 3
  integer(i4b)                           :: aquifer_index_openwq   = 4
  integer(i4b)                           :: OpenWQindex_s
  integer(i4b)                           :: OpenWQindex_r
  integer(i4b)                           :: iy_r
  integer(i4b)                           :: iz_r
  integer(i4b)                           :: iy_s
  integer(i4b)                           :: iz_s
  real(rkind)                            :: wflux_s2r
  real(rkind)                            :: wmass_source

  ! Summa to OpenWQ units
  ! PrecipVars
  real(rkind)                            :: scalarRainfall_summa_m3
  real(rkind)                            :: scalarSnowfall_summa_m3
  real(rkind)                            :: scalarThroughfallRain_summa_m3
  real(rkind)                            :: scalarThroughfallSnow_summa_m3
  ! CanopyVars
  real(rkind)                            :: canopyStorWat_kg_m3
  real(rkind)                            :: scalarCanopySnowUnloading_summa_m3
  real(rkind)                            :: scalarCanopyLiqDrainage_summa_m3
  real(rkind)                            :: scalarCanopyTranspiration_summa_m3
  real(rkind)                            :: scalarCanopyEvaporation_summa_m3
  real(rkind)                            :: scalarCanopySublimation_summa_m3
  ! runoff vars
  real(rkind)                            :: scalarRunoffVol_m3
  real(rkind)                            :: scalarSurfaceRunoff_summa_m3
  real(rkind)                            :: scalarInfiltration_summa_m3
  ! Snow_SoilVars
  real(rkind)                            :: mLayerLiqFluxSnow_summa_m3
  real(rkind)                            :: iLayerLiqFluxSoil_summa_m3
  real(rkind)                            :: mLayerVolFracWat_summa_m3
  real(rkind)                            :: scalarSnowSublimation_summa_m3
  real(rkind)                            :: scalarSfcMeltPond_summa_m3
  real(rkind)                            :: scalarGroundEvaporation_summa_m3
  real(rkind)                            :: scalarExfiltration_summa_m3
  real(rkind)                            :: mLayerBaseflow_summa_m3
  real(rkind)                            :: scalarSoilDrainage_summa_m3
  real(rkind)                            :: mLayerTranspire_summa_m3
  ! AquiferVars
  real(rkind)                            :: scalarAquiferBaseflow_summa_m3
  real(rkind)                            :: scalarAquiferRecharge_summa_m3
  real(rkind)                            :: scalarAquiferStorage_summa_m3
  real(rkind)                            :: scalarAquiferTranspire_summa_m3

  summaVars: associate(&
      timeStruct     => summa1_struc%timeStruct             , &
      fluxStruct     => summa1_struc%fluxStruct             , &
      nGRU           => summa1_struc%nGRU)



  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute

  hru_index = 0

  ! Summa does not have a y-direction, 
  ! so the dimension will always be 1
  iy_r = 1 
  iy_s = 1

  do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount
      hru_index = hru_index + 1
      if (summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex) == isWater) cycle

      ! ####################################################################
      ! Associate relevant variables
      ! ####################################################################

      DomainVars: associate( &
        ! General Summa info
        hru_area_m2 => summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)  &
      )

      PrecipVars: associate( &
        ! Precipitation 
        scalarRainfall_summa_kg_m2_s             => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarRainfall)%dat(1)                  ,&
        scalarSnowfall_summa_kg_m2_s             => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSnowfall)%dat(1)                  ,&
        scalarThroughfallRain_summa_kg_m2_s      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallRain)%dat(1)           ,&
        scalarThroughfallSnow_summa_kg_m2_s      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallSnow)%dat(1)            &
      )

      CanopyVars: associate( &     
        ! Canopy           
        scalarCanopyWat_summa_kg_m2              => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)  ,&
        scalarCanopySnowUnloading_summa_kg_m2_s  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)       ,&
        scalarCanopyLiqDrainage_summa_kg_m2_s    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)         ,&
        scalarCanopyTranspiration_summa_kg_m2_s  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)       ,& 
        scalarCanopyEvaporation_summa_kg_m2_s    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)         ,&
        scalarCanopySublimation_summa_kg_m2_s    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySublimation)%dat(1)          &
      )

      RunoffVars: associate(&
        scalarSurfaceRunoff_m_s                  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)             ,&   
        scalarInfiltration_m_s                   => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarInfiltration)%dat(1)               &
      )

      Snow_SoilVars: associate(&
        ! Snow + Soil - Control Volume
        current_nSnow                             => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)             ,&
        current_nSoil                             => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)             ,&
        nSnow                                     => gru_struc(iGRU)%hruInfo(iHRU)%nSnow                                                  ,&
        nSoil                                     => gru_struc(iGRU)%hruInfo(iHRU)%nSoil                                                  ,& 
        ! Layer depth and water frac
        mLayerDepth_summa_m                       => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(:)      ,&
        mLayerVolFracWat_summa_frac               => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(:) ,&
        ! Snow Fluxes
        scalarSnowSublimation_summa_kg_m2_s       => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSnowSublimation)%dat(1)           ,&
        scalarSfcMeltPond_kg_m2                   => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSfcMeltPond)%dat(1)  ,&
        iLayerLiqFluxSnow_summa_m_s               => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%iLayerLiqFluxSnow)%dat(:)               ,&
        
        ! Soil Fluxes
        scalarGroundEvaporation_summa_kg_m2_s     => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarGroundEvaporation)%dat(1)         ,&
        iLayerLiqFluxSoil_summa_m_s               => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%iLayerLiqFluxSoil)%dat(:)               ,&
        scalarExfiltration_summa_m_s              => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarExfiltration)%dat(1)              ,&
        mLayerBaseflow_summa_m_s                  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerBaseflow)%dat(:)                  ,&
        scalarSoilDrainage_summa_m_s              => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSoilDrainage)%dat(1)              ,&
        mLayerTranspire_summa_m_s                 => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerTranspire)%dat(:)                  &
      )

      AquiferVars: associate(&
        ! Aquifer
        scalarAquiferStorage_summa_m              => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1), &
        scalarAquiferRecharge_summa_m_s           => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferRecharge)%dat(1)              , &        
        scalarAquiferBaseflow_summa_m_s           => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)              , &        
        scalarAquiferTranspire_summa_m_s          => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferTranspire)%dat(1)               &
      )

      ! ####################################################################
      ! Converte associate variable units: from SUMMA to OpenWQ units
      ! Here only scalar/unlayered variables
      ! OpenWQ: Volume (m3), Time (sec)
      ! Where: Vol in kg/m2, then convert to m3 by multipling by (hru_area_m2 / iden_water)
      ! Where: Flux in kg/m2/s, then convert to m3/time_step by multiplying by (hru_area_m2 * data_step / iden_water)
      ! ####################################################################

      ! PrecipVars
      scalarRainfall_summa_m3 = scalarRainfall_summa_kg_m2_s               * hru_area_m2 * data_step / iden_water
      scalarSnowfall_summa_m3 = scalarSnowfall_summa_kg_m2_s               * hru_area_m2 * data_step / iden_water
      scalarThroughfallRain_summa_m3 = scalarThroughfallRain_summa_kg_m2_s * hru_area_m2 * data_step / iden_water ! flux
      scalarThroughfallSnow_summa_m3 = scalarThroughfallSnow_summa_kg_m2_s * hru_area_m2 * data_step / iden_water ! flux

      ! CanopyVars
      canopyStorWat_kg_m3 = scalarCanopyWat_summa_kg_m2                             * hru_area_m2 / iden_water ! vol
      scalarCanopySnowUnloading_summa_m3 = scalarCanopySnowUnloading_summa_kg_m2_s  * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopyLiqDrainage_summa_m3 = scalarCanopyLiqDrainage_summa_kg_m2_s      * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopyTranspiration_summa_m3 = scalarCanopyTranspiration_summa_kg_m2_s  * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopyEvaporation_summa_m3 = scalarCanopyEvaporation_summa_kg_m2_s      * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopySublimation_summa_m3 = scalarCanopySublimation_summa_kg_m2_s      * hru_area_m2 * data_step / iden_water ! flux
      
      ! runoff vars
      scalarSurfaceRunoff_summa_m3  = scalarSurfaceRunoff_m_s * hru_area_m2 * data_step
      scalarInfiltration_summa_m3   = scalarInfiltration_m_s  * hru_area_m2 * data_step


      ! Snow_SoilVars (unlayered variables)
      ! Other variables are layered and added below as needed
      scalarSnowSublimation_summa_m3 = scalarSnowSublimation_summa_kg_m2_s      * hru_area_m2 * data_step / iden_water
      scalarGroundEvaporation_summa_m3 = scalarGroundEvaporation_summa_kg_m2_s  * hru_area_m2 * data_step / iden_water
      scalarSfcMeltPond_summa_m3 = scalarSfcMeltPond_kg_m2                      * hru_area_m2 / iden_water
      scalarExfiltration_summa_m3 = scalarExfiltration_summa_m_s                * hru_area_m2 * data_step
      scalarSoilDrainage_summa_m3 = scalarSoilDrainage_summa_m_s                * hru_area_m2 * data_step


      ! AquiferVars
      scalarAquiferStorage_summa_m3 = scalarAquiferStorage_summa_m        * hru_area_m2
      scalarAquiferRecharge_summa_m3 = scalarAquiferRecharge_summa_m_s    * hru_area_m2 * data_step
      scalarAquiferBaseflow_summa_m3 = scalarAquiferBaseflow_summa_m_s    * hru_area_m2 * data_step
      scalarAquiferTranspire_summa_m3 = scalarAquiferTranspire_summa_m_s  * hru_area_m2 * data_step
      
      ! Reset Runoff (it's not tracked by SUMMA, so need to track it here)
      scalarRunoffVol_m3 = 0._rkind                                                                   !initialization of this variable is required to limit the runoff aggreggation to each hru.

      ! ####################################################################
      ! Apply Fluxes
      ! Call RunSpaceStep
      ! ####################################################################

      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 1 Fluxes involving the canopy
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------      
      if(scalarCanopyWat_summa_kg_m2 /= realMissing) then
        
        ! ====================================================
        ! 1.1 precipitation -> canopy 
        ! ====================================================
        ! *Source*:
        ! PRECIP (external flux, so need call openwq_run_space_in) 
        ! *Recipient*: canopy (only 1 z layer)
        OpenWQindex_r = canopy_index_openwq
        iz_r          = 1 
        ! *Flux*: the portion of rainfall and snowfall not throughfall
        wflux_s2r = (scalarRainfall_summa_m3 - scalarThroughfallRain_summa_m3) &
                    + (scalarSnowfall_summa_m3 - scalarThroughfallSnow_summa_m3)
        ! *Call openwq_run_space_in* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space_in(                                   &
                                      simtime,                                &
                                      'PRECIP',                               &
                                      OpenWQindex_r, hru_index, iy_r, iz_r,   &
                                      wflux_s2r)

        ! ====================================================
        ! 1.2 canopy -> upper snow layer or runoff pool
        ! scalarCanopySnowUnloading + scalarCanopyLiqDrainage
        ! ====================================================
        ! *Flux*
        ! snow uloading + liq drainage
        wflux_s2r = scalarCanopySnowUnloading_summa_m3 &
                    + scalarCanopyLiqDrainage_summa_m3
        ! *Source*
        ! canopy (only 1 z layer)
        OpenWQindex_s = canopy_index_openwq
        iz_s          = 1
        wmass_source = canopyStorWat_kg_m3
        ! *Recipient* depends on snow layers
        if (current_nSnow .gt. 0)then
          OpenWQindex_r = snow_index_openwq
          iz_r = 1 ! upper layer
        else
          OpenWQindex_r = runoff_index_openwq
          iz_r = 1 ! (has only 1 layer)
          scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r;
        end if
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                                      &
                                    simtime,                                  &
                                    OpenWQindex_s, hru_index, iy_s, iz_s,     &
                                    OpenWQindex_r, hru_index, iy_r, iz_r,     &
                                    wflux_s2r,  &
                                    wmass_source)
        
        ! ====================================================
        ! 1.3 canopy -> OUT (lost from model) (Evap + Subl)
        ! ====================================================
        ! *Source*:
        ! canopy (only 1 z layer)
        OpenWQindex_s = canopy_index_openwq
        iz_s          = 1
        wmass_source = canopyStorWat_kg_m3
        ! *Recipient*: 
        ! lost from system
        OpenWQindex_r = -1
        iz_r          = -1
        ! *Flux*
        ! transpiration + evaporation + sublimation
        wflux_s2r =  scalarCanopyEvaporation_summa_m3  &
                      + scalarCanopySublimation_summa_m3
        
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                                       &
                                      simtime,                                 &
                                      OpenWQindex_s, hru_index, iy_s, iz_s,    &
                                      OpenWQindex_r, hru_index, iy_r, iz_r,    &
                                      wflux_s2r,  &
                                      wmass_source)

      endif

      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 2. Snow / runoff
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------
      ! do all the snow fluxes

      ! ====================================================
      ! 2.1 precicipitation -> upper snow/runoff layer
      ! scalarThroughfallRain + scalarThroughfallSnow
      ! ====================================================
      ! *Flux*
      ! throughfall rain and snow
      wflux_s2r = scalarThroughfallRain_summa_m3 &
                  + scalarThroughfallSnow_summa_m3
      if (current_nSnow .gt. 0)then
        ! *Source*:
        ! PRECIP (external flux, so need call openwq_run_space_in)
        ! *Recipient*: 
        ! snow+soil (upper layer)
        OpenWQindex_r = snow_index_openwq
        iz_r          = 1
      else
        OpenWQindex_r = runoff_index_openwq
        iz_r          = 1
        scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r ! Needed because runoff volume is not tracked
      end if
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space_in(                                     &
                                    simtime,                                  &
                                    'PRECIP',                                 &
                                    OpenWQindex_r, hru_index, iy_r, iz_r,     &
                                    wflux_s2r                                 &
                                    )

      ! Below fluxes only occur when there is no snow
      if (current_nSnow .gt. 0)then

        ! ====================================================
        ! 2.2 snow -> OUT (lost from model) (sublimation)
        ! ====================================================
        ! *Source*:
        ! snow (upper layer)
        OpenWQindex_s = snow_index_openwq
        iz_s          = 1
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(1) * hru_area_m2 * mLayerDepth_summa_m(1)
        wmass_source              = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! lost from system
        OpenWQindex_r = -1
        iz_r          = -1
        ! *Flux*
        ! snow sublimation
        wflux_s2r = scalarSnowSublimation_summa_m3
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                       &
          simtime,                                      &
          OpenWQindex_s, hru_index, iy_s, iz_s,         &
          OpenWQindex_r, hru_index, iy_r, iz_r,         &
          wflux_s2r,                                    &
          wmass_source)

        ! ====================================================
        ! 2.3 snow internal fluxes
        ! ====================================================
        do iLayer = 1, nSnow-1 ! last layer of snow becomes different fluxes 
          ! *Source*: 
          ! snow(iLayer)
          OpenWQindex_s = snow_index_openwq
          iz_s          = iLayer
          mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer) * hru_area_m2 * mLayerDepth_summa_m(iLayer)
          wmass_source              = mLayerVolFracWat_summa_m3
          ! *Recipient*: 
          ! snow(iLayer+1)
          OpenWQindex_r = snow_index_openwq
          iz_r          = iLayer + 1
          ! *Flux*
          mLayerLiqFluxSnow_summa_m3  = iLayerLiqFluxSnow_summa_m_s(iLayer) * hru_area_m2 * data_step
          wflux_s2r                   = mLayerLiqFluxSnow_summa_m3 
          ! *Call openwq_run_space* if wflux_s2r not 0
          err=openwq_obj%openwq_run_space(                       &
            simtime,                                      &
            OpenWQindex_s, hru_index, iy_s, iz_s,         &
            OpenWQindex_r, hru_index, iy_r, iz_r,         &
            wflux_s2r,                                    &
            wmass_source)
        end do

        ! ====================================================
        ! 2.4 snow drainage from the last soil layer -> runoff
        ! ====================================================
        ! *Flux*
        mLayerLiqFluxSnow_summa_m3  = iLayerLiqFluxSnow_summa_m_s(nSnow) * hru_area_m2 * data_step
        wflux_s2r                   = mLayerLiqFluxSnow_summa_m3 
        ! *Source*: 
        ! snow(nSnow)
        OpenWQindex_s = snow_index_openwq
        iz_s          = iLayer
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow) * hru_area_m2 * mLayerDepth_summa_m(nSnow)
        wmass_source              = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! runoff (has one layer only)
        OpenWQindex_r = runoff_index_openwq
        iz_r          = 1
        scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r;
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                       &
          simtime,                                      &
          OpenWQindex_s, hru_index, iy_s, iz_s,         &
          OpenWQindex_r, hru_index, iy_r, iz_r,         &
          wflux_s2r,                                    &
          wmass_source)
      end if
  
      ! ====================================================
      ! 2.5 snow without a layer -> runoff 
      ! scalarSfcMeltPond should be 0 if this occurs
      ! ====================================================
      ! *Source*
      ! snow (this is the case of snow without layer)

      ! need the if condition to protect from invalid read
      ! if the size of mLayerVolFracWat_summa_frac matches the number of soil layers
      ! then summa is expecting no snow for this HRU over the simulation of the model
      if ((nSnow .gt. 0)) then
        ! *Flux*
        ! snow uloading + liq drainage
        wflux_s2r = scalarSfcMeltPond_summa_m3
        ! *Source*
        OpenWQindex_s = snow_index_openwq
        iz_s          = 1
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow) * hru_area_m2 * mLayerDepth_summa_m(nSnow)
        wmass_source              = mLayerVolFracWat_summa_m3
        ! *Recipient*
        ! runoff (has one layer only)
        OpenWQindex_r = runoff_index_openwq
        iz_r          = 1
        scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r;
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                       &
          simtime,                                      &
          OpenWQindex_s, hru_index, iy_s, iz_s,         &
          OpenWQindex_r, hru_index, iy_r, iz_r,         &
          wflux_s2r,                                    &
          wmass_source)
      endif
      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 3. runoff
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------   

      ! ====================================================
      ! 3.1 infiltration
      ! runoff -> top layer of the soil
      ! ====================================================
      ! *Flux*
      wflux_s2r = scalarInfiltration_summa_m3
      ! *Source*: 
      ! runoff (has 1 layer only)
      OpenWQindex_s = runoff_index_openwq
      iz_s          = 1
      wmass_source  = scalarRunoffVol_m3
      ! *Recipient*: 
      ! soil upper layer
      OpenWQindex_r = soil_index_openwq
      iz_r          = 1
      scalarRunoffVol_m3 = scalarRunoffVol_m3 - wflux_s2r;
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                       &
        simtime,                                      &
        OpenWQindex_s, hru_index, iy_s, iz_s,         &
        OpenWQindex_r, hru_index, iy_r, iz_r,         &
        wflux_s2r,                                    &
        wmass_source)

      ! ====================================================
      ! 3.2 surface runoff
      ! runoff -> OUT lost from the system
      ! ====================================================
      ! *Source*: 
      ! runoff (has only 1 layer)
      OpenWQindex_s = runoff_index_openwq
      iz_s          = 1
      wmass_source  = scalarRunoffVol_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r          = -1
      ! *Flux*
      ! wflux_s2r = scalarSurfaceRunoff_summa_m3 
      wflux_s2r = scalarRunoffVol_m3
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                       &
        simtime,                                      &
        OpenWQindex_s, hru_index, iy_s, iz_s,         &
        OpenWQindex_r, hru_index, iy_r, iz_r,         &
        wflux_s2r,                                    &
        wmass_source)
      
      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 4. soil
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------  

      ! ====================================================
      ! 4.1 soil fluxes
      ! upper soil -> OUT (lost from system) (ground evaporation)
      ! ====================================================
      ! *Source*: 
      ! upper soil layer
      OpenWQindex_s = soil_index_openwq
      iz_s          = 1
      mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow+1) * hru_area_m2 * mLayerDepth_summa_m(nSnow+1)
      wmass_source              = mLayerVolFracWat_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r  = -1
      iz_r           = -1
      ! *Flux*
      wflux_s2r = scalarGroundEvaporation_summa_m3
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                       &
        simtime,                                      &
        OpenWQindex_s, hru_index, iy_s, iz_s,         &
        OpenWQindex_r, hru_index, iy_r, iz_r,         &
        wflux_s2r,                                    &
        wmass_source)

      ! ====================================================
      ! 4.2 exfiltration
      ! Lost from the system (first soil layer)
      ! ====================================================
      ! *Source*: 
      ! upper soil layer
      OpenWQindex_s = soil_index_openwq
      iz_s          = 1
      mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow+1) * hru_area_m2 * mLayerDepth_summa_m(nSnow+1)
      wmass_source              = mLayerVolFracWat_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r          = -1
      ! *Flux*
      wflux_s2r = scalarExfiltration_summa_m3
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                       &
        simtime,                                      &
        OpenWQindex_s, hru_index, iy_s, iz_s,         &
        OpenWQindex_r, hru_index, iy_r, iz_r,         &
        wflux_s2r,                                    &
        wmass_source)

      ! ====================================================
      ! 4.3 mLayerBaseflow
      ! Lost from the system at each soil layer
      ! ====================================================
      do iLayer = 1, nSoil
      
        ! *Source*:  
        ! each soil layer
        OpenWQindex_s = soil_index_openwq
        !iz_s          = nSnow + iLayer
        iz_s          = iLayer
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer+nSnow) * hru_area_m2 * mLayerDepth_summa_m(iLayer+nSnow)
        wmass_source              = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! lost from system
        OpenWQindex_r = -1
        iz_r          = -1
        ! *Flux*
        mLayerBaseflow_summa_m3 = mLayerBaseflow_summa_m_s(iLayer) * hru_area_m2 * data_step
        if (iLayer == 1)then
          mLayerBaseflow_summa_m3 = mLayerBaseflow_summa_m3 - scalarExfiltration_summa_m3
        endif
        wflux_s2r = mLayerBaseflow_summa_m3
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                       &
          simtime,                                      &
          OpenWQindex_s, hru_index, iy_s, iz_s,         &
          OpenWQindex_r, hru_index, iy_r, iz_r,         &
          wflux_s2r,                                    &
          wmass_source)
      end do

      ! ====================================================
      ! 4.4 transpiration from the soil
      ! Lost from the system
      ! ====================================================
      do iLayer = 1, nSoil
        ! *Source*:  
        ! all soil layers
        OpenWQindex_s = soil_index_openwq
        iz_s          = iLayer
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer+nSnow) * hru_area_m2 * mLayerDepth_summa_m(iLayer+nSnow)
        wmass_source              = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! lost from system
        OpenWQindex_r = -1
        iz_r          = -1
        mLayerTranspire_summa_m3 = mLayerTranspire_summa_m_s(iLayer) * hru_area_m2 * data_step
        ! *Flux*
        wflux_s2r = mLayerTranspire_summa_m3
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                       &
          simtime,                                      &
          OpenWQindex_s, hru_index, iy_s, iz_s,         &
          OpenWQindex_r, hru_index, iy_r, iz_r,         &
          wflux_s2r,                                    &
          wmass_source)
      end do
      
      ! ====================================================
      ! 4.5 soil internal fluxes
      ! ====================================================
      do iLayer = 1, nSoil - 1 ! last layer of soil becomes different fluxes
        ! *Source*:
        ! soil layer iLayer
        OpenWQindex_s = soil_index_openwq
        iz_s          = iLayer
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer+nSnow) * hru_area_m2 * mLayerDepth_summa_m(iLayer+nSnow)
        wmass_source              = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! soi layer iLayer+1
        OpenWQindex_r = soil_index_openwq
        iz_r          = iLayer + 1
        ! *Flux*
        ! flux between soil layer
        iLayerLiqFluxSoil_summa_m3  = iLayerLiqFluxSoil_summa_m_s(iLayer) * hru_area_m2 * data_step
        wflux_s2r                   = iLayerLiqFluxSoil_summa_m3 
        ! *Call openwq_run_space* if wflux_s2r not 0
        err=openwq_obj%openwq_run_space(                       &
          simtime,                                      &
          OpenWQindex_s, hru_index, iy_s, iz_s,         &
          OpenWQindex_r, hru_index, iy_r, iz_r,         &
          wflux_s2r,                                    &
          wmass_source)
      end do
    
      ! ====================================================
      ! 4.6 soil Draianage into the aquifer
      ! ====================================================
      ! *Source*:
      ! lower soil layer
      OpenWQindex_s = soil_index_openwq
      iz_s          = nSoil
      mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSoil) * hru_area_m2 * mLayerDepth_summa_m(nSoil)
      wmass_source              = mLayerVolFracWat_summa_m3
      ! *Recipient*: 
      ! aquifer (has only 1 layer)
      OpenWQindex_r = aquifer_index_openwq
      iz_r          = 1
      ! *Flux*
      ! flux between soil layer (it's -1 because the first layer gets)
      wflux_s2r = scalarSoilDrainage_summa_m3 
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                       &
        simtime,                                      &
        OpenWQindex_s, hru_index, iy_s, iz_s,         &
        OpenWQindex_r, hru_index, iy_r, iz_r,         &
        wflux_s2r,                                    &
        wmass_source)

      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 5 Aquifer Fluxes
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------

      ! ====================================================
      ! 5.1 Aquifer -> OUT (lost from model) (baseflow) 
      ! ====================================================
      ! *Source*: 
      ! aquifer (only 1 z layer)
      OpenWQindex_s = aquifer_index_openwq
      iz_s          = 1
      wmass_source = scalarAquiferStorage_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r          = -1
      ! *Flux*
      wflux_s2r = scalarAquiferBaseflow_summa_m3
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                       &
        simtime,                                      &
        OpenWQindex_s, hru_index, iy_s, iz_s,         &
        OpenWQindex_r, hru_index, iy_r, iz_r,         &
        wflux_s2r,                                    &
        wmass_source)

      ! ====================================================
      ! 5.2 Aquifer -> OUT (lost from model) (transpiration) 
      ! ====================================================
      ! *Source*: 
      ! aquifer (only 1 z layer)
      OpenWQindex_s = aquifer_index_openwq
      iz_s          = 1
      wmass_source = scalarAquiferStorage_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r          = -1
      ! *Flux*
      wflux_s2r = scalarAquiferTranspire_summa_m3
      ! *Call openwq_run_space* if wflux_s2r not 0
      err=openwq_obj%openwq_run_space(                         &
        simtime,                                        &
        OpenWQindex_s, hru_index, iy_s, iz_s,           &
        OpenWQindex_r, hru_index, iy_r, iz_r,           &
        wflux_s2r,                                      & 
        wmass_source)

      end associate AquiferVars
      end associate Snow_SoilVars
      end associate RunoffVars
      end associate CanopyVars
      end associate PrecipVars
      end associate DomainVars
      
    end do
  end do
end associate summaVars
end subroutine openwq_run_space_step


subroutine openwq_run_time_end(summa1_struc)
  USE summa_type, only:summa1_type_dec      ! master summa data type
  USE var_lookup, only:iLookTIME            ! named variables for time data structure
  implicit none

  ! Dummy Varialbes
  type(summa1_type_dec), intent(in)  :: summa1_struc

  ! Local Variables
  integer(i4b)                       :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  integer(i4b)                       :: err ! error control

  summaVars: associate(&
      timeStruct     => summa1_struc%timeStruct       &       
  )

  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute

  err=openwq_obj%openwq_run_time_end(simtime)           ! minute

  end associate summaVars
end subroutine



end module summa_openwq