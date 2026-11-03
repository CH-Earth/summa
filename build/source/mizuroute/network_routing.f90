module network_routing_module

  use nrtype

  use dataTypes,       only: mizu_remap       => remap
  use dataTypes,       only: mizu_runoff      => runoff
  use mizuroute_types, only: mizuroute_topology
  use mizuroute_types, only: river_network_data
  use mizuroute_types, only: spatial_remap_data

  use var_lookup,      only: ixNTOPO            ! index of variables for the network topology
  use public_var,      only: is_flux_wm, realMissing

  implicit none

  private
  public :: network_routing
  public :: combine_water_hrus

contains

  ! -----------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------

  ! area-weighted lake precipitation and evaporation of the water-class HRUs on each host runoff element, zero without one
  subroutine combine_water_hrus(ixElement, area, precip, evap, host_precip, host_evapo)

    implicit none

    integer(i4b), intent(in)  :: ixElement(:)   ! runoff element holding each water-class HRU
    real(dp),     intent(in)  :: area(:)        ! area of each water-class HRU (m2)
    real(dp),     intent(in)  :: precip(:)      ! precipitation on each water-class HRU (m/s)
    real(dp),     intent(in)  :: evap(:)        ! evaporation from each water-class HRU (m/s)
    real(dp),     intent(out) :: host_precip(:) ! precipitation on each runoff element (m/s)
    real(dp),     intent(out) :: host_evapo(:)  ! evaporation on each runoff element (m/s)
    real(dp)                  :: totArea(size(host_precip))
    integer(i4b)              :: iWater

    totArea = 0._dp; host_precip = 0._dp; host_evapo = 0._dp
    do iWater = 1, size(ixElement)
      totArea(ixElement(iWater))     = totArea(ixElement(iWater))     + area(iWater)
      host_precip(ixElement(iWater)) = host_precip(ixElement(iWater)) + area(iWater)*precip(iWater)
      host_evapo(ixElement(iWater))  = host_evapo(ixElement(iWater))  + area(iWater)*evap(iWater)
    end do
    where (totArea > 0._dp)
      host_precip = host_precip/totArea
      host_evapo  = host_evapo/totArea
    end where

  end subroutine combine_water_hrus

  ! -----------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------

  subroutine network_routing(iStep, simDate, sub_idx, river_network, routing_map, do_remapping, &
                             ierr, message)

    use process_remap_module, only: remap_runoff
    use process_remap_module, only: basin2reach
    use lake_route_module,    only: lake_route

    USE globalData,           only: rch_routes  ! routing methods instantiated
    USE globalData,           only: iTime       ! routing sub-step counter (used by lake_route)
    USE globalData,           only: simDatetime ! current date (used by the Hanasaki and HYPE lake models)
    USE datetime_data,        only: datetime    ! date-time constructor
    USE public_var,           only: is_lake_sim ! flag to route lake segments with the lake models
    USE public_var,           only: calendar    ! calendar name

    implicit none

    ! input
    integer(i4b)             , intent(in)      :: iStep          ! host-model time step since the start of the run
    integer(i4b)             , intent(in)      :: simDate(5)     ! host-model time stamp: year, month, day, hour, minute
    integer(i4b)             , intent(in)      :: sub_idx
    type(river_network_data) , intent(inout)   :: river_network
    type(mizu_remap)         , intent(in)      :: routing_map
    logical(lgt)             , intent(in)      :: do_remapping

    ! output
    integer(i4b)             , intent(out)     :: ierr
    character(*)             , intent(out)     :: message

    ! locals
    integer(i4b) :: ix
    integer(i4b) :: iSeg,jSeg,kSeg
    integer(i4b) :: iUps
    integer(i4b) :: n_seg
    real(dp)     :: T0,T1
    real(dp)     :: fracStep
    integer(i4b) :: iSub
    character(len=256)  :: cmessage

    ! initialize error control
    ierr    = 0
    message = 'network_routing/'

    n_seg = river_network%core%topology%n_seg

    ! the host model owns the clock
    simDatetime(1) = datetime(simDate(1), simDate(2), simDate(3), simDate(4), simDate(5), 0._dp, calendar)

    !---------------------------------------------------------------------
    ! remap runoff from the host land model to river-network HRUs
    !---------------------------------------------------------------------

    if (do_remapping) then
      call remap_runoff(river_network%core%runoff,             &   ! input: routed runoff from the host land model
                        routing_map,                           &   ! input: mapping structure for routing
                        river_network%core%runoff%basinRunoff, &   ! output: runoff for basin HRUs
                        ierr, cmessage)                            ! output: error control
      if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if
    else
      river_network%core%runoff%basinRunoff = river_network%core%runoff%sim
    end if

    ! save basin runoff for host-model output
    river_network%driver%basin_runoff(:,sub_idx) = river_network%core%runoff%basinRunoff(:)

    !---------------------------------------------------------------------
    ! convert basin runoff depth to lateral reach inflow [m3/s]
    !---------------------------------------------------------------------

    ! aggregate basin runoff to each stream segment and convert runoff depth to volumetric lateral inflow
    call basin2reach(river_network%core%runoff%basinRunoff,    & ! input: basin runoff (m/s)
                     river_network%core%ntopo,                 & ! input: reach topology
                     river_network%core%param,                 & ! input: reach parameter
                     river_network%driver%reach_inflow,        & ! output: reach inflow (m3/s)
                     ierr, cmessage)                             ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! transfer lateral inflow (routing workspace) to flux data structures
    do iSeg = 1, n_seg
      river_network%core%flux(iSeg)%BASIN_QR(0) = river_network%core%flux(iSeg)%BASIN_QR(1)        ! streamflow from previous step
      river_network%core%flux(iSeg)%BASIN_QR(1) = river_network%driver%reach_inflow(iSeg)          ! streamflow (m3/s)
    end do

    !---------------------------------------------------------------------
    ! lake precipitation and evaporation [m3/s], mapped like runoff
    !---------------------------------------------------------------------

    if (is_lake_sim) then

      ! save the host-model values for host-model output
      river_network%driver%lake_precip(:,sub_idx) = river_network%driver%host_precip(:)
      river_network%driver%lake_evap(:,sub_idx)   = river_network%driver%host_evapo(:)

      ! remap to river-network HRUs; runoff%sim is reused as the remap input, as in standalone mizuRoute
      if (do_remapping) then
        river_network%core%runoff%sim = river_network%driver%host_precip
        call remap_runoff(river_network%core%runoff, routing_map, river_network%core%runoff%basinPrecip, ierr, cmessage)
        if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if
        river_network%core%runoff%sim = river_network%driver%host_evapo
        call remap_runoff(river_network%core%runoff, routing_map, river_network%core%runoff%basinEvapo, ierr, cmessage)
        if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if
      else
        river_network%core%runoff%basinPrecip = river_network%driver%host_precip
        river_network%core%runoff%basinEvapo  = river_network%driver%host_evapo
      end if

      ! aggregate to reaches and convert depth to volume (reach_inflow is reused as workspace)
      call basin2reach(river_network%core%runoff%basinPrecip, river_network%core%ntopo, river_network%core%param, &
                       river_network%driver%reach_inflow, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      river_network%core%flux(:)%basinPrecip = river_network%driver%reach_inflow(:)

      call basin2reach(river_network%core%runoff%basinEvapo, river_network%core%ntopo, river_network%core%param, &
                       river_network%driver%reach_inflow, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      river_network%core%flux(:)%basinEvapo = river_network%driver%reach_inflow(:)

    end if

    !---------------------------------------------------------------------
    ! network routing
    !---------------------------------------------------------------------

    fracStep = 1._dp / real(river_network%driver%time%n_sub, dp)

    ! start the reach water budget of this host step (first routing method)
    if (river_network%driver%wq%active) then
      do iSeg = 1, n_seg
        river_network%driver%wq%vol_start(iSeg) = river_network%core%flux(iSeg)%ROUTE(1)%REACH_VOL(1)
      end do
      river_network%driver%wq%vol_lateral  = 0._dp
      river_network%driver%wq%vol_upstream = 0._dp
      river_network%driver%wq%vol_outflow  = 0._dp
      river_network%driver%wq%vol_wm       = 0._dp
    end if

    ! * loop through routing methods
    do ix = 1, size(rch_routes)
   
      ! alias the local polymorphic routing object for the selected routing method
      associate(rch_route => rch_routes(ix)%rch_route)
  
      ! initialize streamflow for the land model step
      river_network%driver%method(ix)%streamflow(:,sub_idx) = 0._dp

      ! * loop through substeps
      do iSub = 1, river_network%driver%time%n_sub

        T0 = real(iSub-1, dp) * river_network%driver%time%dt_sub
        T1 = real(iSub,   dp) * river_network%driver%time%dt_sub

        ! lake_route initialises lake storage when iTime == 1
        iTime = (iStep-1)*river_network%driver%time%n_sub + iSub

        ! * loop through stream segments
        do iSeg = 1, n_seg
  
          ! process segments in the prescribed upstream-to-downstream routing order
          jSeg = river_network%core%topology%ntopo(iSeg)%var(ixNTOPO%rchOrder)%dat(1)
   
          ! route runoff for segment jSeg using method ix (lake model for lake segments)
          if (is_lake_sim .and. river_network%core%ntopo(jSeg)%isLake) then
            call lake_route(jSeg, ix,                   &
                            river_network%core%ntopo,   &
                            river_network%core%param,   &
                            river_network%core%flux,    &
                            ierr, cmessage)
          else
            call rch_route%route(jSeg,                       &
                                 T0, T1,                     &
                                 river_network%core%ntopo,   &
                                 river_network%core%param,   &
                                 river_network%core%state,   &
                                 river_network%core%flux,    &
                                 ierr, cmessage)
          endif
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   
          ! aggregate streamflow per substep
          river_network%driver%method(ix)%streamflow(jSeg, sub_idx) = river_network%driver%method(ix)%streamflow(jSeg, sub_idx) + &
                                                                      river_network%core%flux(jSeg)%ROUTE(ix)%REACH_Q * fracStep

          ! storage at the end of the host step
          if (is_lake_sim) river_network%driver%method(ix)%volume(jSeg, sub_idx) = river_network%core%flux(jSeg)%ROUTE(ix)%REACH_VOL(1)

          ! accumulate the reach water budget (upstream reaches are already routed for this substep)
          if (ix == 1 .and. river_network%driver%wq%active) then
            associate(wq => river_network%driver%wq, flux => river_network%core%flux, dt_sub => river_network%driver%time%dt_sub)
              wq%vol_outflow(jSeg) = wq%vol_outflow(jSeg) + flux(jSeg)%ROUTE(ix)%REACH_Q * dt_sub
              wq%vol_lateral(jSeg) = wq%vol_lateral(jSeg) + flux(jSeg)%BASIN_QR(1) * dt_sub
              do iUps = 1, size(river_network%core%ntopo(jSeg)%UREACHI)
                kSeg = river_network%core%ntopo(jSeg)%UREACHI(iUps)
                wq%vol_upstream(jSeg) = wq%vol_upstream(jSeg) + flux(kSeg)%ROUTE(ix)%REACH_Q * dt_sub
              end do
              if (is_flux_wm .and. flux(jSeg)%REACH_WM_FLUX /= realMissing) &
                wq%vol_wm(jSeg) = wq%vol_wm(jSeg) + flux(jSeg)%ROUTE(ix)%REACH_WM_FLUX_actual * dt_sub
              wq%vol_end(jSeg) = flux(jSeg)%ROUTE(ix)%REACH_VOL(1)
            end associate
          end if

        end do  ! * loop through stream segments

        !stop 'check'

      end do  ! * loop through substeps
   
      end associate
   
    end do ! * loop through routing methods

  end subroutine network_routing

end module network_routing_module
