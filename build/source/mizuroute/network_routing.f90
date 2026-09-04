module network_routing_module

  use nrtype

  use dataTypes,       only: mizu_remap       => remap
  use dataTypes,       only: mizu_runoff      => runoff
  use mizuroute_types, only: mizuroute_topology
  use mizuroute_types, only: river_network_data
  use mizuroute_types, only: spatial_remap_data

  use var_lookup,      only: ixNTOPO            ! index of variables for the network topology

  implicit none

  private
  public :: network_routing

contains

  ! -----------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------

  subroutine network_routing(sub_idx, river_network, routing_map, do_remapping, &
                             ierr, message)

    use process_remap_module, only: remap_runoff
    use process_remap_module, only: basin2reach

    USE globalData,           only: rch_routes  ! routing methods instantiated

    implicit none

    ! input
    integer(i4b)             , intent(in)      :: sub_idx
    type(river_network_data) , intent(inout)   :: river_network
    type(mizu_remap)         , intent(in)      :: routing_map
    logical(lgt)             , intent(in)      :: do_remapping

    ! output
    integer(i4b)             , intent(out)     :: ierr
    character(*)             , intent(out)     :: message

    ! locals
    integer(i4b) :: ix
    integer(i4b) :: iSeg,jSeg
    integer(i4b) :: n_seg
    real(dp)     :: T0,T1
    real(dp)     :: fracStep
    integer(i4b) :: iSub
    character(len=256)  :: cmessage

    ! initialize error control
    ierr    = 0
    message = 'network_routing/'

    n_seg = river_network%core%topology%n_seg

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
    ! network routing
    !---------------------------------------------------------------------

    fracStep = 1._dp / real(river_network%driver%time%n_sub, dp)

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

        ! * loop through stream segments
        do iSeg = 1, n_seg
  
          ! process segments in the prescribed upstream-to-downstream routing order
          jSeg = river_network%core%topology%ntopo(iSeg)%var(ixNTOPO%rchOrder)%dat(1)
   
          ! route runoff for segment jSeg using method ix
          call rch_route%route(jSeg,                       &
                               T0, T1,                     &
                               river_network%core%ntopo,   &
                               river_network%core%param,   &
                               river_network%core%state,   &
                               river_network%core%flux,    &
                               ierr, cmessage)
          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   
          ! aggregate streamflow per substep
          river_network%driver%method(ix)%streamflow(jSeg, sub_idx) = river_network%driver%method(ix)%streamflow(jSeg, sub_idx) + &
                                                                      river_network%core%flux(jSeg)%ROUTE(ix)%REACH_Q * fracStep

        end do  ! * loop through stream segments

        !stop 'check'

      end do  ! * loop through substeps
   
      end associate
   
    end do ! * loop through routing methods

  end subroutine network_routing

end module network_routing_module
