module wq_exchange_module

  ! Reach water budget and runoff delivery map handed to a water-quality coupler.

  use nrtype
  use dataTypes,       only: mizu_remap => remap
  use mizuroute_types, only: river_network_data
  use var_lookup,      only: ixNTOPO            ! index of variables for the network topology

  implicit none

  private
  public :: init_wq_exchange

contains

  ! -----------------------------------------------------------------------------------------------
  ! Allocate the exchange arrays and build the map from runoff elements to reaches.
  ! map_area(k) is the reach inflow (m3 s-1) produced by unit runoff (1 m s-1) from one element,
  ! obtained by passing unit runoff through the same remapping and aggregation used for routing.
  ! -----------------------------------------------------------------------------------------------
  subroutine init_wq_exchange(river_network, routing_map, do_remapping, ierr, message)

    use process_remap_module, only: remap_runoff
    use process_remap_module, only: basin2reach
    use public_var,           only: iulog

    implicit none

    type(river_network_data), intent(inout) :: river_network
    type(mizu_remap),         intent(in)    :: routing_map
    logical(lgt),             intent(in)    :: do_remapping
    integer(i4b),             intent(out)   :: ierr
    character(*),             intent(out)   :: message

    integer(i4b)          :: n_seg, n_q, iq, iSeg, n_map, n_unmapped
    real(dp), allocatable :: sim_save(:), basin_save(:), reach_q(:)
    integer(i4b), allocatable :: tmp_reach(:)
    real(dp),     allocatable :: tmp_area(:)
    character(len=256)    :: cmessage

    ierr    = 0
    message = 'init_wq_exchange/'

    associate(core => river_network%core, wq => river_network%driver%wq)

    n_seg = core%topology%n_seg
    n_q   = size(core%runoff%sim)

    allocate(wq%order(n_seg), wq%down_index(n_seg), wq%vol_start(n_seg), wq%vol_end(n_seg),         &
             wq%vol_lateral(n_seg), wq%vol_upstream(n_seg), wq%vol_outflow(n_seg), wq%vol_wm(n_seg),    &
             wq%map_start(n_q+1), &
             reach_q(n_seg), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating the exchange arrays'; return; endif

    do iSeg = 1, n_seg
      wq%order(iSeg)      = core%topology%ntopo(iSeg)%var(ixNTOPO%rchOrder)%dat(1)
      wq%down_index(iSeg) = core%ntopo(iSeg)%DREACHI
      wq%vol_start(iSeg)  = core%flux(iSeg)%ROUTE(1)%REACH_VOL(1)
    end do
    wq%vol_end      = wq%vol_start
    wq%vol_lateral  = 0._dp
    wq%vol_upstream = 0._dp
    wq%vol_outflow  = 0._dp
    wq%vol_wm       = 0._dp

    ! unit runoff from one element at a time
    sim_save   = core%runoff%sim
    basin_save = core%runoff%basinRunoff
    allocate(tmp_reach(0), tmp_area(0))
    do iq = 1, n_q
      core%runoff%sim     = 0._dp
      core%runoff%sim(iq) = 1._dp
      if(do_remapping)then
        call remap_runoff(core%runoff, routing_map, core%runoff%basinRunoff, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      else
        core%runoff%basinRunoff = core%runoff%sim
      endif
      call basin2reach(core%runoff%basinRunoff, core%ntopo, core%param, reach_q, ierr, cmessage, limitRunoff=.false.)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      wq%map_start(iq) = size(tmp_reach) + 1
      do iSeg = 1, n_seg
        if(reach_q(iSeg) > 0._dp)then
          tmp_reach = [tmp_reach, iSeg]
          tmp_area  = [tmp_area,  reach_q(iSeg)]
        endif
      end do
    end do
    n_map = size(tmp_reach)
    wq%map_start(n_q+1) = n_map + 1
    ! runoff elements that reach no segment: their water quality leaves the domain with their runoff
    n_unmapped = count([(wq%map_start(iq+1) == wq%map_start(iq), iq=1,n_q)])
    if(n_unmapped > 0) write(iulog,'(A,I0,A)') ' WARNING: water quality: ', n_unmapped, &
      ' runoff element(s) deliver to no river reach (no river HRU maps to them); their solute leaves the domain'
    call move_alloc(tmp_reach, wq%map_reach)
    call move_alloc(tmp_area,  wq%map_area)
    core%runoff%sim         = sim_save
    core%runoff%basinRunoff = basin_save

    wq%active = .true.

    end associate

  end subroutine init_wq_exchange

end module wq_exchange_module
