!-----------------------------------------------------------------------
! Compatibility shim for mizuRoute.
! Provides the reach-scale water-balance calculation required by the
! routing routines without the MPI/global-domain dependencies used by
! the full mizuRoute water_balance module.
!-----------------------------------------------------------------------

MODULE water_balance

USE nrtype
! data type
USE dataTypes, ONLY: STRFLX         ! fluxes in each reach
! global parameters
USE public_var, ONLY: iulog         ! i/o logical unit number
USE public_var, ONLY: dt            ! simulation time step

implicit none

private
public:: comp_reach_wb
public:: comp_reach_mb
!public:: comp_global_wb

CONTAINS

  ! *********************************************************************
  ! public subroutine: compute reach water balance
  ! *********************************************************************
  SUBROUTINE comp_reach_wb(seg_id,     &     ! input: reach/lake id
                           ixRoute,    &     ! input: index of routing method
                           Qupstream,  &     ! input: inflow from upstream
                           Qlat,       &     ! input: lateral flow into reach
                           RCHFLX_in,  &     ! inout: reach flux data structure
                           verbose,    &     !
                           lakeFlag,   &     !
                           tolerance)

  ! Descriptions
  ! Compute water balance per simulation time step and reach/lake
  ! During a time step dt = t1 - t0 [sec]
  ! dVol     [m3]   = Vol(t1) - Vol(t0)
  ! flux_in  [m3/s] = inflow(dt) + lateral(dt) + Prec(dt)
  ! flux_out [m3/s] = outflow(dt) + abstract(dt) + Evap(dt)
  !
  ! Compare dVol vs (flux_in - flux_out)*dt

  implicit none
  ! Argument variables:
  integer(i4b), intent(in)                 :: seg_id         ! input: routing method index
  integer(i4b), intent(in)                 :: ixRoute        ! input: routing method index
  real(dp),     intent(in)                 :: Qupstream      ! input: total inflow from upstream reaches
  real(dp),     intent(in)                 :: Qlat           ! input: lateral flow into reach
  type(STRFLX), intent(inout)              :: RCHFLX_in      ! inout: Reach fluxes data structure
  logical(lgt), intent(in)                 :: lakeFlag       ! input: reach index to be examined
  logical(lgt), intent(in)                 :: verbose        ! input: reach index to be examined
  real(dp),     optional, intent(in)       :: tolerance      ! input: wb error tolerance trigering print out
  ! Local variables:
  real(dp)                                 :: wb_tol         !
  real(dp)                                 :: dVol           !
  real(dp)                                 :: Qin            !
  real(dp)                                 :: Qlateral       !
  real(dp)                                 :: Qout           !
  real(dp)                                 :: precip         !
  real(dp)                                 :: evapo          !
  real(dp)                                 :: Qtake_demand   !
  real(dp)                                 :: Qtake_actual   !

  if (present(tolerance)) then
    wb_tol=tolerance
  else
    wb_tol=2.e-5_dp
  end if

  ! volume change
  dVol     = RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(1)-RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(0)
  ! in flux
  Qin      = Qupstream *dt
  Qlateral = Qlat *dt
  if (lakeFlag) then
    precip   = RCHFLX_in%basinprecip *dt
  else
    precip   = 0._dp
  end if
  ! out flux
  Qout         = -1._dp *RCHFLX_in%ROUTE(ixRoute)%REACH_Q *dt
  Qtake_demand = -1._dp *RCHFLX_in%REACH_WM_FLUX *dt
  Qtake_actual = -1._dp *RCHFLX_in%ROUTE(ixRoute)%REACH_WM_FLUX_actual *dt
  if (lakeFlag) then
    evapo    = -1._dp *RCHFLX_in%basinevapo *dt
  else
    evapo   = 0._dp
  end if

  RCHFLX_in%ROUTE(ixRoute)%WB = dVol - (Qin+ Qlateral+ precip+ Qtake_actual+ Qout+ evapo)

  if (verbose .or. abs(RCHFLX_in%ROUTE(ixRoute)%WB) > wb_tol) then
    write(iulog,'(A)')         ' -------------------------------'
    write(iulog,'(A)')         ' -- reach water balance check --'
    write(iulog,'(A)')         ' -------------------------------'
    write(iulog,'(A,1PG15.7)') '  id                  = ', seg_id
    write(iulog,'(A,1PG15.7)') '  lake                = ', lakeFlag
    write(iulog,'(A)')         '  1 = 5-(6+7+8+9+10+12)'
    write(iulog,'(A,1PG15.7)') '  1 WBerr [m3]        = ', RCHFLX_in%ROUTE(ixRoute)%WB
    write(iulog,'(A,1PG15.7)') '  3 Vol at t0 [m3]    = ', RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(0)
    write(iulog,'(A,1PG15.7)') '  4 Vol at t1 [m3]    = ', RCHFLX_in%ROUTE(ixRoute)%REACH_VOL(1)
    write(iulog,'(A,1PG15.7)') '  5 dVol [m3]         = ', dVol
    write(iulog,'(A,1PG15.7)') '  6 inflow [m3]       = ', Qin
    write(iulog,'(A,1PG15.7)') '  7 lateral flow [m3] = ', Qlateral
    write(iulog,'(A,1PG15.7)') '  8 precip [m3]       = ', precip
    write(iulog,'(A,1PG15.7)') '  9 outflow [m3]      = ', Qout
    write(iulog,'(A,1PG15.7)') ' 10 take-actual [m3]  = ', Qtake_actual
    write(iulog,'(A,1PG15.7)') ' 11 take-demand [m3]  = ', Qtake_demand
    write(iulog,'(A,1PG15.7)') ' 12 evaporation [m3]  = ', evapo
  endif
  if (abs(RCHFLX_in%ROUTE(ixRoute)%WB) > wb_tol) then
    write(iulog,'(A,1PG15.7,1X,A,1X,1PG15.7)') ' WARNING: abs. WB error [m3] = ', abs(RCHFLX_in%ROUTE(ixRoute)%WB), '>',wb_tol
  end if

  END SUBROUTINE comp_reach_wb

  ! *********************************************************************
  ! public subroutine: compute reach water balance
  ! *********************************************************************
  SUBROUTINE comp_reach_mb(seg_id,     &     ! input: reach/lake id
                           ixRoute,    &     ! input: index of routing method
                           Cupstream,  &     ! input: mass flux from upstream
                           Clat,       &     ! input: lateral flow into reach
                           RCHFLX_in,  &     ! inout: reach flux data structure
                           verbose,    &     !
                           lakeFlag,   &     !
                           tolerance)
  ! Descriptions
  ! Compute constituent mass balance per simulation time step and reach/lake
  ! During a time step dt = t1 - t0 [sec]
  ! dMass    [mg]   = mass(t1) - mass(t0)
  ! flux_in  [mg/s] = inflow(dt) + lateral(dt)
  ! flux_out [mg/s] = outflow(dt)
  !
  ! Compare dMass vs (flux_in - flux_out)*dt
  implicit none
  ! Argument variables:
  integer(i4b), intent(in)                 :: seg_id         ! input: routing method index
  integer(i4b), intent(in)                 :: ixRoute        ! input: routing method index
  real(dp),     intent(in)                 :: Cupstream      ! input: total inflow from upstream reaches
  real(dp),     intent(in)                 :: Clat           ! input: lateral flow into reach
  type(STRFLX), intent(inout)              :: RCHFLX_in      ! inout: Reach fluxes data structure
  logical(lgt), intent(in)                 :: lakeFlag       ! input: reach index to be examined
  logical(lgt), intent(in)                 :: verbose        ! input: reach index to be examined
  real(dp),     optional, intent(in)       :: tolerance      ! input: wb error tolerance trigering print out
  ! Local variables:
  real(dp)                                 :: MBerr          ! Mass balance error [mg]
  real(dp)                                 :: mb_tol         ! mass balance error tolerance [mg]
  real(dp)                                 :: dMass          ! mass change per time step [mg]
  real(dp)                                 :: Cin            ! mass flux from upstream to a reach [mg]
  real(dp)                                 :: Clateral       ! lateral mass flux [mg]
  real(dp)                                 :: Cout           ! mass flux out of a reach [mg]

  if (present(tolerance)) then
    mb_tol=tolerance
  else
    mb_tol=2.e-5_dp ! in mg
  end if

  ! mass change
  dMass     = RCHFLX_in%ROUTE(ixRoute)%reach_solute_mass(1)-RCHFLX_in%ROUTE(ixRoute)%reach_solute_mass(0)
  ! in flux
  Cin      = Cupstream *dt
  Clateral = Clat *dt
  ! out flux
  Cout         = -1._dp *RCHFLX_in%ROUTE(ixRoute)%reach_solute_flux *dt

  MBerr = dMass - (Cin+ Clateral+ Cout)

  if (verbose) then
    write(iulog,'(A)')         ' -------------------------------------'
    write(iulog,'(A)')         ' -- reach solute mass balance check --'
    write(iulog,'(A)')         ' -------------------------------------'
    write(iulog,'(A,1PG15.7)') '  id                  = ', seg_id
    write(iulog,'(A,1PG15.7)') '  lake                = ', lakeFlag
    write(iulog,'(A)')         '  1 = 5-(6+7+8)'
    write(iulog,'(A,1PG15.7)') '  1 WBerr [mg]        = ', MBerr
    write(iulog,'(A,1PG15.7)') '  3 Mass at t0 [mg]   = ', RCHFLX_in%ROUTE(ixRoute)%reach_solute_mass(0)
    write(iulog,'(A,1PG15.7)') '  4 Mass at t1 [mg]   = ', RCHFLX_in%ROUTE(ixRoute)%reach_solute_mass(1)
    write(iulog,'(A,1PG15.7)') '  5 dMass [mg]        = ', dMass
    write(iulog,'(A,1PG15.7)') '  6 influx [mg]       = ', Cin
    write(iulog,'(A,1PG15.7)') '  7 lateral flux [mg] = ', Clateral
    write(iulog,'(A,1PG15.7)') '  8 outflux [mg]      = ', Cout
  endif
  if (abs(MBerr) > mb_tol) then
    write(iulog,'(A,1PG15.7,1X,A,1X,1PG15.7)') ' WARNING: abs. MB error [m3] = ', abs(MBerr), '>',mb_tol
  end if

  END SUBROUTINE comp_reach_mb

  !     ! *********************************************************************
  !     ! public subroutine: compute global water balance
  !     ! *********************************************************************
  !     SUBROUTINE comp_global_wb(ixRoute, verbose, ierr, message)
  !    
  !       USE globalData, ONLY: masterproc
  !       USE globalData, ONLY: RCHFLX_trib
  !       USE globalData, ONLY: NETOPO_main
  !       USE globalData, ONLY: NETOPO_trib
  !       USE globalData, ONLY: nRch_mainstem        ! scalar data: number of mainstem reaches
  !       USE globalData, ONLY: nRch_trib            ! scalar data: number of tributary reaches
  !       USE globalData, ONLY: nTribOutlet
  !       USE mpi_utils,  ONLY: shr_mpi_reduce
  !    
  !       ! Descriptions
  !       ! Compute water balance per simulation time step over the whole domain
  !       ! During a time step dt = t1 - t0 [sec]
  !       ! dVol     [m3]   = Vol(t1) - Vol(t0)
  !       ! flux_in  [m3/s] = lateral(dt) + Prec(dt)
  !       ! flux_out [m3/s] = abstract(dt) + Evap(dt)
  !       !                   outflow(dt) where outlet
  !       !
  !       ! SUM(dVol) - SUM(flux_in - flux_out)*dt
  !       ! SUM is SUM over the whole domain
  !    
  !       implicit none
  !       ! Argument variables:
  !       integer(i4b),      intent(in)     :: ixRoute        ! input: routing method index
  !       logical(lgt),      intent(in)     :: verbose        ! input: reach index to be examined
  !       integer(i4b),      intent(out)    :: ierr
  !       character(strLen), intent(out)    :: message        ! error message
  !       ! Local variables:
  !       integer(i4b)                      :: lwr,upr        ! loop index
  !       real(dp)                          :: wb_local(7)
  !       real(dp)                          :: wb_global(7)
  !       real(dp)                          :: wb_mainstem(7)
  !       real(dp)                          :: wb_trib(7)
  !       real(dp)                          :: wb_error
  !       character(strLen)                 :: cmessage       ! error message from subroutine
  !    
  !       ierr=0; message='comp_global_wb/'
  !    
  !       if (masterproc) then
  !         wb_trib     = 0._dp
  !         wb_mainstem = 0._dp
  !         if (nRch_mainstem > 0) then
  !           call accum_water_balance(NETOPO_main(1:nRch_mainstem), RCHFLX_trib(1:nRch_mainstem), wb_mainstem, ierr, cmessage)
  !           if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  !         end if
  !         if (nRch_trib > 0) then
  !           lwr = nRch_mainstem + nTribOutlet + 1
  !           upr = nRch_mainstem + nTribOutlet + nRch_trib
  !           call accum_water_balance(NETOPO_trib, RCHFLX_trib(lwr:upr), wb_trib, ierr, cmessage)
  !           if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  !         end if
  !         wb_local = wb_mainstem + wb_trib
  !       else ! non-main processors
  !         call accum_water_balance(NETOPO_trib, RCHFLX_trib, wb_local, ierr, cmessage)
  !         if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  !       end if
  !    
  !       ! sum water balance components across all the processors
  !       wb_global = 0._dp
  !       call shr_mpi_reduce(wb_local, 'sum', wb_global, ierr, cmessage)
  !       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  !    
  !       wb_error = wb_global(1)-sum(wb_global(2:6))
  !    
  !       if (verbose) then
  !         if (masterproc) then
  !           write(iulog,'(A)')         ' ---------------------------'
  !           write(iulog,'(A)')         ' -- global water balance  --'
  !           write(iulog,'(A)')         ' ---------------------------'
  !           write(iulog,'(A)')         ' 8=1-(2+3+4+5+6)'
  !           write(iulog,'(A,1PG15.7)') ' 1 dVol [m3]              = ', wb_global(1)
  !           write(iulog,'(A,1PG15.7)') ' 2 lateral flow [m3]      = ', wb_global(2)
  !           write(iulog,'(A,1PG15.7)') ' 3 precip [m3]            = ', wb_global(3)
  !           write(iulog,'(A,1PG15.7)') ' 4 waterTake-actual [m3]  = ', wb_global(4)
  !           write(iulog,'(A,1PG15.7)') ' 5 evaporation [m3]       = ', wb_global(5)
  !           write(iulog,'(A,1PG15.7)') ' 6 outflow [m3]           = ', wb_global(6)
  !           write(iulog,'(A,1PG15.7)') ' 7 waterTake-demand [m3]  = ', wb_global(7)
  !           write(iulog,'(A,1PG15.7)') ' 8 WBerr [m3]             = ', wb_error
  !         end if
  !       endif
  !       if (abs(wb_error) > 1._dp) then ! tolerance is 1 [m3]
  !         write(iulog,'(A,1PG15.7,1X,A)') ' WARNING: global WB error [m3] = ', wb_error, '> 1.0 [m3]'
  !       end if
  !       flush(iulog)
  !    
  !       CONTAINS
  !    
  !       SUBROUTINE accum_water_balance(NETOPO_in, RCHFLX_in, water_budget, ierr, message)
  !    
  !         USE dataTypes, ONLY: RCHTOPO   ! data structure - Network topology
  !         USE dataTypes, ONLY: STRFLX    ! data structure - fluxes in each reach
  !    
  !         implicit none
  !         ! Arguments:
  !         type(RCHTOPO),     intent(in)  :: NETOPO_in(:)
  !         type(STRFLX),      intent(in)  :: RCHFLX_in(:)
  !         real(dp),          intent(out) :: water_budget(7)
  !         integer(i4b),      intent(out) :: ierr
  !         character(strLen), intent(out) :: message           ! error message
  !         ! Local variables:
  !         integer(i4b)               :: ix          ! loop index
  !         integer(i4b)               :: nRch_local  ! number of reach
  !    
  !         ierr=0; message='accum_water_balance/'
  !    
  !         nRch_local = size(RCHFLX_in)
  !         if (nRch_local/=size(NETOPO_in)) then
  !           ierr=20; message=trim(message)//'RCHFLX and NETOPO has different sizes'; return
  !         end if
  !    
  !         water_budget = 0._dp
  !         do ix = 1, nRch_local
  !           water_budget(1) = water_budget(1) + (RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_VOL(1)-RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_VOL(0))
  !           ! in flux
  !           water_budget(2) = water_budget(2) + RCHFLX_in(ix)%BASIN_QR(1) *dt
  !           if (NETOPO_in(ix)%isLake) then
  !             water_budget(3) = water_budget(3) + RCHFLX_in(ix)%basinprecip *dt
  !           end if
  !           ! out flux
  !           water_budget(4) = water_budget(4) - RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_WM_FLUX_actual *dt
  !           if (NETOPO_in(ix)%isLake) then
  !             water_budget(5) = water_budget(5) - RCHFLX_in(ix)%basinevapo *dt
  !           end if
  !           if (NETOPO_in(ix)%DREACHI==-1 .and. NETOPO_in(ix)%DREACHK<=0) then ! to-do: better way to detect outlet segments
  !             water_budget(6) = water_budget(6) - RCHFLX_in(ix)%ROUTE(ixRoute)%REACH_Q *dt
  !           end if
  !           ! misc. flux information
  !           water_budget(7) = water_budget(7) - RCHFLX_in(ix)%REACH_WM_FLUX *dt
  !         end do
  !       END SUBROUTINE accum_water_balance
  !    
  !     END SUBROUTINE comp_global_wb
  !    
END MODULE water_balance
