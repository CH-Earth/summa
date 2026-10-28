! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

! **************************************************************************************************
! Simulated series for calibration targets
!
! A calibration target scores a simulated series against observations.  Routed streamflow comes from
! mizuRoute, but every other target - basin storage change against GRACE, glacier mass balance -
! reads an ordinary SUMMA variable that until now only reached the output file.
!
! This module collects any such variable, by name, as a time series over one spatial unit of the
! domain: the name is resolved once against SUMMA's own metadata, and the value is then reduced and
! recorded each model time step.  Resolving against the metadata rather than a list kept here means a
! target can name any SUMMA variable, and names stay honest - an unknown one is refused when the
! calibration starts rather than silently scored against something else.
!
! The spatial unit is resolved the same way, and for the same reason.  A domain mean is right for a
! single basin and wrong for anything else: two basins in one domain hold their own storage, and an
! id the domain does not contain is refused at start-up rather than scored against a neighbour.
! **************************************************************************************************
module simulated_series

USE nr_type,    only: i4b,i8b,rkind,lgt
USE summa_type, only: summa1_type_dec

USE globalData, only: realMissing
USE globalData, only: integerMissing
USE globalData, only: gru_struc              ! gru-hru mapping, with the domain information

USE var_lookup, only: iLookATTR              ! named variables for the local attributes
USE var_lookup, only: iLookBVAR              ! named variables for the basin-average variables
USE var_lookup, only: iLookPROG              ! named variables for the prognostic variables

USE build_options, only: mizuroute_active

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_coupling, only: get_mizuroute_streamflow
#endif

  implicit none
  private

  ! which SUMMA structure a simulated series is read from
  integer(i4b), parameter, public :: ix_series_bvar       = 1   ! basin-average variables (GRU)
  integer(i4b), parameter, public :: ix_series_diag       = 2   ! diagnostic variables   (HRU/domain)
  integer(i4b), parameter, public :: ix_series_prog       = 3   ! prognostic variables   (HRU/domain)
  integer(i4b), parameter, public :: ix_series_flux       = 4   ! model fluxes           (HRU/domain)
  integer(i4b), parameter, public :: ix_series_param      = 5   ! model parameters       (HRU/domain)
  integer(i4b), parameter, public :: ix_series_streamflow = 6   ! simulated streamflow   (GRU/reach)

  ! which spatial unit of the domain a series is reduced over
  integer(i4b), parameter, public :: ix_unit_domain = 1    ! every GRU, area-weighted
  integer(i4b), parameter, public :: ix_unit_gru    = 2    ! one GRU, its HRUs area-weighted
  integer(i4b), parameter, public :: ix_unit_hru    = 3    ! one HRU
  integer(i4b), parameter, public :: ix_unit_reach  = 4    ! one river reach, routed by mizuRoute

  ! One simulated series: a named SUMMA variable, resolved to its structure and index, over one
  ! spatial unit of the domain, and the value of that variable there at every model time step.
  type, public :: sim_series_type
    character(len=64)        :: name = ''                 ! variable name, as the target asked for it
    character(len=64)        :: units = ''                ! its units, from SUMMA's own metadata
    integer(i4b)             :: ix_struct = integerMissing ! structure the variable lives in
    integer(i4b)             :: ix_var    = integerMissing ! index of the variable within it
    integer(i4b)             :: ix_unit   = ix_unit_domain ! spatial unit the value is reduced over
    integer(i8b)             :: unit_id   = integerMissing ! id of that unit, as the target named it
    integer(i4b)             :: ix_gru    = integerMissing ! local index of the selected GRU
    integer(i4b)             :: ix_hru    = integerMissing ! index of the selected HRU within it
    integer(i4b)             :: ix_seg    = integerMissing ! index of the selected reach
    logical(lgt)             :: volumetric = .false.       ! a rate per unit area summed to a volume rate over the unit
    real(rkind), allocatable :: values(:)                  ! its value at each model time step
  end type sim_series_type

  public :: init_simulated_series
  public :: collect_simulated_series
  public :: find_simulated_series
  public :: get_outlet_streamflow
  public :: is_routed_streamflow
  public :: spatial_unit_index

contains

  ! **************************************************************************************************
  ! Report whether a variable name means the routed streamflow mizuRoute produces.
  !
  ! Streamflow is the one simulated series that is not a SUMMA variable, so it is named rather than
  ! resolved, and every caller asks the same way.
  ! **************************************************************************************************
  pure function is_routed_streamflow(name) result(isFlow)
    character(*), intent(in) :: name
    logical(lgt)             :: isFlow

    select case(trim(name))
      case ('streamflow','discharge'); isFlow = .true.
      case default;                    isFlow = .false.
    end select

  end function is_routed_streamflow

  ! **************************************************************************************************
  ! Get streamflow at the basin outlet.
  ! Use mizuRoute when active; otherwise sum routed runoff from all GRUs.
  ! **************************************************************************************************
  subroutine get_outlet_streamflow(modelTimeStep, summa_struct, flowSim)

    integer(i4b),          intent(in)  :: modelTimeStep
    type(summa1_type_dec), intent(in)  :: summa_struct
    real(rkind),           intent(out) :: flowSim

    integer(i4b) :: iGRU

    ! Get routed streamflow from mizuRoute when available and enabled
    if(mizuroute_active)then
      if(summa_struct%config%use_mizuroute)then
        call get_mizuroute_streamflow(modelTimeStep,summa_struct,flowSim)
        return
      endif
    endif

    ! Otherwise calculate streamflow from SUMMA routed runoff
    flowSim = 0._rkind

    do iGRU=1,size(summa_struct%bvarStruct%gru)
      flowSim = flowSim + &
        summa_struct%bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1) * &
        summa_struct%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
    enddo

  end subroutine get_outlet_streamflow

  ! **************************************************************************************************
  ! The spatial unit a configuration word names, or integerMissing when it names none of them.
  ! **************************************************************************************************
  pure function spatial_unit_index(name) result(ix)
    character(*), intent(in) :: name
    integer(i4b)             :: ix

    select case(trim(name))
      case ('domain'); ix = ix_unit_domain
      case ('gru');    ix = ix_unit_gru
      case ('hru');    ix = ix_unit_hru
      case ('reach');  ix = ix_unit_reach
      case default;    ix = integerMissing
    end select

  end function spatial_unit_index

  ! **************************************************************************************************
  ! Resolve a list of variable names and spatial units, and allocate a series for each.
  !
  ! Each name is looked up in SUMMA's metadata, in turn, as a basin-average variable, a diagnostic, a
  ! prognostic, a flux, and a parameter.  A name that matches none of those is refused here, at
  ! start-up, naming the variable that could not be found.
  !
  ! The spatial unit is resolved the same way, from the id the target named to the index the run uses,
  ! so an id the domain does not contain is refused here rather than silently scored against another.
  ! A reach is the exception: it belongs to the river network rather than to SUMMA, so it is left for
  ! the caller to resolve against the mizuRoute topology.
  !
  ! Streamflow at the basin outlet is collected using get_outlet_streamflow.
  ! The existing spatial-unit identifier ix_unit_domain is used for this target.
  ! Streamflow from a GRU without mizuRoute is calculated from its routed runoff and area.
  ! **************************************************************************************************
  subroutine init_simulated_series(names,unitKind,unitId,numtim,series,err,message)
    USE get_ixname_module, only: get_ixBvar,get_ixDiag,get_ixProg,get_ixFlux,get_ixParam
    USE globalData, only: bvar_meta,diag_meta,prog_meta,flux_meta,mpar_meta
    implicit none
    character(*),                           intent(in)  :: names(:)    ! variable names to collect
    integer(i4b),                           intent(in)  :: unitKind(:) ! spatial unit of each of them
    integer(i8b),                           intent(in)  :: unitId(:)   ! id of that unit
    integer(i4b),                           intent(in)  :: numtim      ! number of model time steps
    type(sim_series_type), allocatable,     intent(out) :: series(:)   ! the resolved series
    integer(i4b),                           intent(out) :: err         ! error code
    character(*),                           intent(out) :: message     ! error message
    integer(i4b) :: iSeries
    integer(i4b) :: ix
    character(len=256) :: cmessage

    err=0
    message='init_simulated_series/'

    allocate(series(size(names)),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate the simulated series'
      return
    endif

    do iSeries=1,size(names)
      series(iSeries)%name=trim(names(iSeries))
      series(iSeries)%ix_unit=unitKind(iSeries)
      series(iSeries)%unit_id=unitId(iSeries)

      ! the spatial unit, from the id the target named to the index the run uses
      call resolve_spatial_unit(series(iSeries),err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! Streamflow at a river reach
      if(series(iSeries)%ix_unit == ix_unit_reach)then
        if(is_routed_streamflow(names(iSeries)))then
          series(iSeries)%ix_struct=ix_series_streamflow
          series(iSeries)%units='m3 s-1'
        else
          message=trim(message)//'"'//trim(names(iSeries))//'" is asked for on a reach, which carries '// &
                  'routed streamflow only'
          err=20; return
        endif
      
      ! Streamflow at the basin outlet
      else if(is_routed_streamflow(names(iSeries)) .and. series(iSeries)%ix_unit == ix_unit_domain)then
        series(iSeries)%ix_struct=ix_series_streamflow
        series(iSeries)%units='m3 s-1'
      
      ! Streamflow from a GRU without mizuRoute
      else if(is_routed_streamflow(names(iSeries)) .and. series(iSeries)%ix_unit == ix_unit_gru)then
        series(iSeries)%ix_struct=ix_series_bvar
        series(iSeries)%ix_var=iLookBVAR%averageRoutedRunoff
        series(iSeries)%volumetric=.true.
        series(iSeries)%units='m3 s-1'
      
      else
        ! resolve the name against SUMMA's metadata, most likely structure first
        ix=get_ixBvar(trim(names(iSeries)))
        if(ix/=integerMissing)then
          series(iSeries)%ix_struct=ix_series_bvar
        else
          ix=get_ixDiag(trim(names(iSeries)))
          if(ix/=integerMissing)then
            series(iSeries)%ix_struct=ix_series_diag
          else
            ix=get_ixProg(trim(names(iSeries)))
            if(ix/=integerMissing)then
              series(iSeries)%ix_struct=ix_series_prog
            else
              ix=get_ixFlux(trim(names(iSeries)))
              if(ix/=integerMissing)then
                series(iSeries)%ix_struct=ix_series_flux
              else
                ix=get_ixParam(trim(names(iSeries)))
                if(ix/=integerMissing) series(iSeries)%ix_struct=ix_series_param
              endif
            endif
          endif
        endif

        if(ix==integerMissing)then
          message=trim(message)//'"'//trim(names(iSeries))//'" is not a SUMMA variable: a calibration '// &
                  'target must name "streamflow" or a variable SUMMA knows'
          err=20; return
        endif
        series(iSeries)%ix_var=ix

        ! a basin variable is held once per GRU, so there is no one HRU's value to take
        if(series(iSeries)%ix_struct == ix_series_bvar .and. series(iSeries)%ix_unit == ix_unit_hru)then
          message=trim(message)//'"'//trim(names(iSeries))//'" is held per GRU, not per HRU: score it over a gru'
          err=20; return
        endif

        ! units come from the same metadata as the variable, so the check against the observation
        ! file's units means something rather than being waved through
        select case(series(iSeries)%ix_struct)
          case (ix_series_bvar);  series(iSeries)%units=trim(bvar_meta(ix)%varunit)
          case (ix_series_diag);  series(iSeries)%units=trim(diag_meta(ix)%varunit)
          case (ix_series_prog);  series(iSeries)%units=trim(prog_meta(ix)%varunit)
          case (ix_series_flux);  series(iSeries)%units=trim(flux_meta(ix)%varunit)
          case (ix_series_param); series(iSeries)%units=trim(mpar_meta(ix)%varunit)
        end select
      endif

      allocate(series(iSeries)%values(numtim),source=realMissing,stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate values for simulated series "'//trim(names(iSeries))//'"'
        return
      endif
    enddo

  end subroutine init_simulated_series

  ! **************************************************************************************************
  ! Resolve the id of a series' spatial unit to the index the run holds it at.
  !
  ! GRU and HRU ids are the ones the domain numbers its units by, not positions in this run, so a
  ! subset run finds its units wherever they happen to sit.  A reach is left alone: the river network
  ! is mizuRoute's, and only the caller holding the topology can resolve it.
  ! **************************************************************************************************
  subroutine resolve_spatial_unit(series,err,message)
    implicit none
    type(sim_series_type), intent(inout) :: series
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message
    integer(i4b) :: iGRU,jHRU

    err=0
    message='resolve_spatial_unit/'

    select case(series%ix_unit)

      case (ix_unit_domain, ix_unit_reach)
        ! the whole domain needs no index, and a reach is the caller's to resolve

      case (ix_unit_gru)
        do iGRU=1,size(gru_struc)
          if(gru_struc(iGRU)%gru_id == series%unit_id)then
            series%ix_gru=iGRU
            exit
          endif
        enddo

      case (ix_unit_hru)
        outer: do iGRU=1,size(gru_struc)
          do jHRU=1,gru_struc(iGRU)%hruCount
            if(gru_struc(iGRU)%hruInfo(jHRU)%hru_id == series%unit_id)then
              series%ix_gru=iGRU
              series%ix_hru=jHRU
              exit outer
            endif
          enddo
        enddo outer

      case default
        message=trim(message)//'simulated series "'//trim(series%name)//'" asks for an unknown spatial unit'
        err=20; return

    end select

    if(series%ix_unit == ix_unit_gru .or. series%ix_unit == ix_unit_hru)then
      if(series%ix_gru == integerMissing)then
        write(message,'(a,i0,a)') trim(message)//'simulated series "'//trim(series%name)//'" asks for '// &
              merge('gru','hru',series%ix_unit == ix_unit_gru)//' ',series%unit_id,', which this domain does not hold'
        err=20; return
      endif
    endif

  end subroutine resolve_spatial_unit

  ! **************************************************************************************************
  ! Return the index of the series of a variable over a spatial unit, or integerMissing when it was not
  ! collected.  The unit is part of the key: one variable over two basins is two different series.
  ! **************************************************************************************************
  pure function find_simulated_series(series,name,unitKind,unitId) result(iSeries)
    type(sim_series_type), intent(in) :: series(:)
    character(*),          intent(in) :: name
    integer(i4b),          intent(in) :: unitKind
    integer(i8b),          intent(in) :: unitId
    integer(i4b)                      :: iSeries
    integer(i4b) :: i

    iSeries=integerMissing
    do i=1,size(series)
      if(trim(series(i)%name)/=trim(name)) cycle
      if(series(i)%ix_unit/=unitKind) cycle
      if(unitKind/=ix_unit_domain .and. series(i)%unit_id/=unitId) cycle
      iSeries=i
      return
    enddo

  end function find_simulated_series

  ! **************************************************************************************************
  ! Record the value of every collected series over its spatial unit for one model time step.
  !
  ! Basin-average variables are one value per GRU; the rest are per HRU and domain.  Both are reduced
  ! to one number the same way, as an area-weighted mean over the GRUs and HRUs the unit covers - every
  ! one for the domain, one GRU's for a GRU, one HRU for an HRU - so a target compares a series against
  ! an observation of the same place whichever structure its variable came from.  Storage and flux
  ! densities (per unit area) are what this suits.  A reach belongs to the network, and the caller records it.
  ! **************************************************************************************************
  subroutine collect_simulated_series(modelTimeStep,summa_struct,series,err,message)
    implicit none
    integer(i4b),          intent(in)    :: modelTimeStep  ! index of the model time step
    type(summa1_type_dec), intent(in)    :: summa_struct   ! top-level SUMMA data structure
    type(sim_series_type), intent(inout) :: series(:)      ! series to record into
    integer(i4b),          intent(out)   :: err            ! error code
    character(*),          intent(out)   :: message        ! error message
    integer(i4b) :: iSeries
    integer(i4b) :: iGRU,jHRU,iDOM
    integer(i4b) :: gruFirst,gruLast                       ! GRUs the unit covers
    integer(i4b) :: hruFirst,hruLast                       ! HRUs it covers within each of them
    real(rkind)  :: total                                  ! area-weighted sum over the unit
    real(rkind)  :: area                                   ! total weighting area
    real(rkind)  :: wt                                     ! weight of one contribution
    real(rkind)  :: domTotal                               ! total domain area of an HRU (m2)
    real(rkind)  :: domWt                                  ! share of one domain in its HRU (-)

    err=0
    message='collect_simulated_series/'
    if(modelTimeStep < 1) return

    do iSeries=1,size(series)
      if(modelTimeStep > size(series(iSeries)%values)) cycle

      ! Collect streamflow at the basin outlet
      if(series(iSeries)%ix_struct == ix_series_streamflow .and. &
         series(iSeries)%ix_unit == ix_unit_domain)then

        call get_outlet_streamflow(modelTimeStep, summa_struct, &
                                  series(iSeries)%values(modelTimeStep))
        cycle

      endif

      ! Collect streamflow at a requested river reach
      if(series(iSeries)%ix_unit == ix_unit_reach)then

        if(mizuroute_active)then
          if(summa_struct%config%use_mizuroute)then
            call get_mizuroute_streamflow(modelTimeStep, summa_struct, &
                                         series(iSeries)%values(modelTimeStep), &
                                         series(iSeries)%ix_seg)
          endif
        endif

        cycle

      endif

      ! Collect ordinary SUMMA variables
      total=0._rkind
      area =0._rkind

      if(series(iSeries)%ix_unit == ix_unit_domain)then
        gruFirst=1
        gruLast =summa_struct%nGRU_local
      else
        gruFirst=series(iSeries)%ix_gru
        gruLast =series(iSeries)%ix_gru
      endif

      select case(series(iSeries)%ix_struct)

        ! basin-average variables: one value per GRU, weighted by basin area
        case (ix_series_bvar)
          do iGRU=gruFirst,gruLast
            wt=summa_struct%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
            if(wt <= 0._rkind) wt=1._rkind
            total=total+wt*summa_struct%bvarStruct%gru(iGRU)%var(series(iSeries)%ix_var)%dat(1)
            area =area +wt
          enddo

        ! everything else is per HRU and domain, weighted by HRU area and, within an HRU, by domain area
        case default
          do iGRU=gruFirst,gruLast
            if(series(iSeries)%ix_unit == ix_unit_hru)then
              hruFirst=series(iSeries)%ix_hru
              hruLast =series(iSeries)%ix_hru
            else
              hruFirst=1
              hruLast =gru_struc(iGRU)%hruCount
            endif
            do jHRU=hruFirst,hruLast
              wt=summa_struct%attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
              if(wt <= 0._rkind) wt=1._rkind
              domTotal=0._rkind
              do iDOM=1,gru_struc(iGRU)%hruInfo(jHRU)%domCount
                domTotal=domTotal+max(summa_struct%progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1),0._rkind)
              enddo
              do iDOM=1,gru_struc(iGRU)%hruInfo(jHRU)%domCount
                ! a domain without area holds no state worth counting, unless no domain has area
                if(domTotal > 0._rkind)then
                  domWt=max(summa_struct%progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1),0._rkind)/domTotal
                else
                  domWt=1._rkind/real(gru_struc(iGRU)%hruInfo(jHRU)%domCount,rkind)
                endif
                if(domWt <= 0._rkind) cycle
                total=total+wt*domWt*series_value(summa_struct,series(iSeries),iGRU,jHRU,iDOM)
                area =area +wt*domWt
              enddo
            enddo
          enddo

      end select

      if(series(iSeries)%volumetric)then
        series(iSeries)%values(modelTimeStep)=total
      else if(area > 0._rkind)then
        series(iSeries)%values(modelTimeStep)=total/area
      else
        series(iSeries)%values(modelTimeStep)=realMissing
      endif
    enddo

  end subroutine collect_simulated_series

  ! **************************************************************************************************
  ! Read one HRU/domain value of a series from the structure it lives in.
  ! **************************************************************************************************
  function series_value(summa_struct,series,iGRU,jHRU,iDOM) result(value)
    implicit none
    type(summa1_type_dec), intent(in) :: summa_struct
    type(sim_series_type), intent(in) :: series
    integer(i4b),          intent(in) :: iGRU,jHRU,iDOM
    real(rkind)                       :: value

    select case(series%ix_struct)
      case (ix_series_diag)
        value=summa_struct%diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case (ix_series_prog)
        value=summa_struct%progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case (ix_series_flux)
        value=summa_struct%fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case (ix_series_param)
        value=summa_struct%mparStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case default
        value=realMissing
    end select

  end function series_value

end module simulated_series
