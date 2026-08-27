module mizuroute_coupling

  USE nr_type, only: i4b, rkind
  USE summa_type, only:summa1_type_dec

  ! mizuRoute public interface
  use public_var, only: integerMissing

  implicit none
  private

  public :: init_mizuroute_from_summa
  public :: route_mizuroute_from_summa
  public :: define_mizuroute_output_from_summa
  public :: write_mizuroute_output_from_summa

  ! *****************************************************************************
  ! SUMMA--mizuRoute coupling interface
  ! *****************************************************************************
  !
  ! This module provides the thin interface between SUMMA and mizuRoute.
  ! It translates between the SUMMA host-model interface and explicitly named
  ! mizuRoute interfaces, but should not contain implementation logic belonging
  ! to either model.
  !
  ! IMPORTANT:
  ! Do not USE modules with names shared by SUMMA and mizuRoute in this module
  ! (e.g., globalData, var_lookup, read_param_module, popMetadat_module).
  ! Both models define modules with these names, and exposing them in the
  ! coupling layer can result in ambiguous or incorrect module resolution.
  !
  ! Native SUMMA logic should remain in SUMMA routines, and native mizuRoute
  ! logic should remain in routines compiled as part of the MIZUROUTE target.
  ! This coupling module should only:
  !
  !   - access the SUMMA structure exposed at the coupling interface;
  !   - transfer data between SUMMA and mizuRoute structures; and
  !   - call uniquely named mizuRoute initialization, routing, and I/O routines.
  !
  ! Keeping this interface thin isolates the SUMMA and mizuRoute module
  ! namespaces and avoids the need to rename shared module names in either model.
  ! *****************************************************************************

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------
  ! Initialize mizuRoute within the SUMMA data structures
  !-----------------------------------------------------------------------
  subroutine init_mizuroute_from_summa(summaStruct, ierr, message)
 
  USE init_mizuRoute, only: init_mizuroute_domain 

  type(summa1_type_dec), intent(inout) :: summaStruct
  integer,               intent(out)   :: ierr
  character(*),          intent(out)   :: message
  
  real(rkind)                          :: length_conv
  real(rkind)                          :: time_conv
  integer(i4b)                         :: nSpace(1:2) = integerMissing
  integer(i4b)                         :: n_write
  character(len=256)                   :: cmessage
  
  ierr = 0
  message = 'init_mizuroute_from_summa/'
  
  associate(info   => summaStruct%mizu_info,   &
            domain => summaStruct%mizu_domain)
  
  ! ---- transfer information from summa ----
  
  ! general info
  info%is_print     = .true.
  info%do_mizuroute = .true.
  info%do_remapping = allocated(info%remap%remap_file)
  
  ! time information
  n_write           = summaStruct%n_write
  info%dt_landmodel = summaStruct%data_step
  
  ! SUMMA provides runoff on a one-dimensional basin (GRU) domain 
  nSpace(1) = summaStruct%nGRU_local
  nSpace(2) = integerMissing
  
  info%is_gridded = (nSpace(2) /= integerMissing)
  
  ! ---- initialize unit conversions (multipliers) ----
  length_conv = 1._rkind   ! no conversion needed: summa runoff length = m
  time_conv   = 1._rkind   ! no conversion needed: summa runoff time = s-1
  
  ! ---- general routine that can work with all host land models
  call init_mizuroute_domain(info, domain, nSpace, n_write, length_conv, time_conv, &
                             ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end associate
  
  end subroutine init_mizuroute_from_summa
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  ! Network routing in mizuRoute
  !-----------------------------------------------------------------------
  subroutine route_mizuroute_from_summa(modelTimeStep, summaStruct, ierr, message)
 
  USE network_routing_module, only: network_routing

  integer(i4b),          intent(in)    :: modelTimeStep
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message
  
  integer(i4b)       :: iGRU
  integer(i4b)       :: idx_buff
  character(len=256) :: cmessage
  
  ierr    = 0
  message = 'route_mizuroute_from_summa/'
  
  associate(info   => summaStruct%mizu_info,   &
            domain => summaStruct%mizu_domain)
  
    ! Determine the index of the output buffer (if writePerStep n_write=1)
    idx_buff = merge(1, modelTimeStep, summaStruct%n_write == 1)
  
    ! Transfer summa routed runoff into the mizuRoute runoff structure
    domain%river_network%runoff%sim(:) = summaStruct%routedRunoff(:)
  
    ! Route the complete runoff field
    call network_routing(idx_buff,             &
                         domain%river_network, &
                         domain%remap%routing, &
                         info%do_remapping,    &
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end associate
  
  end subroutine route_mizuroute_from_summa
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------  

  !-----------------------------------------------------------------------
  ! Define mizuRoute output based on the SUMMA model structure
  !-----------------------------------------------------------------------
  subroutine define_mizuroute_output_from_summa(ncid, summaStruct, ierr, message)
  
  USE mizuroute_output_module, only: define_mizuroute_output

  integer(i4b),          intent(in)  :: ncid
  type(summa1_type_dec), intent(in)  :: summaStruct
  integer(i4b),          intent(out) :: ierr
  character(*),          intent(out) :: message
  
  character(len=256) :: cmessage
  
  ierr = 0
  message = 'define_mizuroute_output_from_summa/'
  
  call define_mizuroute_output(ncid,                      &
                               summaStruct%mizu_info,     &
                               summaStruct%mizu_domain,   &
                               ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine define_mizuroute_output_from_summa

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Write mizuRoute output from the SUMMA model structure
  !-----------------------------------------------------------------------
  subroutine write_mizuroute_output_from_summa(ncid, istart, numtim, &
                                               summaStruct, ierr, message)

  USE mizuroute_output_module, only: write_mizuroute_output

  integer(i4b),          intent(in)    :: ncid
  integer(i4b),          intent(in)    :: istart
  integer(i4b),          intent(in)    :: numtim
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message

  character(len=256) :: cmessage

  ierr = 0
  message = 'write_mizuroute_output_from_summa/'

  call write_mizuroute_output(ncid,                    &
                              istart,                   &
                              numtim,                   &
                              summaStruct%mizu_info,    &
                              summaStruct%mizu_domain,  &
                              ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

end subroutine write_mizuroute_output_from_summa


end module mizuroute_coupling
