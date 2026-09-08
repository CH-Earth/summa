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
  public :: get_mizuroute_streamflow

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

  ! *****************************************************************************
  ! Spatial organization of runoff passed from SUMMA to mizuRoute
  ! *****************************************************************************
  !
  ! SUMMA uses a hierarchical spatial structure in which one or more HRUs are
  ! contained within each GRU:
  !
  !       SUMMA HRU_1 ----\
  !       SUMMA HRU_2 -----+--> SUMMA GRU
  !             ...       /
  !       SUMMA HRU_n ----/
  !
  ! SUMMA can optionally simulate lateral flow among HRUs within a GRU. Runoff
  ! from the HRUs is then aggregated to the GRU level. SUMMA subsequently
  ! applies its runoff-routing calculation, which represents the time delay
  ! associated with routing through the unresolved river network.
  !
  ! The current coupling therefore supplies one routed runoff value per
  ! SUMMA GRU:
  !
  !       SUMMA HRUs
  !            |
  !            | optional lateral flow among HRUs
  !            v
  !       aggregate runoff to SUMMA GRU
  !            |
  !            | SUMMA time-delay routing through
  !            | the unresolved river network
  !            v
  !       SUMMA GRU runoff
  !            |
  !            | ------------------
  !            | coupling interface
  !            | ------------------
  !            v
  !       mizuRoute host-model runoff
  !            |
  !            | optional spatial remapping
  !            v
  !       mizuRoute river-network HRU runoff
  !            |
  !            | aggregate runoff to reaches
  !            v
  !       lateral inflow to river reaches
  !            |
  !            | route through explicit river network
  !            v
  !       routed streamflow
  !
  !
  ! Coupling interface
  ! ------------------
  !
  ! The coupling interface contains the runoff values and corresponding SUMMA
  ! GRU identifiers:
  !
  !       coupling(:)%id      SUMMA GRU IDs
  !       coupling(:)%qsim    routed SUMMA GRU runoff [m s-1]
  !
  ! These are transferred to the native mizuRoute runoff structure:
  !
  !       runoff%hru_id(:)    IDs of the host-model runoff elements
  !       runoff%sim(:)       runoff on those elements [m s-1]
  !
  ! Thus:
  !
  !       coupling(:)%id   --> runoff%hru_id(:)
  !       coupling(:)%qsim --> runoff%sim(:)
  !
  ! The name runoff%hru_id can be confusing in the coupled configuration.
  ! It does not imply that these IDs are the mizuRoute river-network HRU IDs.
  ! Rather, runoff%hru_id identifies the spatial elements on which runoff is
  ! supplied by the host land model. For the current SUMMA coupling these
  ! elements are GRUs, so runoff%hru_id contains SUMMA GRU IDs.
  !
  !
  ! Spatial remapping
  ! -----------------
  !
  ! mizuRoute may receive runoff on a spatial domain that differs from the
  ! HRUs associated with its river network. When this occurs, runoff is
  ! spatially remapped before river-network routing:
  !
  !       runoff%hru_id / runoff%sim
  !       host-model runoff elements
  !                    |
  !                    | spatial remapping
  !                    v
  !       runoff%basinRunoff
  !       runoff on river-network HRUs
  !
  ! Within the remapping machinery, the host-model runoff elements are called
  ! qHRUs ("runoff HRUs") to distinguish them from the destination
  ! river-network HRUs. The remapping structure therefore contains:
  !
  !       remap%qhru_id(:)    IDs of source runoff elements (qHRUs)
  !       remap%hru_id(:)     IDs of destination river-network HRUs
  !
  ! The qHRU terminology is specific to the remapping relationship.
  ! remap%qhru_id refers back to the host-model IDs stored in runoff%hru_id.
  ! The IDs in remap%qhru_id are matched against runoff%hru_id to locate the
  ! corresponding values in runoff%sim. Those values are then remapped onto
  ! the river-network HRUs identified by remap%hru_id.
  !
  ! For example, if SUMMA supplies runoff from a single GRU with ID 101 that
  ! contributes to multiple mizuRoute river-network HRUs:
  !
  !       runoff%hru_id       = [101]
  !       runoff%sim          = [runoff from SUMMA GRU 101]
  !
  !       remap%qhru_id       = [101, 101, 101, ...]
  !       remap%hru_id        = [river-network HRU IDs ...]
  !
  ! Each occurrence of qhru_id=101 therefore refers to the same source runoff
  ! value stored in runoff%sim(1). The remapping weights distribute that runoff
  ! onto the corresponding river-network HRUs, producing runoff%basinRunoff.
  !
  ! Spatial remapping is optional. If the host land model already supplies
  ! runoff on the mizuRoute river-network HRUs, the remapping step is skipped
  ! and runoff%sim is used directly as runoff%basinRunoff.
  !
  !
  ! River-network routing
  ! ---------------------
  !
  ! After the optional spatial-remapping step, the coupled workflow is:
  !
  !       runoff%sim
  !       runoff supplied by SUMMA
  !              |
  !              | optional spatial remapping
  !              v
  !       runoff%basinRunoff
  !       runoff on river-network HRUs
  !              |
  !              | basin2reach:
  !              | aggregate river-network HRUs to reaches and
  !              | convert runoff depth to lateral reach inflow
  !              v
  !       lateral inflow to river reaches
  !              |
  !              | route_network
  !              v
  !       routed streamflow
  !
  ! Spatial remapping and river-network routing are therefore distinct
  ! operations. Spatial remapping only reconciles the spatial discretization
  ! of runoff supplied by the host land model with the HRUs associated with
  ! the mizuRoute river network. It can be omitted when those spatial
  ! discretizations already coincide.
  !
  ! The coupled implementation does not include the mizuRoute runoff-routing
  ! routine for the unresolved river network. That component is deliberately
  ! excluded from the set of mizuRoute source files compiled and linked into
  ! SUMMA, because the corresponding time-delay routing is always performed by
  ! SUMMA before runoff crosses the coupling interface.
  !
  ! In the coupled configuration, mizuRoute therefore begins with the runoff
  ! supplied by SUMMA and is responsible only for spatial remapping (when
  ! required), aggregation of runoff to river reaches, and routing through the
  ! explicit river network.
  !
  ! Possible HRU-level coupling
  ! ---------------------------
  !
  ! Although the current implementation passes runoff at SUMMA GRU resolution,
  ! the coupling interface is not fundamentally restricted to GRUs. SUMMA
  ! runoff could instead be supplied directly at HRU resolution, with SUMMA
  ! HRU IDs stored in runoff%hru_id. In that configuration the SUMMA HRUs would
  ! be the source qHRUs and mizuRoute could perform the required spatial
  ! aggregation and routing through the unresolved river network. HRU-level
  ! coupling is not currently implemented.
  !
  ! *****************************************************************************

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------
  ! Initialize mizuRoute within the SUMMA data structures
  !-----------------------------------------------------------------------
  subroutine init_mizuroute_from_summa(summaStruct, ierr, message)

  USE public_var,     only: iulog

  USE nr_utils,       only: match_index
  USE init_mizuRoute, only: init_mizuroute_domain 

  type(summa1_type_dec), intent(inout) :: summaStruct
  integer,               intent(out)   :: ierr
  character(*),          intent(out)   :: message
  
  real(rkind)                          :: length_conv
  real(rkind)                          :: time_conv
  integer(i4b)                         :: iGRU
  integer(i4b)                         :: nSpace(1:2) = integerMissing
  integer(i4b)                         :: n_write
  character(len=256)                   :: cmessage
  
  ierr = 0
  message = 'init_mizuroute_from_summa/'
  
  associate(info     => summaStruct%config%mizu_info,   &
            domain   => summaStruct%mizu_domain  )

  ! -----------------------------------------------------------------------
  ! Define host-model information required by mizuRoute
  ! -----------------------------------------------------------------------
 
  ! general info
  info%is_print     = .true.
  info%do_mizuroute = .true.
  info%do_remapping = allocated(info%remap%remap_file)
 
  ! logging
  iulog = summaStruct%config%iulog_summa

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
 
  ! -----------------------------------------------------------------------
  ! Initialize the mizuRoute domain
  !
  ! This performs the major mizuRoute initialization operations, including:
  !   - reading routing configuration and parameter information
  !   - reading the river-network topology
  !   - constructing the river-network data structures
  !   - allocating the runoff and river-routing data structures
  !   - populating the host-model runoff IDs
  !   - reading spatial-remapping information, when required
  !   - constructing the indices required for spatial remapping
  ! -----------------------------------------------------------------------
  
  call init_mizuroute_domain(summaStruct%instance_parallel%rank, &
                             info, domain, nSpace, n_write,      &
                             summaStruct%coupling(:)%id,         &
                             length_conv, time_conv,             &
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
  
  associate(info   => summaStruct%config%mizu_info,   &
            domain => summaStruct%mizu_domain)
  
    ! Determine the index of the output buffer (if writePerStep n_write=1)
    idx_buff = merge(1, modelTimeStep, summaStruct%n_write == 1)
  
    ! Transfer summa routed runoff into the mizuRoute runoff structure
    domain%river_network%core%runoff%sim(:) = summaStruct%coupling(:)%qsim

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
  
  call define_mizuroute_output(ncid,                          &
                               summaStruct%config%mizu_info,  &
                               summaStruct%mizu_domain,       &
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

  call write_mizuroute_output(ncid,                         &
                              istart,                       &
                              numtim,                       &
                              summaStruct%config%mizu_info, &
                              summaStruct%mizu_domain,      &
                              ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine write_mizuroute_output_from_summa

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Get mizuRoute streamflow
  !-----------------------------------------------------------------------
  subroutine get_mizuroute_streamflow(modelTimeStep, summaStruct, simFlow)

  integer(i4b),          intent(in)  :: modelTimeStep
  type(summa1_type_dec), intent(in)  :: summaStruct
  real(rkind),           intent(out) :: simFlow

  integer(i4b) :: idx_buff
  integer(i4b) :: ixSeg

  idx_buff = merge(1, modelTimeStep, summaStruct%n_write == 1)
  ixSeg    = summaStruct%config%mizu_info%ntopo%ixSegOut

  simFlow = &
    summaStruct%mizu_domain%river_network%driver%method(1)%streamflow(ixSeg,idx_buff)

  end subroutine get_mizuroute_streamflow

end module mizuroute_coupling
