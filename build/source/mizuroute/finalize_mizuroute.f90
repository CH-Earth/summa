module finalize_mizuroute_module

  use nr_type, only: i4b

  implicit none
  private

  public :: finalize_mizuroute

contains


  ! **************************************************************************************************
  ! Finalize the coupled mizuRoute model.
  !
  ! Release mizuRoute resources that persist between model evaluations.
  ! Additional cleanup operations can be added here as required.
  ! **************************************************************************************************
  subroutine finalize_mizuroute(err,message)
    
    ! mizuroute global data (shim limiting to data required by summa)
    use globalData, only: rch_routes     ! instantiated routing methods 

    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    character(len=256) :: cmessage

    err = 0
    message = 'finalize_mizuroute/'

    ! deallocate mizuRoute metadata
    call finalize_mizuroute_metadata(err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! deallocate instantiated routing methods
    if(allocated(rch_routes)) deallocate(rch_routes)

    ! additional mizuRoute cleanup can be added here

  end subroutine finalize_mizuroute


  ! **************************************************************************************************
  ! Deallocate allocatable components of persistent mizuRoute metadata structures.
  ! **************************************************************************************************
  subroutine finalize_mizuroute_metadata(err,message)

    ! mizuroute global data (shim, limiting to data required by summa)
    use globalData, only: meta_rflx
    use globalData, only: meta_hflx
    use globalData, only: meta_basinQ
    use globalData, only: meta_irf_bas
    use globalData, only: meta_bas_solute
    use globalData, only: meta_irf
    use globalData, only: meta_kwt
    use globalData, only: meta_kw
    use globalData, only: meta_mc
    use globalData, only: meta_dw
    use globalData, only: meta_solute

    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    err = 0
    message = 'finalize_mizuroute_metadata/'

    call deallocate_varinfo(meta_rflx)
    call deallocate_varinfo(meta_hflx)
    call deallocate_varinfo(meta_basinQ)
    call deallocate_varinfo(meta_irf_bas)
    call deallocate_varinfo(meta_bas_solute)
    call deallocate_varinfo(meta_irf)
    call deallocate_varinfo(meta_kwt)
    call deallocate_varinfo(meta_kw)
    call deallocate_varinfo(meta_mc)
    call deallocate_varinfo(meta_dw)
    call deallocate_varinfo(meta_solute)

  end subroutine finalize_mizuroute_metadata


  ! **************************************************************************************************
  ! Deallocate allocatable components of a mizuRoute variable-metadata array.
  ! **************************************************************************************************
  subroutine deallocate_varinfo(meta)

    use mizuroute_types, only: var_info_new => mizu_var_info_new

    type(var_info_new), intent(inout) :: meta(:)

    integer(i4b) :: iVar

    do iVar=1,size(meta)

      if(allocated(meta(iVar)%varDim)) &
        deallocate(meta(iVar)%varDim)

    enddo

  end subroutine deallocate_varinfo


end module finalize_mizuroute_module
