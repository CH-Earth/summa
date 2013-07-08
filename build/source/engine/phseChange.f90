module phseChange_module
USE nrtype
! physical constants
USE multiconst,only:&
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water     ! intrinsic density of water    (kg m-3)
implicit none
private
public::phseChange
contains

 ! ************************************************************************************************
 ! new subroutine: compute phase change impacts on matric head and volumetric liquid water and ice
 ! ************************************************************************************************
 subroutine phseChange(mLayerTempNew,       & ! intent(in): new temperature vector (K)
                       mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                       mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                       mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                       mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                       mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                       mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                       err,message)           ! intent(out): error control
 ! utility routines
 USE snow_utils_module,only:fracliquid    ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:volFracLiq    ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead    ! compute the matric head based on volumetric liquid water content
 ! data structures
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! input variables
 real(dp),intent(in)           :: mLayerTempNew(:)         ! new estimate of temperature (K)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! before phase change: matric head (m)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! before phase change: volumetric fraction of liquid water (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! before phase change: volumetric fraction of ice (-)
 ! output variables
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! after phase change: matric head (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! after phase change: volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! after phase change: volumetric fraction of ice (-)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
 real(dp),pointer              :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),pointer              :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),pointer              :: theta_sat                ! soil porosity (-)
 real(dp),pointer              :: theta_res                ! soil residual volumetric water content (-)
 real(dp),pointer              :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),pointer              :: kappa                    ! constant in the freezing curve function (m K-1)
 ! local pointers to model variables
 real(dp),pointer              :: mLayerVolFracIce(:)      ! volumetric fraction of ice at the start of the iteration (-)
 real(dp),pointer              :: mLayerTcrit(:)           ! critical soil temperature above which all water is unfrozen (K)
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 real(dp)                      :: fLiq                     ! fraction of liquid water (-)
 real(dp)                      :: dLiq                     ! change in volumetric liiquid water over the iteration (-)
 real(dp)                      :: theta                    ! liquid water equivalent of total water (-)
 integer(i4b)                  :: nSnow                    ! number of snow layers
 integer(i4b)                  :: iLayer                   ! index of model layer
 logical(lgt)                  :: printflag                ! flag to print debug information
 ! initialize error control
 err=0; message="phsechange/"

 ! initialize print flag
 printflag=.false.

 ! assign pointers to model parameters
 snowfrz_scale    => mpar_data%var(iLookPARAM%snowfrz_scale)        ! scaling parameter for the snow freezing curve (K-1)
 vGn_alpha        => mpar_data%var(iLookPARAM%vGn_alpha)            ! van Genutchen "alpha" parameter (m-1)
 vGn_n            => mpar_data%var(iLookPARAM%vGn_n)                ! van Genutchen "n" parameter (-)
 theta_sat        => mpar_data%var(iLookPARAM%theta_sat)            ! soil porosity (-)
 theta_res        => mpar_data%var(iLookPARAM%theta_res)            ! soil residual volumetric water content (-)
 vGn_m            => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)    ! van Genutchen "m" parameter (-)
 kappa            => mvar_data%var(iLookMVAR%scalarKappa)%dat(1)    ! constant in the freezing curve function (m K-1)

 ! assign pointers to index variables
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat  ! volumetric fraction of ice at the start of the iteration (-)
 mLayerTcrit      => mvar_data%var(iLookMVAR%mLayerTcrit)%dat       ! critical soil temperature above which all water is unfrozen (K) 
 layerType        => indx_data%var(iLookINDEX%layerType)%dat        ! layer type (ix_soil or ix_snow)

 ! identify the number of snow layers
 nSnow = count(layerType==ix_snow)

 ! update volumetric liquid and ice content (-)
 do iLayer=1,size(layerType)  ! (process snow and soil separately)
  ! compute liquid water equivalent of total water (liquid plus ice)
  theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
  !if(iLayer==nSnow+1) print*, 'in phseChange', mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer), theta, theta_sat
  select case(layerType(iLayer))

   ! ** snow
   case(ix_snow)
    ! compute the volumetric fraction of liquid water and ice (-)
    fLiq = fracliquid(mLayerTempNew(iLayer),snowfrz_scale)
    mLayerVolFracLiqNew(iLayer) = fLiq*theta
    mLayerVolFracIceNew(iLayer) = (theta - mLayerVolFracLiqNew(iLayer))*(iden_water/iden_ice)
    !write(*,'(a,1x,i4,1x,4(f20.10,1x))') 'in phase change: iLayer, fLiq, theta, mLayerVolFracLiqNew(iLayer), mLayerVolFracIceNew(iLayer) = ', &
    !                                                       iLayer, fLiq, theta, mLayerVolFracLiqNew(iLayer), mLayerVolFracIceNew(iLayer)

    ! avoid excessive change in a single iteration
    dLiq = mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer)
    if(abs(dLiq) > 0.05_dp)then
     !print*, 'modifying dLiq; dLiq = ', dLiq
     mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer) + dLiq*0.5_dp
     mLayerVolFracIceNew(iLayer) = (theta - mLayerVolFracLiqNew(iLayer))*(iden_water/iden_ice)
    endif

   ! ** soil
   case(ix_soil)
    ! constrain theta
    theta = min(theta,theta_sat)
    !print*, 'iLayer, mLayerVolFracIce(iLayer), mLayerTempNew(iLayer), mLayerTcrit(iLayer-nSnow) = ',&
    !         iLayer, mLayerVolFracIce(iLayer), mLayerTempNew(iLayer), mLayerTcrit(iLayer-nSnow)
    ! check that total volumetric water (liquid + ice) does not exceed soil porosity
    !if(theta > theta_sat)then; err=20; message=trim(message)//'volumetric (liquid + ice) content exceeds soil porosity'; return; endif 
    ! compute the matric head (m) volumetric fraction of liquid water and ice (-)
    if(mLayerTempNew(iLayer)<mLayerTcrit(iLayer-nSnow))then
     mLayerMatricHeadNew(iLayer-nSnow) = kappa*(mLayerTempNew(iLayer) - Tfreeze)
     mLayerVolFracLiqNew(iLayer)       = volFracLiq(mLayerMatricHeadNew(iLayer-nSnow),&
                                                    vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     mLayerVolFracIceNew(iLayer)       = (theta - mLayerVolFracLiqNew(iLayer))*(iden_water/iden_ice)
    else
     ! update matric head when all water is **unfrozen** -- if matric head > 0 at iter=m then no change in matric head
     !if(mLayerMatricHeadIter(iLayer-nSnow) > 0._dp)then ! saturated at the start of the iteration
     if(mLayerVolFracIce(iLayer) < tiny(theta))then ! no ice at the start of the iteration
      mLayerMatricHeadNew(iLayer-nSnow) = mLayerMatricHeadIter(iLayer-nSnow)
     else
      ! some water is frozen at the start of the iteration
      if(theta < epsilon(theta))then; err=20; message=trim(message)//'zero volumetric water content'; return; endif      
      mLayerMatricHeadNew(iLayer-nSnow) = matricHead(min(theta,theta_sat),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     endif
     ! update liquid water and ice content 
     mLayerVolFracLiqNew(iLayer)       = theta
     mLayerVolFracIceNew(iLayer)       = 0._dp
    endif

   ! ** check errors
   case default; err=10; message=trim(message)//'unknown case for model layer'; return

  endselect

  ! print results
  !if(iLayer > nSnow .and. iLayer < nSnow+3) &
  ! write(*,'(a,i4,1x,10(f20.10,1x))') 'in phase change: temp, liquid (iter, new), ice (iter, new, diff)', &
  !  iLayer, mLayerTempNew(iLayer), mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracIceIter(iLayer), mLayerVolFracIceNew(iLayer), &
  !  mLayerVolFracIceNew(iLayer) - mLayerVolFracIceIter(iLayer)


  ! sanity check
  if(mLayerVolFracIceNew(iLayer) < 0._dp)then
   write(message,'(a,i0,a,e20.10,a)')trim(message)//"volumetric ice content < 0 [iLayer=",iLayer,&
                                     &"; mLayerVolFracIceNew(iLayer)=",mLayerVolFracIceNew(iLayer),"]"
   err=10; return
  endif

 end do ! (looping through layers)
 endsubroutine phseChange

end module phseChange_module
