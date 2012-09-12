module build_dXdt_module
USE nrtype
implicit none
private
public::build_dXdt
contains

 ! **********************************************************************************************************
 ! new subroutine: build the state equations
 ! **********************************************************************************************************
 subroutine build_dXdt(dX_dt,err,message)
 USE data_struc,only:indx_data,ix_snow         ! model indices
 USE data_struc,only:mvar_data,stateLookup     ! model state variables
 USE var_lookup,only:iLookMVAR,iLookINDEX      ! named variables to define elements in the structure
 implicit none
 ! define dummy variables
 real(dp),intent(out)         :: dX_dt(:)      ! temporal derivative of model state vector
 integer(i4b),intent(out)     :: err
 character(*),intent(out)     :: message
 ! local pointers to data in structures
 real(dp),pointer             :: mLayerRadCondFlux(:)     ! total radiation and heat conduction flux (J m-2 s-1)
 real(dp),pointer             :: mLayerEnthalpyLiqFlux(:) ! energy flux associated with vertical flux of liquid water (J m-2 s-1)
 real(dp),pointer             :: mLayerMassLiqFlux(:)     ! total mass flux (kg m-2 s-1)
 ! local pointers to model indices
 integer(i4b),pointer         :: layerType(:)  ! type of the layer (ix_soil or ix_snow)
 integer(i4b),pointer         :: nLayTotl      ! total number of layers
 ! local variables
 integer(i4b)                 :: nLaySnow      ! number of snow layers
 integer(i4b)                 :: nState        ! number of state variables in the x-vector
 logical(lgt)                 :: maskEnthalpy  (size(dX_dt))
 logical(lgt)                 :: maskBulkDenWat(size(dX_dt))
 logical(lgt)                 :: maskDepth     (size(dX_dt))
 logical(lgt)                 :: maskAlbedo    (size(dX_dt))
 ! initialize error control
 err=0; message="build_dXdt/"
 ! initialize dX_dt to stupidly huge values (to force crash if elements not populated)
 dX_dt(:) = huge(dp)
 ! assign pointers to the model variables
 mLayerRadCondFlux     => mvar_data%var(iLookMVAR%mLayerRadCondFlux)%dat
 mLayerEnthalpyLiqFlux => mvar_data%var(iLookMVAR%mLayerEnthalpyLiqFlux)%dat
 mLayerMassLiqFlux     => mvar_data%var(iLookMVAR%mLayerMassLiqFlux)%dat
 ! assign pointers to total number of layers and the layer type
 nLayTotl  => indx_data%var(iLookINDEX%nLayers)%dat(1)
 layerType => indx_data%var(iLookINDEX%layerType)%dat
 ! get the number of snow layers
 nLaySnow = count(layerType==ix_snow)
 ! define the size of the x-vector
 nState   = nLayTotl*2 + nLaySnow + min(1,nLaySnow)*1
 if(size(stateLookup)/=nState)then; err=20; message=trim(message)//"unexpectedNumberVectorElements/"; return; endif
 ! get the masks
 maskEnthalpy   = (stateLookup==iLookMVAR%mLayerEnthalpy)
 maskBulkDenWat = (stateLookup==iLookMVAR%mLayerBulkDenWat)
 maskDepth      = (stateLookup==iLookMVAR%mLayerDepth)
 maskAlbedo     = (stateLookup==iLookMVAR%scalarAlbedo)
 ! check the masks
 if(count(maskEnthalpy  )/=nLayTotl)         then; err=30; message=trim(message)//"unexpectedSize[enthalpy]/"; return; endif
 if(count(maskBulkDenWat)/=nLayTotl)         then; err=30; message=trim(message)//"unexpectedSize[bulkDenWat]/"; return; endif
 if(count(maskDepth     )/=nLaySnow)         then; err=30; message=trim(message)//"unexpectedSize[depth]/"; return; endif
 if(count(maskAlbedo    )/=min(1,nLaySnow)*1)then; err=30; message=trim(message)//"unexpectedSize[albedo]/"; return; endif
 ! build the state equations
 if(count(maskEnthalpy  )>0) dX_dt = unpack(mLayerRadCondFlux + mLayerEnthalpyLiqFlux, maskEnthalpy, dX_dt)
 if(count(maskBulkDenWat)>0) dX_dt = unpack(mLayerMassLiqFlux, maskBulkDenWat, dX_dt)
 if(count(maskDepth     )>0) stop 'FORTRAN STOP: snow not implemented yet'
 if(count(maskAlbedo    )>0) stop 'FORTRAN STOP: snow not implemented yet'
 !print*, '1 ', mLayerRadCondFlux + mLayerEnthalpyLiqFlux
 !print*, '2 ', mLayerMassLiqFlux
 !print*, '3 ', dX_dt
 end subroutine build_dXdt

end module build_dXdt_module
