module MapStr2vec_module
USE nrtype
implicit none
private
public::str2vector
public::vector2str
contains

 ! **********************************************************************************************************
 ! new subroutine: extract data from structures and put it in a vector
 ! **********************************************************************************************************
 subroutine str2vector(stateVec,err,message)
 USE data_struc,only:indx_data,ix_soil,ix_snow     ! model indices
 USE data_struc,only:mvar_data,stateLookup         ! model state variables
 USE var_lookup,only:iLookMVAR,iLookINDEX          ! named variables to define elements in the structure
 implicit none
 ! define dummy variables
 real(dp),pointer,intent(out) :: stateVec(:)
 integer(i4b),intent(out)     :: err
 character(*),intent(out)     :: message
 ! local variables
 integer(i4b),pointer         :: layerType(:)      ! type of the layer (ix_soil or ix_snow)
 integer(i4b),pointer         :: nLayTotl          ! total number of layers
 integer(i4b)                 :: nLaySoil          ! number of soil layers
 integer(i4b)                 :: nLaySnow          ! number of snow layers
 integer(i4b)                 :: nState            ! number of state variables in the x-vector

 ! initialize error control
 err=0; message="f-fuse-str2vector/"
 ! get the total number of layers and the layer type
 nLayTotl  => indx_data%var(iLookINDEX%nLayers)%dat(1)
 layerType => indx_data%var(iLookINDEX%layerType)%dat
 ! check the size of data structures match the number of layers
 if(size(mvar_data%var(iLookMVAR%mLayerTemp      )%dat)/=nLayTotl .or. &
    size(mvar_data%var(iLookMVAR%mLayerBulkDenIce)%dat)/=nLayTotl .or. &
    size(mvar_data%var(iLookMVAR%mLayerBulkDenLiq)%dat)/=nLayTotl .or. &
    size(mvar_data%var(iLookMVAR%mLayerMatricHead)%dat)/=nLayTotl .or. &
    size(mvar_data%var(iLookMVAR%mLayerDepth     )%dat)/=nLayTotl)then
  err=20; message=trim(message)//"unexpectedNumberVectorElements/"; return
 endif
 ! get the number of snow layers
 nLaySoil = count(layerType==ix_soil)
 nLaySnow = count(layerType==ix_snow)
 ! define the size of the x-vector
 nState   = nLaySoil*4 + &    ! temperature, bulk density of ice, bulk density of liquid water, and matric head
            nLaySnow*4 + &    ! temperature, bulk density of ice, bulk density of liquid water, and depth
            min(1,nLaySnow)*1 ! snow albedo (if snow exists)
 ! allocate space for the x-vector
 if (associated(stateVec))    deallocate(stateVec)
 if (associated(stateLookup)) deallocate(stateLookup)
 allocate(stateVec(nState),stateLookup(nState),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocate/stateVec"; return; endif
 ! populate state vectors and look-up tables
 if(nLaySnow==0)then ! ***** no snow layers *****
                     ! --------------------------
  ! populate the state vector
  stateVec    = (/mvar_data%var(iLookMVAR%mLayerTemp      )%dat(1:nLaySoil), &  ! temperature
                  mvar_data%var(iLookMVAR%mLayerBulkDenIce)%dat(1:nLaySoil), &  ! bulk density of ice
                  mvar_data%var(iLookMVAR%mLayerBulkDenLiq)%dat(1:nLaySoil), &  ! bulk density of liquid water
                  mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(1:nLaySoil)/)   ! matric head
  ! populate a look-up table with named variables describing the data in the state vector
  stateLookup = (/spread(iLookMVAR%mLayerTemp,      1,nLaySoil),  &  ! temperature
                  spread(iLookMVAR%mLayerBulkDenIce,1,nLaySoil),  &  ! bulk density of ice
                  spread(iLookMVAR%mLayerBulkDenLiq,1,nLaySoil),  &  ! bulk density of liquid water
                  spread(iLookMVAR%mLayerMatricHead,1,nLaySoil)/)    ! matric head
 else                ! ***** snow layers exist *****
                     ! -----------------------------
  ! populate the state vector
  stateVec    = (/mvar_data%var(iLookMVAR%mLayerTemp      )%dat(1:nLayTotl), &                    ! temperature (snow and soil)
                  mvar_data%var(iLookMVAR%mLayerBulkDenIce)%dat(1:nLayTotl), &                    ! bulk density of ice (snow and soil)
                  mvar_data%var(iLookMVAR%mLayerBulkDenLiq)%dat(1:nLayTotl), &                    ! bulk density of liquid water (snow and soil)
                  mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(1:nLaySoil), &                    ! matric head (soil only)
                  mvar_data%var(iLookMVAR%mLayerDepth     )%dat(nLaySoil+1:nLaySoil+nLaySnow), &  ! layer depth (snow only)
                  mvar_data%var(iLookMVAR%scalarAlbedo    )%dat(1)/)                              ! snow albedo (scalar, only if snow is present)
  ! populate a look-up table with named variables describing the data in the state vector
  stateLookup = (/spread(iLookMVAR%mLayerTemp,      1,nLayTotl),  &  ! temperature
                  spread(iLookMVAR%mLayerBulkDenIce,1,nLaySnow),  &  ! bulk density of ice
                  spread(iLookMVAR%mLayerBulkDenLiq,1,nLaySnow),  &  ! bulk density of liquid water
                  spread(iLookMVAR%mLayerMatricHead,1,nLaySoil),  &  ! matric head
                  spread(iLookMVAR%mLayerDepth,     1,nLaySnow),  &  ! depth
                         iLookMVAR%scalarAlbedo /)                   ! albedo
 endif
 print*,'stateLookup = ',stateLookup
 print*,'stateVec    = ',stateVec
 stop
 end subroutine str2vector

 ! **********************************************************************************************************
 ! new subroutine: assign the data structures to the data in the vector
 ! **********************************************************************************************************
 subroutine vector2str(stateVec,err,message)
 USE data_struc,only:indx_data,ix_soil,ix_snow       ! model indices
 USE data_struc,only:mpar_data,mvar_data,stateLookup ! model parameters and state variables
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX ! named variables to define elements in the structure
 implicit none
 ! define dummy variables
 real(dp),intent(in)         :: stateVec(:)
 integer(i4b),intent(out)    :: err
 character(*),intent(out)    :: message
 ! local pointers to data in structures
 integer(i4b),pointer        :: layerType(:)         ! type of the layer (ix_soil or ix_snow)
 integer(i4b),pointer        :: nLayTotl             ! total number of layers
 integer(i4b)                :: nLaySoil             ! number of soil layers
 integer(i4b)                :: nLaySnow             ! number of snow layers
 integer(i4b)                :: nState               ! number of state variables in the x-vector
 logical(lgt)                :: maskTemp      (size(stateVec))
 logical(lgt)                :: maskBulkDenIce(size(stateVec))
 logical(lgt)                :: maskBulkDenLiq(size(stateVec))
 logical(lgt)                :: maskMatricHead(size(stateVec))
 logical(lgt)                :: maskDepth     (size(stateVec))
 logical(lgt)                :: maskAlbedo    (size(stateVec))
 ! initialize error control
 err=0; message="fuse-vector2str/"
 ! get the total number of layers and the layer type
 nLayTotl  => indx_data%var(iLookINDEX%nLayers)%dat(1)
 layerType => indx_data%var(iLookINDEX%layerType)%dat
 ! get the number of soil and snow layers
 nLaySoil = count(layerType==ix_soil)
 nLaySnow = count(layerType==ix_snow)
 ! define the size of the x-vector
 nState   = nLaySoil*4 + &    ! temperature, bulk density of ice, bulk density of liquid water, and matric head
            nLaySnow*4 + &    ! temperature, bulk density of ice, bulk density of liquid water, and depth
            min(1,nLaySnow)*1 ! snow albedo (if snow exists)
 ! check the x-vector is appropriately sized
 if(size(stateVec)/=nState)then; err=20; message=trim(message)//"unexpectedNumberVectorElements/"; return; endif
 ! get the masks
 maskTemp       = (stateLookup==iLookMVAR%mLayerTemp)
 maskBulkDenIce = (stateLookup==iLookMVAR%mLayerBulkDenIce)
 maskBulkDenLiq = (stateLookup==iLookMVAR%mLayerBulkDenLiq)
 maskMatricHead = (stateLookup==iLookMVAR%mLayerMatricHead)
 maskDepth      = (stateLookup==iLookMVAR%mLayerDepth)
 maskAlbedo     = (stateLookup==iLookMVAR%scalarAlbedo)
 ! check the masks
 if(count(maskTemp      )/=nLayTotl)         then; err=30; message=trim(message)//"unexpectedSize[Temp]/"; return; endif
 if(count(maskBulkDenIce)/=nLayTotl)         then; err=30; message=trim(message)//"unexpectedSize[BulkDenIce]/"; return; endif
 if(count(maskBulkDenLiq)/=nLayTotl)         then; err=30; message=trim(message)//"unexpectedSize[BulkDenLiq]/"; return; endif
 if(count(maskMatricHead)/=nLaySoil)         then; err=30; message=trim(message)//"unexpectedSize[matricHead]/"; return; endif
 if(count(maskDepth     )/=nLaySnow)         then; err=30; message=trim(message)//"unexpectedSize[Depth]/"; return; endif
 if(count(maskAlbedo    )/=min(1,nLaySnow)*1)then; err=30; message=trim(message)//"unexpectedSize[Albedo]/"; return; endif
 ! assign the data structures to the data in the vector
 if(count(maskTemp      )>0) mvar_data%var(iLookMVAR%mLayerTemp      )%dat(1:nLayTotl) = pack(stateVec,maskTemp)
 if(count(maskBulkDenIce)>0) mvar_data%var(iLookMVAR%mLayerBulkDenIce)%dat(1:nLayTotl) = pack(stateVec,maskBulkDenIce)
 if(count(maskBulkDenLiq)>0) mvar_data%var(iLookMVAR%mLayerBulkDenLiq)%dat(1:nLayTotl) = pack(stateVec,maskBulkDenLiq)
 if(count(maskMatricHead)>0) mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(1:nLaySoil) = pack(stateVec,maskMatricHead)
 if(count(maskDepth     )>0) mvar_data%var(iLookMVAR%mLayerDepth     )%dat(nLaySoil+1:nLaySoil+nLaySnow) = pack(stateVec,maskDepth)
 if(count(maskAlbedo    )>0) mvar_data%var(iLookMVAR%scalarAlbedo    )%dat(1:1)        = pack(stateVec,maskAlbedo)
 print*,maskDepth
 print*,mvar_data%var(iLookMVAR%mLayerDepth     )%dat
 end subroutine vector2str

end module MapStr2Vec_module
