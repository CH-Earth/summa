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

module summa_defineOutput                     ! used to define model output files

! named variables to define new output files
USE globalData, only: noNewFiles              ! no new output files
USE globalData, only: newFileEveryOct1        ! create a new file on Oct 1 every year (start of the USA water year)

! metadata structures
USE globalData,only:attr_meta                 ! attributes metadata structure
USE globalData,only:type_meta                 ! veg/soil type metadata structure
USE globalData,only:mpar_meta                 ! local parameter metadata structure
USE globalData,only:bpar_meta                 ! basin parameter metadata structure
USE globalData,only:grid_meta                 ! grid metadata structure

! named variables
USE var_lookup,only:iLookTIME                 ! named variables for time data structure
USE var_lookup,only:iLookFREQ                 ! named variables for the frequency structure

! safety: set private unless specified otherwise
implicit none
private
public::summa_defineOutputFiles
contains

 ! used to define model output files
 subroutine summa_defineOutputFiles(modelTimeStep, using_buffer, summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nr_type                                                 ! variable types, etc.
 USE summa_type, only:summa1_type_dec                        ! master summa data type
 ! functions and subroutines
 USE def_output_module,only:def_output                       ! module to define model output
 USE modelwrite_module,only:writeParam,writeGridParam        ! module to write model parameters
 ! global data structures
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:structInfo                              ! information on the data structures
 ! file information
 USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
 USE globalData,only:output_fileSuffix                       ! suffix for the output file
 USE globalData,only:newOutputFile                           ! define option for new output files
 USE globalData,only:fileout                                 ! name of the output file
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(in)               :: modelTimeStep        ! time step index
 logical(lgt),intent(in)               :: using_buffer         ! flag for will do buffered write
 type(summa1_type_dec),intent(inout)   :: summa1_struc         ! master summa data structure
 integer(i4b),intent(out)              :: err                  ! error code
 character(*),intent(out)              :: message              ! error message
 ! local variables
 character(LEN=256)                    :: cmessage             ! error message of downwind routine
 integer(i4b)                          :: iGRU,iHRU,iDOM,iGrid ! indices of GRUs, HRUs, domains, and grids
 integer(i4b)                          :: iStruct              ! index of model structure
 ! version information generated during compiling
 INCLUDE 'summaversion.inc'
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&
  timeStruct           => summa1_struc%timeStruct        , & ! x%var(:)               -- model time data
  attrStruct           => summa1_struc%attrStruct        , & ! x%gru(:)%hru(:)%var(:) -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct        , & ! x%gru(:)%hru(:)%var(:) -- local classification of soil veg etc. for each HRU
  idStruct             => summa1_struc%idStruct          , & ! x%gru(:)%hru(:)%var(:) -- local classification of soil veg etc. for each HRU
  mparStruct           => summa1_struc%mparStruct        , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model parameters
  bparStruct           => summa1_struc%bparStruct        , & ! x%gru(:)%var(:)                   -- basin-average parameters
  gridStruct           => summa1_struc%gridStruct        , & ! x%gru(:)%var(:)%var(:)%dat2(:,:)  -- basin grid parameters and variables
  nGRU                 => summa1_struc%nGRU              , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU                & ! number of global hydrologic response units
 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_defineOutputFiles/'

 ! *****************************************************************************
 ! *** define the name of the model output file
 ! *****************************************************************************

 ! define full name of output file
 if(modelTimeStep==1)then
  select case(newOutputFile)
   case(noNewFiles);          ! do nothing, just ensure validity of outputfile option
   case(newFileEveryOct1);
   case default; err=20; message=trim(message)//'unable to identify the option to define new output files'; return
  end select

  fileout = trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//trim(output_fileSuffix)

 ! model time step > 1: define name of output file : new simulations
 else
  write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX),&
        timeStruct%var(iLookTIME%iyyy),'-',timeStruct%var(iLookTIME%iyyy)+1,&
        trim(output_fileSuffix)
 endif

 ! *****************************************************************************
 ! *** define the model output file and write parameters
 ! *****************************************************************************

 ! define the file
 call def_output(using_buffer,summaVersion,buildTime,gitBranch,gitHash,nGRU,nHRU,fileout,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! write parameters with no time dimension
 do iGRU=1,nGRU

  ! write HRU parameters, all written to timestep frequency file
  do iHRU=1,gru_struc(iGRU)%hruCount
   do iStruct=1,size(structInfo)
    select case(trim(structInfo(iStruct)%structName))
     case('attr'); call writeParam(0,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,attrStruct%gru(iGRU)%hru(iHRU),attr_meta,err,cmessage)
     case('type'); call writeParam(0,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,typeStruct%gru(iGRU)%hru(iHRU),type_meta,err,cmessage)
     case('mpar')
      do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
       call writeParam(iDOM,gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,mparStruct%gru(iGRU)%hru(iHRU)%dom(iDOM),mpar_meta,err,cmessage)
      end do
    end select
    if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
   end do  ! (looping through structures)
  end do  ! (looping through HRUs)

  ! write GRU parameters, written to timestep (bpar) or annual (grid) frequency files
  do iStruct=1,size(structInfo)
   select case(trim(structInfo(iStruct)%structName))
    case('bpar'); call writeParam(0,iGRU,bparStruct%gru(iGRU),bpar_meta,err,cmessage) ! write to timestep frequency file
    case('grid')
      if(.not.using_buffer)then ! these params write to the annual frequency file, which is not made if using a buffered write
        do iGrid=1,gru_struc(iGRU)%nGrid ! will not write if there are no grids
          call writeGridParam(iGRU,iGrid,gridStruct%gru(iGRU)%grid(iGrid),grid_meta,err,cmessage)
        end do
      end if
   end select
   if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
  end do  ! (looping through structures)

 end do  ! (looping through GRUs)

 ! end associate statements
 end associate summaVars

 end subroutine summa_defineOutputFiles
end module summa_defineOutput
