# SUMMA Coding Conventions

## Variable names

 1. Use self-describing variable names, even if they have a large number of characters. For example, `canopyEvap` is preferable to `ce`

 1. Start with a lowercase and then for each word break use a capital letter. This is known in the trade as *camelBack*

 1. Define every single variable, including units, with each variable on a separate line. That is
    ```fortran
    real(dp),intent(out) :: canopyEvap ! canopy evaporation (kg m-2 s-1)
    ```
 1. For local instantiations of a global variable, specify the variable dimension as the leading part of the variable name. For example, `scalarAlbedo`, `mLayerDepth`, `iLayerHydCond`, where `scalar` defines a scalar variable, `mLayer` defines variables defined at the mid-point of vertical layers, and `iLayer` defines a variable defined at the interfaces of model layers.

 1. All variables need to have an explicit type declaration statement. This should be enforced through the use of

    ```fortran
    implicit none
    ```

## Hard-coded numbers, physical constants and parameters

 1. Do not use hard-coded numbers in equations to represent physical constants or model parameters.

 1. All physical constants are defined in `/build/source/dshare/multiconst.f90` and should *not* be redefined elsewhere in the code. If you add physical constants, they must be added to `/build/source/dshare/multiconst.f90`.

 1. All model parameters are ideally read in from model configuration files to allow model users to experiment with parameter values without requiring code changes. During model development, it it sometimes easier to skip this step. In that case, at the very least, model parameters must be identified at the top of a routine as part of the variable definition section, e.g.
    ```fortran
    real(dp),parameter              :: facTrustDec=0.25_dp          ! factor decrease in the trust region
    ```
 All the coding conventions for variable names apply in that case.

 When development is _code complete_, all model parameters must be converted to inputs for a pull request to be accepted.

## Commenting

 1. Include a comment block at the start of each subroutine or function, providing a brief description of what the subroutine does (including references). To indicate that this is the start of a subroutine or function, the comment should start with `! <scope> subroutine <subroutine name>: <description>` or `! <scope> function <function name>: <description>`, where `<scope>` is `private|public|internal`. The line above and after this one lines has the form `! *****`. For example:

    ```fortran
     ! ************************************************************************************************
     ! public subroutine init_metad: initialize metadata structures
     ! ************************************************************************************************
    ```

    If there is a more detailed description (which is encouraged), then continue with comments after the above header and end with another line of the form `! *****`, for example:

    ```fortran
     ! ************************************************************************************************
     ! public subroutine groundwatr: compute the groundwater sink term in Richards' equation
     ! ************************************************************************************************
     !
     ! Method
     ! ------
     !
     ! Here we assume that water avaialble for shallow groundwater flow includes is all water above
     ! "field capacity" below the depth zCrit, where zCrit is defined as the lowest point in the soil
     ! profile where the volumetric liquid water content is less than field capacity.
     !
     ! We further assume that transmssivity (m2 s-1) for each layer is defined asuming that the water
     ! available for saturated flow is located at the bottom of the soil profile. Specifically:
     !  trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
     !  trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
     ! where zActive(iLayer) is the effective water table thickness for all layers up to and including
     ! the current layer (working from the bottom to the top).
     !
     ! The outflow from each layer is then (m3 s-1)
     !  mLayerOutflow(iLayer) = trSoil(iLayer)*tan_slope*contourLength
     ! where contourLength is the width of a hillslope (m) parallel to a stream
     !
     ! ************************************************************************************************
    ```

 1. Do not commit code with large sections of code commented out. The version control system should be used for tracking different versions. There is no need to do it manually by commenting out code, which is much more likely to lead to confusion. The only exception is for some debug statements that may be uncommented and reused during debugging.

 1. Comment often. Strive for a comment on every line of code. It takes much less time to "*comment as you go*" than try and comment afterwards.

 3. When there is potential for confusion, define units for the LHS of the assignment

## Multiple statements on a single line

Do not use multiple statements on a single line (separated by a `;`) except in two cases:

 1. Error handling, for example

    ```fortran
    if(err>0)then; message=trim(message)//trim(cmessage); return; endif
    ```

 1. Case statements, for example

    ```fortran
    case('scalarBartDummy'                ); get_ixmvar = iLookMVAR%scalarBartDummy                  ! dummy variable for bart (-)
    case('scalarCosZenith'                ); get_ixmvar = iLookMVAR%scalarCosZenith                  ! cosine of the solar zenith angle (0-1)
    ```

## Indenting

Use one space to indent

 * subroutines within a module
 * if statements (if then else)
 * do constructs

## Subroutines

 1. Organize subroutines into modules (even if the module contains just one subroutine). This avoids the need to write an explicit interface for each subroutine.

    ```fortran
module vegLiqFlux_module
USE nrtype
implicit none
private
public::vegLiqFlux
contains


 ! ************************************************************************************************
 ! public subroutine vegLiqFlux: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine vegLiqFlux(&
                       ! input
                       computeVegFlux,               & ! intent(in): flag to denote if computing energy flux over vegetation
                       scalarCanopyLiqTrial,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       scalarRainfall,               & ! intent(in): rainfall rate (kg m-2 s-1)
                       ! output
                       scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                       err,message)                    ! intent(out): error control
    ```

 1. Make subroutines and functions as private as possible. Only make subroutines and functions public if absolutely necessary. Subroutines and functions in another module are accessed with a `USE` statement with the `only` attribute. For example

    ```fortran
 USE vegliqflux_module,only:vegliqflux                ! compute liquid water fluxes through vegetation
    ```

 1. The `intent` argument must be used when passing variables to subroutines.

## Whitespace

 1. Two whitelines before the comment block that starts a new subroutine or function. No whitelines between the comment block and the start of the subroutine or function.

 1. Two whitelines before the `end module` statement.

 1. In all other cases use a single whiteline to indicate sections of your code.

 1. Strip all trailing whitespace (spaces at the end of a line or on a whiteline). Most editors have a setting that will automatically do this when a file is saved. This will make for much cleaner commits.

## Licensing

 1. SUMMA is distributed under the [GPL-V3 license](http://www.gnu.org/licenses/gpl-3.0.html). This means that all code modifications must be shared with the community. All code files should include a header that contains the following (this should be copied verbatim other than the text in `{}`). A header template is provided in the top-level SUMMA directory as `header.license`.

    ```fortran
    ! SUMMA - Structure for Unifying Multiple Modeling Alternatives
    ! Copyright (C) 2014-{year}  NCAR/RAL
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
    ```
