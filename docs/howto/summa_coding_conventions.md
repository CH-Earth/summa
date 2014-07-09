*SUMMA Coding Conventions*

**Variable names**

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

**Commenting**

 1. Include a comment block at the start of each subroutine, providing a brief description of what the subroutine does (including references).

 1. Do not commit code with large sections of code commented out. The version control system should be used for tracking different versions. There is no need to do it manually by commenting out code, which is much more likely to lead to confusion.

 1. Comment often. Strive for a comment on every line of code. It takes much less time to "*comment as you go*" than try and comment afterwards.

 3. When there is potential for confusion, define units for the LHS of the assignment

**Multiple statements on a single line**

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

**Indenting**

Use one space to indent

 * subroutines within a module
 * if statements (if then else)
 * do constructs

**Subroutines**

 1. Organize subroutines into modules (even if the module contains just one subroutine). This avoids the need to write an explicit interface for each subroutine.

    ```fortran
module vegLiqFlux_module
USE nrtype
implicit none
private
public::vegLiqFlux
contains

 ! ************************************************************************************************
 ! new subroutine: compute water balance for the vegetation canopy
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

 1. Make subroutines as private as possible. Only make subroutines public if absolutely necessary. Subroutines in another module are accessed with a `USE` statement with the `only` attribute. For example

    ```fortran
 USE vegliqflux_module,only:vegliqflux                ! compute liquid water fluxes through vegetation
    ```

 1. The `intent` argument must be used when passing variables to subroutines.

 **Licensing**

 1. SUMMA is distributed under the GPL-V3 license. This means that all code modifications must be shared with the community. All code files should include a header that contains the following (this should be copied verbatim other than the text in `{}`).

    ```fortran
    ! SUMMA - Structure for Unifying Multiple Modeling Alternatives
    ! Copyright (C) {year}  {name of author}
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

