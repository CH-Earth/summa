# Developer Guidelines for Contributing New Modular Components

New modular components may be added by using similar existing modular components as a template. The following steps (if applicable) may be used as a guideline. This process is illustrated using the addition of a new surface hydrology flux parameterization as an example.

## Identify a similar model component

* Identify the appropriate subdirectory within SUMMA's `source` directory such as:
    * `driver`: high-level program operations including the main driver
    * `dshare`: modules related to data storage and access
    * `engine`: low-level operations for physical and numerical processes
        * e.g., applies to flux calculations
* Identify the appropriate source file and module
    * source files have self-explanatory names
        * e.g., `soilLiqFlux.f90` corresponds to operations for liquid water fluxes in soil
    * each source file generally contains one module
        * e.g., `soilLiqFlux.f90` contains `soilLiqFlux_module`
* Identify the appropriate procedure
    * isolate the module procedure
        * e.g., within `soilLiqFlux_module`, the `surfaceFlux` module subroutine handles operations for surface hydrology fluxes
    * isolate the internal procedure
        * e.g., within the `contains` block of `surfaceFlux`, we have `update_surfaceFlux_prescribedHead` containing operations for specifying a prescribed pressure head surface boundary condition
        * `update_surfaceFlux_prescribedHead` may be used as a template for our example contribution
    * note that procedure names in SUMMA are organized using the terms *initialize*, *update*, and *finalize* to categorize operation types
        * *initialize* procedures are used for initial setup steps (initialization of variables, memory allocation, etc.)
        * *update* procedures are used for major computational operations (e.g., flux calculations)
        * *finalize* procedures are for post-processing and final error control checks

## Determine input and output variables
* Found by examining dummy variables in argument lists
    * Note that internal procedures inherit the dummy variables from the applicable module procedure by default
    * e.g., for the `update_surfaceFlux_prescribedHead` internal subroutine, the argument list of the `surfaceFlux` module subroutine applies: `subroutine surfaceFlux(io_soilLiqFlux,in_surfaceFlux,io_surfaceFlux,out_surfaceFlux)`
* Dummy variables may be objects with multiple data and procedure components
    * Such objects are declared using derived types (most commonly defined in `data_types.f90`)
    * Objects may be used to concisely interface data between the procedure and the caller
    * for SUMMA objects, the nomenclature `in_foobar`, `io_foobar`, and `out_foobar` is used for objects that interface input, input-output, and output data between the `foobar` procedure and its caller, respectively
* The `intent` attribute within dummy variable declarations indicates usage for input, input-output, or output
    * e.g., within `surfaceFlux` we have `type(in_type_surfaceFlux) ,intent(in)    :: in_surfaceFlux`, indicating the `in_surfaceFlux` object is for input data only
    * as noted above, the `in_surfaceFlux` object interfaces input data between the `surfaceFlux` module subroutine and its caller (the `soilLiqFlux` module subroutine)

## Create a skeleton of the new procedure
* Choose a self explanatory name for the new procedure
    * e.g. `update_surfaceFlux_example_flux`
* e.g., at the conclusion of this step, we would have a skeleton within the `contains` block of `surfaceFlux` similar to the following:

```fortran
subroutine update_surfaceFlux_example_flux
! main computations for the calculation of an example flux

end subroutine update_surfaceFlux_example_flux
```

* For new module procedures, the argument list from the template routine should be adjusted to match the needs of the new procedure
* For internal procedures (such as the example above), adding new data to interface objects from the corresponding module procedure may be required (see next step)

## Update derived type definitions for interface objects
* It may be desirable to add data components to existing objects related to the template procedure
    * e.g., adding a new numerical constant to be used in calculating a surface hydrology flux would require interfacing that data to the `update_surfaceFlux_example_flux` subroutine, which can be done using the `in_surfaceFlux` object
* Derived type definitions for interface objects are found in `source/dshare/data_types.f90`
    * e.g., for the `in_surfaceFlux` object we have the `in_type_surfaceFlux` derived type in the `data_types` module:
         
        ```fortran
        type, public :: in_type_surfaceFlux ! intent(in) data
          ! input: model control
          logical(lgt) :: firstSplitOper   ! flag indicating if desire to compute infiltration
          logical(lgt) :: deriv_desired    ! flag to indicate if derivatives are desired
          integer(i4b) :: bc_upper         ! index defining the type of boundary conditions
          integer(i4b) :: nRoots           ! number of layers that contain roots
          integer(i4b) :: ixIce            ! index of lowest ice layer
          integer(i4b) :: nSoil            ! number of soil layers
          ! [...] ! additional data components here
         contains
          procedure :: initialize => initialize_in_surfaceFlux
        end type in_type_surfaceFlux
        ```


* Adding a new numerical constant (say `example_flux_constant`) may be done as follows:

    ```fortran
    type, public :: in_type_surfaceFlux ! intent(in) data
      ! input: model control
      logical(lgt) :: firstSplitOper   ! flag indicating if desire to compute infiltration
      logical(lgt) :: deriv_desired    ! flag to indicate if derivatives are desired
      integer(i4b) :: bc_upper         ! index defining the type of boundary conditions
      integer(i4b) :: nRoots           ! number of layers that contain roots
      integer(i4b) :: ixIce            ! index of lowest ice layer
      integer(i4b) :: nSoil            ! number of soil layers
      ! input: values for the example flux
      real(rkind)  :: example_flux_constant ! numerical constant for example flux
      ! [...] ! additional data components here
     contains
      procedure :: initialize => initialize_in_surfaceFlux
    end type in_type_surfaceFlux
    ```

    * note that SUMMA uses the following `kind` parameters: `lgt` for logical variables, `i4b` for integer variables, and `rkind` for real variables

* Additionally, we have procedure components for *initialize* and *finalize* operations for data interfacing
    * e.g., `call in_surfaceFlux % initialize` points to the `initialize_in_surfaceFlux` class procedure (in the `contains` block of the `data_types` module) for initializing data components:

        ```fortran
        soutine initialize_in_surfaceFlux(in_surfaceFlux,nRoots,ixIce,nSoil,nGlce,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,&
                                           &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                           &iLayerHeight,dHydCond_dTemp,iceImpedeFac)
         class(in_type_surfaceFlux),intent(out) :: in_surfaceFlux ! input object for surfaceFlux
         ! [...] ! additional variable declarations here

         associate(&
          ! model control
          firstSplitOper         => in_soilLiqFlux % firstSplitOper,                      & ! flag to compute infiltration
          deriv_desired          => in_soilLiqFlux % deriv_desired,                       & ! flag indicating if derivatives are desired
          ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision & ! index defining the type of boundary conditions
         &)
          ! intent(in): model control
          in_surfaceFlux % firstSplitOper = firstSplitOper          ! flag indicating if desire to compute infiltration
          in_surfaceFlux % deriv_desired  = deriv_desired           ! flag indicating if derivatives are desired
          in_surfaceFlux % bc_upper       = ixBcUpperSoilHydrology  ! index defining the type of boundary conditions (Neumann or Dirichlet)
          in_surfaceFlux % nRoots         = nRoots                  ! number of layers that contain roots
          in_surfaceFlux % ixIce          = ixIce                   ! index of lowest ice layer
          in_surfaceFlux % nSoil          = nSoil                   ! number of soil layers
          in_surfaceFlux % nGlce          = nGlce                   ! number of glacier ice layers
         end associate

         ! [...] ! additional associate blocks here
        end subroutine initialize_in_surfaceFlux
        ```
        
    * new data components, such as `example_flux_constant`, must be applied within the procedure components:

        ```fortran
        subroutine initialize_in_surfaceFlux(in_surfaceFlux,nRoots,ixIce,nSoil,nGlce,ibeg,iend,in_soilLiqFlux,io_soilLiqFlux,&
                                           &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                           &iLayerHeight,dHydCond_dTemp,iceImpedeFac,example_flux_constant)
         class(in_type_surfaceFlux),intent(out) :: in_surfaceFlux ! input object for surfaceFlux
         ! [...] ! additional variable declarations here
         real(rkind),intent(in) :: example_flux_constant ! declaration for new constant

         associate(&
          ! model control
          firstSplitOper         => in_soilLiqFlux % firstSplitOper,                      & ! flag to compute infiltration
          deriv_desired          => in_soilLiqFlux % deriv_desired,                       & ! flag indicating if derivatives are desired
          ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision & ! index defining the type of boundary conditions
         &)
          ! intent(in): model control
          in_surfaceFlux % firstSplitOper = firstSplitOper          ! flag indicating if desire to compute infiltration
          in_surfaceFlux % deriv_desired  = deriv_desired           ! flag indicating if derivatives are desired
          in_surfaceFlux % bc_upper       = ixBcUpperSoilHydrology  ! index defining the type of boundary conditions (Neumann or Dirichlet)
          in_surfaceFlux % nRoots         = nRoots                  ! number of layers that contain roots
          in_surfaceFlux % ixIce          = ixIce                   ! index of lowest ice layer
          in_surfaceFlux % nSoil          = nSoil                   ! number of soil layers
          in_surfaceFlux % nGlce          = nGlce                   ! number of glacier ice layers
         end associate
         
         ! [...] ! additional associate blocks here
        
         ! assignment statements for the example flux
         in_surfaceFlux % example_flux_constant = example_flux_constant ! numerical constant for example flux
        
        end subroutine initialize_in_surfaceFlux
        ```

    * for the above example, we have added a dummy variable for the new example flux constant and an assignment statement to initialize the new data component `in_surfaceFlux % example_flux_constant`
    * note that the corresponding call to `in_surfaceFlux` within the `soilLiqFlux` subroutine would need to be updated to include the additional argument `example_flux_constant`

## Add operations to the skeleton procedure
* add main operations within the skeleton procedure created in the above steps to complete the new modular component
    * for the example flux parameterization (using a toy model of constant infiltration), we have:
     
        ````fortran
        subroutine update_surfaceFlux_example_flux
         ! main computations for the calculation of an example flux
         associate(&
          ! input: flux at the upper boundary
          scalarRainPlusMelt => in_surfaceFlux % scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
          ! input: numerical constants
          example_flux_constant => in_surfaceFlux % example_flux_constant
          ! output: runoff and infiltration
          scalarSurfaceRunoff       => out_surfaceFlux % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
          scalarSurfaceInfiltration => io_surfaceFlux % scalarSurfaceInfiltration   & ! surface infiltration (m s-1)
         &)
          scalarSurfaceInfiltration = example_flux_constant                          ! toy model of constant infiltration
          scalarSurfaceRunoff       = scalarRainPlusMelt - scalarSurfaceInfiltration ! compute surface runoff
         end associate
        end subroutine update_surfaceFlux_example_flux
        ````