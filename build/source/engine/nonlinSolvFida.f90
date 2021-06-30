module nonlinSolvFida_module


  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use nrtype



    type nonlinSolvFn
        procedure(sysFnFida), pointer, nopass :: fcn
    end type nonlinSolvFn

    abstract interface

    integer(c_int) function sysFnFida(ycor, F, mem) result(ierr) bind(C,name='sysFnFida')
         use, intrinsic :: iso_c_binding
         use fida_mod                      ! Fortran interface to IDA
         use fnvector_serial_mod           ! Fortran interface to serial N_Vector
         use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
         use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
         use fsunmatrix_band_mod           ! Fortran interface to banded SUNMatrix
         use fsunlinsol_band_mod           ! Fortran interface to banded SUNLinearSolver
         use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver
         use fsundials_matrix_mod          ! Fortran interface to generic SUNMatrix
         use fsundials_nvector_mod         ! Fortran interface to generic N_Vector
         use fsundials_linearsolver_mod    ! Fortran interface to generic SUNLinearSolver
         use fsundials_nonlinearsolver_mod ! Fortran interface to generic SUNNonlinearSolver
         use nrtype
         use type4IDA
            implicit none
            type(N_Vector) :: ycor
            type(N_Vector) :: F
            type(c_ptr)    :: mem        
         end function

    end interface
  
end module nonlinSolvFida_module
