
module model_calibration

  USE nr_type, only: rkind
  USE parameter_search, only: parameter_spec,parameter_search_info
  USE parameter_fields, only: parameter_field

  implicit none
  private

  ! **************************************************************************************************
  ! Model-independent context for constructing calibrated parameter fields.
  !
  ! spec describes the complete parameter registry and dependencies.
  ! search defines the scalar calibration variables and their search bounds.
  ! fields contain the initialized and effective physical parameter values.
  ! values contains the untransformed scalar variables for the current trial.
  !
  ! The model-specific adapter extracts initialized parameter fields from
  ! the model and writes the effective fields back to the model.
  ! **************************************************************************************************
  type, public :: model_calibration_context
    type(parameter_spec)                  :: spec
    type(parameter_search_info)           :: search
    type(parameter_field), allocatable    :: fields(:)
    real(rkind), allocatable              :: values(:)
  end type model_calibration_context

end module model_calibration

