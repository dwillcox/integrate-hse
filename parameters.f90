module parameters
  use eos
  use CVODE_PARAMETERS
  
  implicit none

  character(len=*), parameter :: parameter_file_name = 'hse.par'
  integer, parameter :: parameter_file_unit = 1
  
contains
  subroutine init_parameters()
    open(unit=parameter_file_unit, file=parameter_file_name, recl=1024, delim='APOSTROPHE')

    call cvode_init(parameter_file_unit)
    call eos_init(parameter_file_unit)
    
    close(unit=parameter_file_unit)
  end subroutine init_parameters

end module parameters

