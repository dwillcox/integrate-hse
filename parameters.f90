module parameters
  implicit none

  character(len=*), parameter :: parameter_file_name = 'hse.par'
  integer, parameter :: parameter_file_unit = 1
  
contains
  subroutine open_parameters()
    open(unit=parameter_file_unit, recl=1024, delim='APOSTROPHE')
  end subroutine open_parameters

  subroutine close_parameters()
    close(unit=parameter_file_unit)
  end subroutine close_parameters
end module parameters

