module parameters
  implicit none

  character(len=*), parameter :: parameter_file_name = 'hse.par'
  
  type :: run_pars
     ! Input Physics Parameters
     double precision :: poly_K
     double precision :: poly_gamma
     double precision :: rho0
     double precision :: DR_SAVE ! Save solution each DR_SAVE (cm)
     integer :: NDR_SAVE         ! (Such that NDR_SAVE*DR_SAVE=R_TARGET=30d5)
  end type run_pars

  type(run_pars), save :: runtime_pars
  
contains
  subroutine get_run_parameters()
    character(len=50) :: label
    integer :: ios
    
    open(unit = 1, file=parameter_file_name)
    read(1, *, iostat=ios) label, runtime_pars%poly_K
    call ios_check_input(ios)
    read(1, *, iostat=ios) label, runtime_pars%poly_gamma
    call ios_check_input(ios)
    read(1, *, iostat=ios) label, runtime_pars%rho0
    call ios_check_input(ios)
    read(1, *, iostat=ios) label, runtime_pars%DR_SAVE
    call ios_check_input(ios)
    read(1, *, iostat=ios) label, runtime_pars%NDR_SAVE
    call ios_check_input(ios)
    close(unit=1)
    call status_run_parameters()
  end subroutine get_run_parameters

  subroutine ios_check_input(ios)
    integer, intent(in) :: ios
    if (ios /= 0) then
       write(*,*) 'IOS Error reading parameter file hse.par!'
       close(unit=1)
       stop
    end if
  end subroutine ios_check_input

  subroutine status_run_parameters()
    write(*,*) 'Input parameters are - '
    write(*,'(A,ES25.14)') 'poly_K: ', runtime_pars%poly_K
    write(*,'(A,ES25.14)') 'poly_gamma: ', runtime_pars%poly_gamma
    write(*,'(A,ES25.14)') 'rho0: ', runtime_pars%rho0
    write(*,'(A,ES25.14)') 'DR_SAVE: ', runtime_pars%DR_SAVE
    write(*,'(A,I25)') 'NDR_SAVE: ', runtime_pars%NDR_SAVE
  end subroutine status_run_parameters
end module parameters
