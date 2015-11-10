program integrator
  use polytrope_eos
  use cvode_indices
  use cvodehse
  
  implicit none

  ! Physics Parameters
  double precision, parameter :: poly_K = 1.982d-6
  double precision, parameter :: poly_gamma = 2.75d0

  ! Integration Data
  double precision, parameter :: R0 = 0.001d0
  double precision, parameter :: rho0 = 5d14
  double precision, dimension(2) :: Y0
  double precision, parameter :: RFIN = 30d5 ! Try to integrate to 30km
  double precision               :: R ! Holds the final R integrated to
  double precision, dimension(2) :: Y ! Holds the final values Y integrated to

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CVODE
  integer, parameter :: KEY = 1 ! Use CVODE
  ! NOTE: Below, integer*8 necessary to match the C int type expected
  integer*8, parameter :: NEQ = 2 ! Size of ODE system 
  integer ::            IER     ! Error flag

  integer, parameter ::   METH = 1 ! Use non-stiff Adams integration
  integer, parameter :: ITMETH = 2 ! Use Newton iteration
  integer, parameter ::  IATOL = 1 ! Use scalar absolute tolerance
  double precision, parameter :: RTOL = 1.0d-8 ! Relative tolerance
  double precision, parameter :: ATOL = 1.0d-8 ! Absolute tolerance

  integer, parameter :: MAX_NSTEPS = 10000 ! Maximum number of steps to solution
  
  integer, parameter :: SET_JAC_FLAG = 1 ! Use my supplied Jacobian
  integer, parameter :: ITASK = 1 ! Normal mode (overshoot and interpolate)

  integer, dimension(21) ::         IOUT ! Optional integer outputs
  double precision, dimension(6) :: ROUT ! Optional real outputs

  ! User integer parameters
  integer, dimension(1), parameter :: IPAR = (/ 0 /)
  ! User real parameters
  double precision, dimension(1), parameter :: RPAR = (/ 0.0d0 /) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Initialize the EOS
  call eos_init(poly_K, poly_gamma)

  ! Print Physical Parameters
  write(*,*) 'Integrating starting with:'
  write(*,*) 'R0 = ', R0
  write(*,*) 'Central Density, rho0 = ', rho0
  
  ! Print EOS parameters
  write(*,*) 'Initialized EOS with K = ', poly_K, ' and gamma = ', poly_gamma

  ! Initialize the Integration
  ! Mass (r=0) = 0.0
  Y0(iMASS) = 0.0d0
  ! Pressure (r=0) = EOS(rho0)
  call eos_rho(rho0) ! Compute EOS quantities at density rho0
  Y0(iPRES) = eos_data%p ! EOS Pressure for density rho0

  ! Print Solution IC's
  write(*,*) 'Initial Conditions:'
  write(*,*) 'Central Mass: ', Y0(iMASS)
  write(*,*) 'Central Pressure: ', Y0(iPRES)
  
  ! Setup CVODE
  call FNVINITS(KEY, NEQ, IER)
  call FCVMALLOC(R0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, &
       IOUT, ROUT, IPAR, RPAR, IER)
  call FCVSETIIN('MAX_NSTEPS', MAX_NSTEPS, IER)
  call FCVLAPACKDENSE(NEQ, IER)
  call FCVLAPACKDENSESETJAC(SET_JAC_FLAG, IER)
  
  ! Do the integration
  call FCVODE(RFIN, R, Y, ITASK, IER)

  ! Deallocate CVODE Memory
  call FCVFREE

  ! Print output
  write(*,*) 'Integrated to radius: ', R
  write(*,*) 'Final Mass: ', Y(iMASS)
  write(*,*) 'Final Pressure: ', Y(iPRES)
  
end program integrator
