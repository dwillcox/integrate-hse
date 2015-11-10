program integrator
  use polytrope_eos
  use cvode_indices
  use cvodehse
  use data_wrangler
  use parameters
  
  implicit none

  ! Integration Data Structures and Parameters
  double precision, parameter :: R0 = 0.001d0
  double precision, dimension(2) :: Y0
  double precision :: R_TARGET ! Try to integrate to RFIN = NDRV*DRSV
  double precision :: RK ! Stores the radius at the current interval
  integer :: K ! Stores which interval we're in (count to NDR_SAVE)
  integer :: KFIN ! Stores how many intervals we went (K = 1 to KFIN)
  double precision :: DR_SAVE ! Save solution each DR_SAVE (cm)
  integer :: NDR_SAVE ! (Such that NDR_SAVE*DR_SAVE=R_TARGET=30d5)
  double precision               :: R ! Holds the R integrated to
  double precision, dimension(2) :: Y ! Holds the solution Y integrated to
  double precision, dimension(2) :: DKY ! Holds the p=0 root solution for Y
  double precision, dimension(:,:), allocatable :: hse_profile

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

  integer, parameter :: NRTFN = 1 ! Number of root functions to solve during CVODE
  
  ! User integer parameters
  integer, dimension(1), parameter :: IPAR = (/ 0 /)
  ! User real parameters
  double precision, dimension(1), parameter :: RPAR = (/ 0.0d0 /) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read runtime parameters
  call get_run_parameters()
  DR_SAVE = runtime_pars%DR_SAVE
  NDR_SAVE = runtime_pars%NDR_SAVE

  ! Allocate Profile array
  allocate( hse_profile(NPROFILE, NDR_SAVE) )
  
  ! Initialize the EOS
  call eos_init(runtime_pars%poly_K, runtime_pars%poly_gamma)

  ! Print Physical Parameters
  write(*,*) 'Integrating starting with:'
  write(*,'(A,ES25.14)') 'R0 = ', R0
  write(*,'(A,ES25.14)') 'Central Density, rho0 = ', runtime_pars%rho0
  
  ! Print EOS parameters
  write(*,'(A,ES25.14,A,ES25.14)') 'Initialized EOS with K = ', runtime_pars%poly_K, ' and gamma = ', runtime_pars%poly_gamma

  ! Initialize the Integration
  ! Mass (r=0) = 0.0
  Y0(jmass) = 0.0d0
  ! Pressure (r=0) = EOS(rho0)
  call eos_rho(runtime_pars%rho0) ! Compute EOS quantities at density rho0
  Y0(jpres) = eos_data%p ! EOS Pressure for density rho0
  ! Integration target R_TARGET = NDR_SAVE*DR_SAVE
  R_TARGET = DBLE(NDR_SAVE)*DR_SAVE

  ! Print Solution IC's
  write(*,*) 'Initial Conditions:'
  write(*,'(A,ES25.14)') 'Central Mass: ', Y0(jmass)
  write(*,'(A,ES25.14)') 'Central Pressure: ', Y0(jpres)
  
  ! Setup CVODE
  call FNVINITS(KEY, NEQ, IER)
  call FCVMALLOC(R0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, &
       IOUT, ROUT, IPAR, RPAR, IER)
  if (IER /= 0) then
     write(*,*) 'Error occurred in FCVMALLOC!'
     stop
  end if
  call FCVROOTINIT(NRTFN, IER)
  if (IER /= 0) then
     write(*,*) 'Error occurred in FCVROOTINIT!'
     stop
  end if
  call FCVSETIIN('MAX_NSTEPS', MAX_NSTEPS, IER)
  call FCVLAPACKDENSE(NEQ, IER)
  call FCVLAPACKDENSESETJAC(SET_JAC_FLAG, IER)
  
  ! Do the integration
  do K = 1, NDR_SAVE
     RK = min(R0 + DBLE(K)*DR_SAVE, R_TARGET)
     call FCVODE(RK, R, Y, ITASK, IER)

     if (IER == 2) then
        ! The root p=0 was found at R
        ! Y = Y(RK) not Y(R) so we have to find Y(R)
        call FCVODE(R, R, Y, ITASK, IER)
        if (IER == 0) then
           call store_solution(Y, hse_profile, K, R)
           KFIN = K
           write(*,*) 'p=0 Root Found! Integration halting.'
        else
           write(*,*) 'Error: p=0 Root find was unsuccessful!'
        end if
        exit
     else
        ! Store solution values if this was a successful integration
        call store_solution(Y, hse_profile, K, R)
        KFIN = K
     end if

  end do

  if (KFIN == NDR_SAVE) then
     write(*,*) 'WARNING: p=0 was not reached! The profile may be incomplete!'
  end if
  
  ! Deallocate CVODE Memory
  call FCVROOTFREE
  call FCVFREE

  ! Write profile to output file
  call write_solution(hse_profile, KFIN)
  
  ! Print output to console
  write(*,'(A,ES25.14)') 'Integrated to radius: ', R
  write(*,'(A,ES25.14)') 'Final Mass: ', Y(jmass)
  write(*,'(A,ES25.14)') 'Final Pressure: ', Y(jpres)
  
end program integrator

