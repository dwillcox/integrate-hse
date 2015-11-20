module CVODE_PARAMETERS

  implicit none

  type :: cv_pars_t
     double precision :: central_density
     double precision :: DR_SAVE ! Save solution each DR_SAVE (cm)
     integer :: NDR_SAVE ! (Such that NDR_SAVE*DR_SAVE=R_TARGET=30d5)
     double precision :: R_TOL = 1.0d-12 ! Relative tolerance
     double precision :: A_TOL = 1.0d-12 ! Absolute tolerance
     integer :: MAX_NSTEPS = 10000 ! Maximum number of steps to solution
     integer :: SET_JAC_FLAG = 1 ! Use my supplied Jacobian
     integer :: ITASK = 1 ! Normal mode (overshoot and interpolate)
  end type cv_pars_t
  
  type :: cv_data_t
     ! Integration Data Structures and Parameters
     double precision :: R0 = 0.0d0
     double precision, dimension(2) :: Y0
     double precision :: R_TARGET ! Try to integrate to RFIN = NDRV*DRSV
     double precision :: RK ! Stores the radius at the current interval
!     integer :: K ! Stores which interval we're in (count to NDR_SAVE)
     integer :: KFIN ! Stores how many intervals we went (K = 1 to KFIN)
     double precision               :: R ! Holds the R integrated to
     double precision, dimension(2) :: Y ! Holds the solution Y integrated to
     double precision, dimension(2) :: DKY ! Holds the p=0 root solution for Y
     double precision, dimension(:,:), allocatable :: hse_profile

     ! CVODE
     integer :: KEY = 1 ! Use CVODE
     ! NOTE: Below, integer*8 necessary to match the C int type expected
     integer*8 :: NEQ = 2 ! Size of ODE system 
     integer ::   IER     ! Error flag
     integer ::   METH = 1 ! Use non-stiff Adams integration
     integer :: ITMETH = 2 ! Use Newton iteration
     integer ::  IA_TOL = 1 ! Use scalar absolute tolerance
     integer, dimension(21) ::         IOUT ! Optional integer outputs
     double precision, dimension(6) :: ROUT ! Optional real outputs
     integer :: NRTFN = 1 ! Number of root functions to solve during CVODE
     ! User integer parameters
     integer, dimension(1) :: IPAR = (/ 0 /)
     ! User real parameters
     double precision, dimension(1) :: RPAR = (/ 0.0d0 /) 
  end type cv_data_t

  type(cv_data_t), save :: cv_data
  type(cv_pars_t), save :: cv_pars

contains
  subroutine cvode_init(pfile_unit)
    integer, intent(in) :: pfile_unit
    ! Read runtime parameters
    namelist /intparams/ cv_pars
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=intparams)
  end subroutine cvode_init

end module CVODE_PARAMETERS
