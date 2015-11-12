module physical_constants

  implicit none

  type :: physical_constants_t
     double precision, parameter :: PI = 3.141592653589793238462643383279502884197d0
     double precision, parameter :: cm_per_fm = 1.0d-13 

     ! From PDG 2014
     double precision, parameter :: c  = 299792458d2 ! cm/s
     double precision, parameter :: kb_MeV_per_K = 8.6173324d-11 ! MeV/K
     double precision, parameter :: G  = 6.67384d-8 ! cm^3/g/s^2
     double precision, parameter :: MeV_per_erg = 1.602176565d-6

     ! From AME 2012 Atomic Mass Evaluation
     double precision, parameter :: amu = 1660538.921d-30 ! g
  end type physical_constants_t

  type(physical_constants_t), save :: pc

end module physical_constants
