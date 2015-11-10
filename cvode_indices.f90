module CVODE_INDICES
  implicit none

  ! Indices for solution Y
  integer, parameter :: jmass = 1
  integer, parameter :: jpres = 2

    ! Indices for hse_profile
  integer, parameter :: NPROFILE = 5
  integer, parameter :: imass = 1
  integer, parameter :: ipres = 2
  integer, parameter :: idens = 3
  integer, parameter :: ieint = 4
  integer, parameter :: irads = 5
end module CVODE_INDICES
