module ls_wrap_eos
  ! Wrapper for the Lattimer-Swesty EOS
  use eos_types
  use physical_constants
  use parameters
  include 'lseos_v2.7.f'
  include 'eos_m4a.inc'
  
  implicit none

  type :: ls_internal_vars
     double precision :: ye      ! Electron fraction
     double precision :: t       ! Temperature (MeV)
     double precision :: p       ! Pressure (MeV/fm^3)
     !density of nucleons in nuclei, good initial guess 0.155
     double precision :: nuclei_nucleon_density 
     double precision :: eta_po  ! good initial guess -15
     double precision :: eta_no  ! good initial guess -10
     double precision :: baryon_density           ! Baryon density (1/fm^3)
     double precision :: exterior_proton_fraction ! Exterior proton fraction
     double precision, parameter :: XPREV=0.0d0   ! EOS dummy variable
     integer, parameter :: IFLAG =1  ! Sets EOS input variable
     integer, parameter :: EOSFLG=0  ! Returns which EOS was used internally
     integer, parameter :: FORFLG=0  ! Sets whether EOS should determine internal EOS to use
     integer, parameter :: SF    =1  ! EOS Error flag, 1 if successful
  end type ls_internal_vars

  type(eos_data), save :: ls_vars
  type(ls_internal_vars), save :: ls_ivars
  
contains

  subroutine ls_init()
    namelist /lsparams/ ls_ivars
    rewind(unit=parameter_file_unit)
    read(unit=parameter_file_unit, nml=lsparams)

    namelist /lseosinit/ ls_vars ! Initialize the density and temperature
    rewind(unit=parameter_file_unit)
    read(unit=parameter_file_unit, nml=eosinitv)
  end subroutine ls_init

  subroutine quant_cgs_to_nuc()
    ls_ivars%t = pc%kb_MeV_per_K * ls_vars%t
    ls_ivars%baryon_density = ls_vars%rho/pc%amu
    ls_ivars%p = ls_vars%p * pc%MeV_per_erg * pc%cm_per_fm**3
  end subroutine quant_cgs_to_nuc
  
  subroutine quant_nuc_to_cgs()
    ls_vars%t = ls_ivars%t/pc%kb_MeV_per_K
    ls_vars%rho = ls_ivars%baryon_density * pc%amu
    ls_vars%p = ls_ivars%p/(pc%MeV_per_erg * pc%cm_per_fm**3)
    ls_vars%e = UTOT/pc%amu
  end subroutine quant_nuc_to_cgs

  subroutine ls_p(pres)
    double precision, intent(in) :: pres
    
  end subroutine ls_p

  subroutine ls_rho(dens)
    double precision, intent(in) :: dens
  end subroutine ls_rho

  subroutine ls_compute()
    !DON
  end subroutine ls_compute
  
end module ls_wrap_eos
