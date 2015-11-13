module polytrope_eos
  use eos_types
  
  implicit none

  type :: poly_internal_vars
     double precision :: K, gamma
  end type poly_internal_vars

  type(eos_data), save :: poly_vars
  type(poly_internal_vars), save :: poly_ivars
  
contains

  subroutine poly_init()
    namelist /polyparams/ poly_ivars
    rewind(unit=parameter_file_unit)
    read(unit=parameter_file_unit, nml=polyparams)

    namelist /polyeosinit/ poly_vars ! Initialize the density and temperature
    rewind(unit=parameter_file_unit)
    read(unit=parameter_file_unit, nml=eosinitv)
  end subroutine poly_init

  subroutine poly_p(pres)
    double precision, intent(in) :: pres
    poly_vars%p = pres
    poly_vars%rho = (poly_vars%p/poly_ivars%K)**(1.0d0/poly_ivars%gamma)
    call poly_compute()
  end subroutine poly_p

  subroutine poly_rho(dens)
    double precision, intent(in) :: dens
    poly_vars%rho = dens
    poly_vars%p = poly_ivars%K * (poly_vars%rho**poly_ivars%gamma)
    call poly_compute()
  end subroutine poly_rho

  subroutine poly_compute()
    poly_vars%e = (poly_vars%p/poly_vars%rho)/(poly_ivars%gamma - 1.0d0)
    poly_vars%drho_dp = (1.0d0/(poly_ivars%K*poly_ivars%gamma))* &
         poly_vars%rho**(1.0d0-poly_ivars%gamma)
    poly_vars%de_dp = (-poly_vars%p/(poly_ivars%gamma - 1.0d0)/ &
         (poly_vars%rho**2))*poly_vars%drho_dp + &
         (1.0d0/poly_vars%rho/(poly_ivars%gamma - 1.0d0))
  end subroutine poly_compute
  
end module polytrope_eos
