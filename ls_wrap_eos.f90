module ls_wrap_eos
  ! Wrapper for the Lattimer-Swesty EOS
  use eos_types
  use physical_constants

  include 'lseos_v2.7.f'
  include 'eos_m4c.inc'
  
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

     ! Relative tolerance of pressure solution
     double precision :: pressure_solve_rtol = 1.0d-12
     integer :: pressure_solve_max_nsteps = 1000
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

    ! The initial density and temperature are central values, so find central pressure now
    call ls_rho(ls_vars%rho)
  end subroutine ls_init

  subroutine quant_cgs_to_nuc()
    ls_ivars%t = pc%kb_MeV_per_K * ls_vars%t
    ls_ivars%baryon_density = ls_vars%rho/pc%amu
    ls_ivars%p = ls_vars%p * pc%MeV_per_erg * pc%cm_per_fm**3
  end subroutine quant_cgs_to_nuc
  
  subroutine quant_nuc_to_cgs()
    double precision :: chi_pt, chi_pn
    double precision :: dp_dt
    chi_pt = (pc%cm_per_fm**(-3)) * pc%kb_MeV_per_K/pc%MeV_per_erg
    chi_pn = pc%amu * pc%MeV_per_erg
    ls_vars%t   = ls_ivars%t/pc%kb_MeV_per_K
    ls_vars%rho = ls_ivars%baryon_density * pc%amu
    ls_vars%p = ls_ivars%p/(pc%MeV_per_erg * pc%cm_per_fm**3)
    ls_vars%e = UTOT/pc%amu
    ls_vars%drho_dp = chi_pn/DPDN
    dp_dt = DPDT * chi_pt
    ls_vars%de_dp   = (-ls_vars%t * dp_dt + ls_vars%p) * ls_vars%drho_dp/(ls_vars%rho**2)
  end subroutine quant_nuc_to_cgs

  subroutine ls_p(p_target)
    double precision, intent(in) :: p_target
    double precision :: dp_drho_i, p_i, rho_i, f_i, df_drho_i
    integer :: i
    ! Find the density which will give this pres
    ! This subroutine assumes ls_compute() has been called before.
    ! First, get the pressure at the current density
    rho_i = ls_vars%rho
    p_i = ls_vars%p
    dp_drho_i = 1.0d0/ls_vars%drho_dp
    f_i = (p_i-p_target)/p_target
    df_drho_i = dp_drho_i/p_target
    
    ! Check to see if we're already at the solution (saves an EOS call)
    if ( abs(f_i) <= ls_ivars%pressure_solve_rtol ) then
       ! Current density is the solution
       return
    else
       rho_i = rho_i - f_i/df_drho_i
    end if

    ! Newton-Raphson
    do i = 1, ls_ivars%pressure_solve_max_nsteps
       call p_target_fun(p_i, p_target, f_i, df_drho_i, rho_i)
       ! Check to see if we're at the solution
       if ( abs(f_i) <= ls_ivars%pressure_solve_rtol ) then
          ! Current density is the solution
          return
       else
          rho_i = rho_i - f_i/df_drho_i
       end if
    end do

    ! Error: reached pressure_solve_max_nsteps
    print(*,*) 'ERROR: reached pressure_solve_max_nsteps'
    stop
  end subroutine ls_p

  subroutine p_target_fun(p, p_target, f, fp, rho)
    double precision, intent(out) :: p, f, fp
    double precision, intent(in)    :: p_target, rho
    call ls_rho(rho)
    p = ls_vars%p
    f = (p-p_target)/p_target
    fp = 1.0d0/ls_vars%drho_dp/p_target
  end subroutine p_target_fun

  subroutine ls_rho(dens)
    double precision, intent(in) :: dens
    ! Find the pressure given this density
    ls_vars%rho = dens
    call ls_compute()
  end subroutine ls_rho

  subroutine ls_compute()
    double precision, dimension(4) :: input_variables

    call quant_cgs_to_nuc()
    
    input_variables(1) = ls_ivars%t
    input_variables(2) = ls_ivars%nuclei_nucleon_density
    input_variables(3) = ls_ivars%eta_po
    input_variables(4) = ls_ivars%eta_no
    call INVEOS( &
         input_variables, &
         0.0d0, &
         ls_ivars%ye, &
         ls_ivars%baryon_density, &
         ls_ivars%IFLAG, &
         ls_ivars%EOSFLG, &
         ls_ivars%FORFLG, &
         ls_ivars%SF, &
         ls_ivars%XPREV, &
         ls_ivars%exterior_proton_fraction, &
         )
    
    ls_ivars%t                      = input_variables(1)
    ls_ivars%nuclei_nucleon_density = input_variables(2)
    ls_ivars%eta_po                 = input_variables(3)
    ls_ivars%eta_no                 = input_variables(4)
    
    ls_ivars%p = PTOT

    call quant_nuc_to_cgs()
  end subroutine ls_compute
  
end module ls_wrap_eos
