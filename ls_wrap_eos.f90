module ls_wrap_eos
  ! Wrapper for the Lattimer-Swesty EOS
  use eos_m4c_data
  use lattimer_swesty_eos
  use eos_types
  use physical_constants

!  include 'lseos_v2.7.f'
!  include 'eos_m4c.inc'
  
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
     double precision :: XPREV=0.0d0   ! EOS dummy variable
     integer :: IFLAG =1  ! Sets EOS input variable
     integer :: EOSFLG=0  ! Returns which EOS was used internally
     integer :: FORFLG=0  ! Sets whether EOS should determine internal EOS to use
     integer :: SF    =1  ! EOS Error flag, 1 if successful

     ! Relative tolerance of pressure solution
     double precision :: pressure_solve_rtol = 1.0d-12
     integer :: pressure_solve_max_nsteps = 1000
  end type ls_internal_vars

  type(eos_data), target, save :: ls_vars
  type(ls_internal_vars), save :: ls_ivars

  ! interface
  !    subroutine INVEOS( &
  !        input_vars, &
  !        told, &
  !        ye, &
  !        baryon_density, &
  !        IFLAG, &
  !        EOSFLG, &
  !        FORFLG, &
  !        SF, &
  !        XPREV, &
  !        exterior_proton_fraction)
  !      double precision, dimension(4), intent(inout) :: input_vars
  !      double precision, intent(in) :: told, ye, baryon_density, XPREV
  !      double precision, intent(inout) :: exterior_proton_fraction
  !      integer, intent(in) :: IFLAG, FORFLG
  !      integer, intent(inout) :: EOSFLG, SF
  !    end subroutine INVEOS
  ! end interface

contains

  subroutine ls_init(pfile_unit)
    integer, intent(in) :: pfile_unit
    
    namelist /lsparams/ ls_ivars
    namelist /lseosinit/ ls_vars ! Initialize the density and temperature
    
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=lsparams)

    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=lseosinit)

    ! Initialize the Lattimer-Swesty EOS Module Variables
    call LOADMX()
    
    ! The initial density and temperature are central values, so find central pressure now
    call ls_rho(ls_vars%rho)
  end subroutine ls_init

  subroutine quant_cgs_to_nuc()
    ls_ivars%t = pc%kb_MeV_per_K * ls_vars%t
    write(*,*) 'Using density: ', ls_vars%rho
    ls_ivars%baryon_density = (ls_vars%rho/pc%amu)*pc%cm_per_fm**3
    write(*,*) 'Calculated baryon density: ', ls_ivars%baryon_density
    ls_ivars%p = ls_vars%p * pc%cm_per_fm**3/pc%erg_per_MeV
  end subroutine quant_cgs_to_nuc
  
  subroutine quant_nuc_to_cgs()
    double precision :: chi_pt, chi_pn
    double precision :: dp_dt
    chi_pt = (pc%cm_per_fm**(-3)) * pc%kb_MeV_per_K * pc%erg_per_MeV
    chi_pn = pc%amu/pc%erg_per_MeV
    ls_vars%t   = ls_ivars%t/pc%kb_MeV_per_K
    ls_vars%rho = ls_ivars%baryon_density*pc%amu/(pc%cm_per_fm**3)
    ls_vars%p = ls_ivars%p*pc%erg_per_MeV/(pc%cm_per_fm**3)
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

    write(*,*) '!!!!!!!!!!!!Now in ls_p... p_target: ',p_target
    write(*,*) 'p_i: ',p_i
    
    ! Check to see if we're already at the solution (saves an EOS call)
    if ( abs(f_i) <= ls_ivars%pressure_solve_rtol ) then
       ! Current density is the solution
       write(*,*) 'Returning current density'
       return
    else
       rho_i = rho_i - f_i/df_drho_i
    end if

    ! Newton-Raphson
    do i = 1, ls_ivars%pressure_solve_max_nsteps
       write(*,*) 'nr iteration: ',i
       write(*,*) 'rho_i: ',rho_i
       call p_target_fun(p_i, p_target, f_i, df_drho_i, rho_i)
       write(*,*) 'got back'
       ! Check to see if we're at the solution
       if ( abs(f_i) <= ls_ivars%pressure_solve_rtol ) then
          ! Current density is the solution
          return
       else
          write(*,*) 'rho_i prev: ',rho_i
          rho_i = rho_i - f_i/df_drho_i
          write(*,*) 'rho_i next: ',rho_i
          write(*,*) 'f_i: ', f_i
          write(*,*) 'dfdrho_i: ',df_drho_i
       end if
    end do

    ! Error: reached pressure_solve_max_nsteps
    write(*,*) 'ERROR: reached pressure_solve_max_nsteps'
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
    write(*,*) 'Density: ',dens
    call ls_compute()
  end subroutine ls_rho

  subroutine ls_compute()
    double precision, dimension(4) :: input_variables
    double precision :: t_old_pass

    call quant_cgs_to_nuc()

    t_old_pass = ls_ivars%t
    input_variables(1) = ls_ivars%t
    input_variables(2) = ls_ivars%nuclei_nucleon_density
    input_variables(3) = ls_ivars%eta_po
    input_variables(4) = ls_ivars%eta_no

    write(*,*) 'Calling INVEOS'
    write(*,*) 'input_variables(1): ',input_variables(1)
    write(*,*) 'input_variables(2): ',input_variables(2)
    write(*,*) 'input_variables(3): ',input_variables(3)
    write(*,*) 'input_variables(4): ',input_variables(4)
    write(*,*) 't_old_pass        : ',t_old_pass
    write(*,*) 'ye        : ',ls_ivars%ye
    write(*,*) 'baryon density        : ',ls_ivars%baryon_density
    write(*,*) 'IFLAG        : ',ls_ivars%IFLAG
    write(*,*) 'EOSFLG        : ',ls_ivars%EOSFLG
    write(*,*) 'FORFLG        : ',ls_ivars%FORFLG
    write(*,*) 'SF        : ',ls_ivars%SF
    write(*,*) 'XPREV        : ',ls_ivars%XPREV
    write(*,*) 'exterior_proton_fraction: ',ls_ivars%exterior_proton_fraction
    call INVEOS( &
         input_variables, &
         t_old_pass, &
         ls_ivars%ye, &
         ls_ivars%baryon_density, &
         ls_ivars%IFLAG, &
         ls_ivars%EOSFLG, &
         ls_ivars%FORFLG, &
         ls_ivars%SF, &
         ls_ivars%XPREV, &
         ls_ivars%exterior_proton_fraction)
    write(*,*) 'Returned from INVEOS'
    ls_ivars%t                      = input_variables(1)
    ls_ivars%nuclei_nucleon_density = input_variables(2)
    ls_ivars%eta_po                 = input_variables(3)
    ls_ivars%eta_no                 = input_variables(4)
!    write(*,*) 'Got back eta_po: ',input_variables(3)
!    write(*,*) 'Got back eta_no: ',input_variables(4)
    ls_ivars%p = PTOT

    call quant_nuc_to_cgs()
  end subroutine ls_compute
  
end module ls_wrap_eos
