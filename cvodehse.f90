module cvodehse

  implicit none

contains

  subroutine INIT_PARAMETERS(IPAR, RPAR)
    use cvode_indices
    integer, dimension(*), intent(inout) :: IPAR
    double precision, dimension(*), intent(inout) :: RPAR
  end subroutine INIT_PARAMETERS
  
  subroutine FCVFUN(R, Y, YDOT, IPAR, RPAR, IER) bind(C, name='fcvfun_')
    use cvode_indices
    use polytrope_eos
    use physical_constants

    double precision, dimension(*), intent(in) :: Y, RPAR
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(out) :: YDOT
    double precision, intent(in) :: R
    integer, intent(out) :: IER

    double precision :: M, p, rho, e

    IER = -1
    
    M = Y(jmass)
    p = Y(jpres)

    call eos_p(p)
    rho = eos_data%rho
    e = eos_data%e
    
    YDOT(jmass) = 4.0d0*PI*(R**2)*rho*(1.0d0 + e/c**2)
    YDOT(jpres) = -G/(R**2)*(rho*(1.0d0 + e/c**2) + p/c**2)* &
         (M + 4.0d0*PI*(R**3)*p/c**2)/(1.0d0 - G*M/(R*c**2))

    IER = 0 ! Successful
  end subroutine FCVFUN

  subroutine FCVDJAC(NEQ, R, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER) bind(C, name='fcvdjac_')
    use cvode_indices
    use physical_constants
    use polytrope_eos

    integer, intent(in) :: NEQ ! number of ODEs
    double precision, intent(in) :: R ! independent variable
    double precision, dimension(*), intent(in) :: Y, FY ! y and its derivative
    double precision, dimension(NEQ,*), intent(out) :: DJAC ! dense Jacobian
    double precision, intent(in) :: H ! current stepsize
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    double precision, dimension(NEQ), intent(in) :: WK1, WK2, WK3
    integer, intent(out) :: IER

    double precision :: K, gamma, drho_dp, de_dp, e, rho
    double precision :: M, p
    double precision :: chi, xi ! Scratch variables

    IER = -1
    
    M = Y(jmass)
    p = Y(jpres)
    
    call eos_p(p)

    K = eos_data%K
    gamma = eos_data%gamma
    drho_dp = eos_data%drho_dp
    de_dp = eos_data%de_dp
    e = eos_data%e
    rho = eos_data%rho
    
    !! Compute the (Dense) Jacobian Matrix and store it column-wise in DJAC
    ! DJAC(1,1) = d(dM/dr)/dM
    DJAC(1,1) = 0.0d0
    
    ! DJAC(1,2) = d(dM/dr)/dp
    DJAC(1,2) = 4.0d0*PI*(R**2)*(K**(-1.0d0/gamma) * p**(-1.0d0 + 1.0d0/gamma) / gamma + &
         1.0d0/(gamma - 1.0d0)/c**2)

    chi = 1.0d0 - G*M/(R*c**2)
    ! DJAC(2,1) = d(dp/dr)/dM
    DJAC(2,1) = -G/(R**2) * (rho*(1.0d0 + e/c**2) + p/c**2) * &
         (1.0d0/chi + (M + 4.0d0*PI*(R**3)*p/c**2)*G/(r*c**2)/(chi**2))

    xi = (M + 4.0d0*PI*(R**3)*p/c**2)/chi
    ! DJAC(2,2) = d(dp/dr)/dp
    DJAC(2,2) = -G/(R**2) * ( &
         (1.0d0 + e/c**2) * xi * drho_dp + &
         (rho/c**2) * xi * de_dp + &
         (xi/c**2 + (rho*(1.0d0 + e/c**2) + p/c**2)*(4.0d0*PI*R**3)/(c**2)/chi) &
         )

    IER = 0 ! Success
  end subroutine FCVDJAC

  subroutine FCVROOTFN(R, Y, G, IPAR, RPAR, IER) bind(C, name='fcvrootfn_')
    use cvode_indices

    double precision, intent(in) :: R
    double precision, dimension(*), intent(in) :: Y
    double precision, dimension(*), intent(inout) :: G
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    integer, intent(out) :: IER

    IER = -1
    G(iproot) = Y(jpres)
    IER = 0
  end subroutine FCVROOTFN

end module cvodehse
