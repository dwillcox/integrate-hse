module cvodehse

  implicit none

contains
  
  subroutine FCVFUN(R, Y_VEC, YDOT_VEC, IPAR, RPAR, IER) bind(C, name='fcvfun_')
    use cvode_indices
    use eos
    use physical_constants

    double precision, dimension(*), intent(in) :: Y_VEC, RPAR
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(out) :: YDOT_VEC
    double precision, intent(in) :: R
    integer, intent(out) :: IER

    double precision :: M, p, rho, e

    IER = -1
    
    M = Y_VEC(jmass)
    p = Y_VEC(jpres)

    call eos_p(p)
    rho = eos_vars%rho
    e = eos_vars%e

    if (R == 0.0d0) then
       YDOT_VEC(jmass) = 0.0d0
       YDOT_VEC(jpres) = 0.0d0
    else
       YDOT_VEC(jmass) = 4.0d0 * pc%PI * (R**2) * rho * (1.0d0 + e/(pc%c**2))
       YDOT_VEC(jpres) = -pc%G/(R**2)*(rho*(1.0d0 + e/(pc%c**2)) + p/(pc%c**2))* &
            (M + 4.0d0 * pc%PI * (R**3) * p/(pc%c**2))/(1.0d0 - pc%G * M/(R*(pc%c**2)))
    end if
    
    IER = 0 ! Successful
  end subroutine FCVFUN

  subroutine FCVDJAC(NEQ, R, Y_VEC, FY_VEC, DJAC, H_STEP, IPAR, RPAR, WK1, WK2, WK3, IER) bind(C, name='fcvdjac_')
    use cvode_indices
    use physical_constants
    use eos

    integer, intent(in) :: NEQ ! number of ODEs
    double precision, intent(in) :: R ! independent variable
    double precision, dimension(*), intent(in) :: Y_VEC, FY_VEC ! y and its derivative
    double precision, dimension(NEQ,*), intent(out) :: DJAC ! dense Jacobian
    double precision, intent(in) :: H_STEP ! current stepsize
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    double precision, dimension(NEQ), intent(in) :: WK1, WK2, WK3
    integer, intent(out) :: IER

    double precision :: drho_dp, de_dp, e, rho
    double precision :: M, p
    double precision :: chi, xi ! Scratch variables

    IER = -1
    
    M = Y_VEC(jmass)
    p = Y_VEC(jpres)
    
    call eos_p(p)

    ! K = eos_vars%K
    ! gamma = eos_vars%gamma
    drho_dp = eos_vars%drho_dp
    de_dp = eos_vars%de_dp
    e = eos_vars%e
    rho = eos_vars%rho

    if (R == 0.0d0) then
       DJAC(1,1) = 0.0d0
       DJAC(1,2) = 0.0d0
       DJAC(2,1) = 0.0d0
       DJAC(2,2) = 0.0d0
    else
       ! !! Compute the (Dense) Jacobian Matrix and store it column-wise in DJAC
       ! ! DJAC(1,1) = d(dM/dr)/dM
       ! DJAC(1,1) = 0.0d0

       ! ! DJAC(1,2) = d(dM/dr)/dp
       ! DJAC(1,2) = 4.0d0*(pc%PI)*(R**2)*(K**(-1.0d0/gamma) * &
       !      p**(-1.0d0 + 1.0d0/gamma) / gamma + &
       !      1.0d0/(gamma - 1.0d0)/(pc%c**2))

       ! chi = 1.0d0 - (pc%G)*M/(R*(pc%c**2))
       ! ! DJAC(2,1) = d(dp/dr)/dM
       ! DJAC(2,1) = -(pc%G)/(R**2) * (rho*(1.0d0 + e/(pc%c**2)) + p/(pc%c**2)) * &
       !      (1.0d0/chi + (M + 4.0d0*(pc%PI)*(R**3)*p/(pc%c**2))*(pc%G)/ &
       !      (r*(pc%c**2))/(chi**2))

       ! xi = (M + 4.0d0*(pc%PI)*(R**3)*p/(pc%c**2))/chi
       ! ! DJAC(2,2) = d(dp/dr)/dp
       ! DJAC(2,2) = -(pc%G)/(R**2) * ( &
       !      (1.0d0 + e/(pc%c**2)) * xi * drho_dp + &
       !      (rho/(pc%c**2)) * xi * de_dp + &
       !      (xi/(pc%c**2) + (rho*(1.0d0 + e/(pc%c**2)) + p/(pc%c**2))* &
       !      (4.0d0*(pc%PI)*R**3)/(pc%c**2)/chi) &
       !      )
    end if

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
