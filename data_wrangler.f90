module data_wrangler

  use cvode_indices
  use polytrope_eos

  implicit none

contains

  subroutine store_solution(ysol, hse_profile, k, r)
    double precision, dimension(:), intent(in)   :: ysol 
    double precision, dimension(:,:), intent(inout) :: hse_profile
    integer, intent(in) :: k
    double precision, intent(in) :: r

    hse_profile(imass, k) = ysol(jmass)
    hse_profile(ipres, k) = ysol(jpres)

    call eos_p(ysol(jpres))
    hse_profile(idens, k) = eos_data%rho
    hse_profile(ieint, k) = eos_data%e
    hse_profile(i_ddens_dp, k) = eos_data%drho_dp
    hse_profile(i_deint_dp, k) = eos_data%de_dp

    hse_profile(irads, k) = r
  end subroutine store_solution
  
  subroutine write_solution(hse_profile, KFIN)
    double precision, dimension(:,:), intent(in) :: hse_profile
    integer, intent(in) :: KFIN
    integer :: K, J

    ! Output Profile Misc
    character(len=50) :: rfmt
    character(len=50) :: hfmt
    character(len=*), parameter  :: profile_file_name = 'hse_profile.dat'
    
    !! Save profile to file
    ! Set format for writing column entries
    write(hfmt,'(A,I5,A)') '(', NPROFILE, '(1x,A25))'
    write(rfmt,'(A,I5,A)') '(', NPROFILE, '(1x,ES25.14))'
    open(unit=10, file=profile_file_name, action='write', &
         status='replace', recl=(25*NPROFILE+10))
    write(10, fmt=hfmt) 'Mass', 'Pressure', 'Density', 'Einternal', 'D_RHO_D_P', 'D_Eint_D_P', 'Radius'
    do K = 1, KFIN
       write(10, fmt=rfmt) (hse_profile(J, K), J=1, NPROFILE)
    end do
    close(unit=10)
  end subroutine write_solution

end module data_wrangler
