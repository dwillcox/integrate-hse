module eos
  use eos_types
  use ls_wrap_eos
  use polytrope_eos

  implicit none

  procedure(sub), pointer :: eos_set_pars, eos_p, eos_rho, eos_compute
  type(eos_data), pointer :: eos_vars

contains
  subroutine eos_init(WHICH_EOS)
    integer, intent(in) :: WHICH_EOS

    if (WHICH_EOS == 1) then
       ! Polytropic EOS
       eos_vars => poly_vars
       eos_set_pars => poly_init
       eos_p => poly_p
       eos_rho => poly_rho
       eos_compute => poly_compute
    else if (WHICH_EOS == 2) then
       ! Lattimer-Swesty EOS
       eos_vars => ls_vars
       eos_set_pars => ls_init
       eos_p => ls_p
       eos_rho => ls_rho
       eos_compute => ls_compute
    end if
    call eos_set_pars()
  end subroutine eos_init
  
end module eos
