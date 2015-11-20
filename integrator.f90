program integrator
  use eos
  use cvode_indices
  use cvode_parameters
  use cvodehse
  use data_wrangler
  use parameters
  
  implicit none

  integer KOUNTER
  
  ! Initialize parameters
  call init_parameters()
  cv_pars%central_density = eos_vars%rho

  ! Allocate Profile array
  allocate( cv_data%hse_profile(NPROFILE, cv_pars%NDR_SAVE) )
  
  ! Print Physical Parameters
  write(*,*) 'Integrating starting with:'
  write(*,'(A,ES25.14)') 'R0 = ', cv_data%R0
  write(*,'(A,ES25.14)') 'Central Density = ', cv_pars%central_density
  
  ! Initialize the Integration
  ! Mass (r=0) = 0.0
  cv_data%Y0(jmass) = 0.0d0
  ! Pressure (r=0) = EOS(central_density)
  call eos_rho(cv_pars%central_density) ! Compute EOS quantities at density central_density
  cv_data%Y0(jpres) = eos_vars%p ! EOS Pressure for density central_density
  ! Integration target R_TARGET = NDR_SAVE*DR_SAVE
  cv_data%R_TARGET = DBLE(cv_pars%NDR_SAVE)*cv_pars%DR_SAVE

  ! Print Solution IC's
  write(*,*) 'Initial Conditions:'
  write(*,'(A,ES25.14)') 'Central Mass: ', cv_data%Y0(jmass)
  write(*,'(A,ES25.14)') 'Central Pressure: ', cv_data%Y0(jpres)
  
  ! Setup CVODE
  call FNVINITS(cv_data%KEY, cv_data%NEQ, cv_data%IER)
  call FCVMALLOC(cv_data%R0, cv_data%Y0, cv_data%METH, cv_data%ITMETH, &
       cv_data%IA_TOL, cv_pars%R_TOL, cv_pars%A_TOL, &
       cv_data%IOUT, cv_data%ROUT, cv_data%IPAR, cv_data%RPAR, cv_data%IER)
  if (cv_data%IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVMALLOC!'
     stop
  end if
  call FCVROOTINIT(cv_data%NRTFN, cv_data%IER)
  if (cv_data%IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVROOTINIT!'
     stop
  end if
  call FCVSETIIN('MAX_NSTEPS', cv_pars%MAX_NSTEPS, cv_data%IER)
  call FCVLAPACKDENSE(cv_data%NEQ, cv_data%IER)
  call FCVLAPACKDENSESETJAC(cv_pars%SET_JAC_FLAG, cv_data%IER)
  
  ! Do the integration
  do KOUNTER = 1, cv_pars%NDR_SAVE
     cv_data%RK = min(cv_data%R0 + DBLE(KOUNTER)*cv_pars%DR_SAVE, cv_data%R_TARGET)
     call FCVODE(cv_data%RK, cv_data%R, cv_data%Y, cv_pars%ITASK, cv_data%IER)

     if (cv_data%IER == 2) then
        ! The root p=0 was found at R
        ! Y = Y(RK) not Y(R) so we have to find Y(R)
        call FCVODE(cv_data%R, cv_data%R, cv_data%Y, cv_pars%ITASK, cv_data%IER)
        if (cv_data%IER == 0) then
           call store_solution(cv_data%Y, cv_data%hse_profile, KOUNTER, cv_data%R)
           cv_data%KFIN = KOUNTER
           write(*,*) 'p=0 Root Found! Integration halting.'
        else
           write(*,*) 'Error: p=0 Root find was unsuccessful!'
        end if
        exit
     else
        ! Store solution values if this was a successful integration
        call store_solution(cv_data%Y, cv_data%hse_profile, KOUNTER, cv_data%R)
        cv_data%KFIN = KOUNTER
     end if

  end do

  if (cv_data%KFIN == cv_pars%NDR_SAVE) then
     write(*,*) 'WARNING: p=0 was not reached! The profile may be incomplete!'
  end if
  
  ! Deallocate CVODE Memory
  call FCVROOTFREE
  call FCVFREE

  ! Write profile to output file
  call write_solution(cv_data%hse_profile, cv_data%KFIN)
  
  ! Print output to console
  write(*,'(A,ES25.14)') 'Integrated to radius: ', cv_data%R
  write(*,'(A,ES25.14)') 'Final Mass: ', cv_data%Y(jmass)
  write(*,'(A,ES25.14)') 'Final Pressure: ', cv_data%Y(jpres)
  
end program integrator
