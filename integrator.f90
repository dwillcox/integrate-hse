program integrator
  use eos
  use cvode_indices
  use cvodehse
  use data_wrangler
  use parameters
  
  implicit none

  ! Initialize parameters
  call init_parameters()
  central_density = eos_vars%rho

  ! Allocate Profile array
  allocate( hse_profile(NPROFILE, NDR_SAVE) )
  
  ! Print Physical Parameters
  write(*,*) 'Integrating starting with:'
  write(*,'(A,ES25.14)') 'R0 = ', R0
  write(*,'(A,ES25.14)') 'Central Density = ', central_density
  
  ! Initialize the Integration
  ! Mass (r=0) = 0.0
  Y0(jmass) = 0.0d0
  ! Pressure (r=0) = EOS(central_density)
  call eos_rho(central_density) ! Compute EOS quantities at density central_density
  Y0(jpres) = eos_data%p ! EOS Pressure for density central_density
  ! Integration target R_TARGET = NDR_SAVE*DR_SAVE
  R_TARGET = DBLE(NDR_SAVE)*DR_SAVE

  ! Print Solution IC's
  write(*,*) 'Initial Conditions:'
  write(*,'(A,ES25.14)') 'Central Mass: ', Y0(jmass)
  write(*,'(A,ES25.14)') 'Central Pressure: ', Y0(jpres)
  
  ! Setup CVODE
  call FNVINITS(KEY, NEQ, IER)
  call FCVMALLOC(R0, Y0, METH, ITMETH, IA_TOL, R_TOL, A_TOL, &
       IOUT, ROUT, IPAR, RPAR, IER)
  if (IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVMALLOC!'
     stop
  end if
  call FCVROOTINIT(NRTFN, IER)
  if (IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVROOTINIT!'
     stop
  end if
  call FCVSETIIN('MAX_NSTEPS', MAX_NSTEPS, IER)
  call FCVLAPACKDENSE(NEQ, IER)
  call FCVLAPACKDENSESETJAC(SET_JAC_FLAG, IER)
  
  ! Do the integration
  do K = 1, NDR_SAVE
     RK = min(R0 + DBLE(K)*DR_SAVE, R_TARGET)
     call FCVODE(RK, R, Y, ITASK, IER)

     if (IER == 2) then
        ! The root p=0 was found at R
        ! Y = Y(RK) not Y(R) so we have to find Y(R)
        call FCVODE(R, R, Y, ITASK, IER)
        if (IER == 0) then
           call store_solution(Y, hse_profile, K, R)
           KFIN = K
           write(*,*) 'p=0 Root Found! Integration halting.'
        else
           write(*,*) 'Error: p=0 Root find was unsuccessful!'
        end if
        exit
     else
        ! Store solution values if this was a successful integration
        call store_solution(Y, hse_profile, K, R)
        KFIN = K
     end if

  end do

  if (KFIN == NDR_SAVE) then
     write(*,*) 'WARNING: p=0 was not reached! The profile may be incomplete!'
  end if
  
  ! Deallocate CVODE Memory
  call FCVROOTFREE
  call FCVFREE

  ! Write profile to output file
  call write_solution(hse_profile, KFIN)
  
  ! Print output to console
  write(*,'(A,ES25.14)') 'Integrated to radius: ', R
  write(*,'(A,ES25.14)') 'Final Mass: ', Y(jmass)
  write(*,'(A,ES25.14)') 'Final Pressure: ', Y(jpres)
  
end program integrator
