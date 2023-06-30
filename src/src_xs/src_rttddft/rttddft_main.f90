! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created Apr 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module is the kernel of a RT-TDDFT calculation.
!> It contains the subroutine `coordinate_rttddft_calculation`, which manages
!> a RT-TDDFT calculation. Also here are implemented the following subroutines:
!> print_jpa, print_nexc, print_total_energy, printTiming, uppot, uprho.
module rttddft_main
  use asserts, only: assert
  use errors_warnings, only: terminate_if_false
  use mod_charge_and_moment, only: chgval
  use mod_kpoint, only: nkpt
  use mod_lattice, only: omega
  use mod_mpi_env, only: mpiinfo
  use modinput, only: input, input_type
  use modmpi, only: rank, procs, ierr, mpi_env_k, distribute_loop, mpiglobal, barrier
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_CurrentDensity, only: UpdateCurrentDensity
  use rttddft_GlobalVariables
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_Energy, only: TotalEnergy, obtain_energy_rttddft
  use rttddft_io, only: open_file_info, write_file_info, write_file_info_header, &
    close_file_info, write_wavefunction, open_files_jpa, close_files_jpa, &
    write_jpa, open_file_etot, close_file_etot, write_total_energy, &
    open_file_nexc, close_file_nexc, write_nexc, open_file_timing, write_timing, close_file_timing
  use rttddft_VectorPotential, only: Calculate_Vector_Potential, Evolve_A_ind => Solve_ODE_Vector_Potential
  use rttddft_Wavefunction, only: UpdateWavefunction

  implicit none

  private
  integer                 :: iprint
  real(dp)                :: timei, timef

  type(TimingRTTDDFT),&
              allocatable :: timingstore(:)

  public :: coordinate_rttddft_calculation

contains

  !> This subroutine manages a RT-TDDFT calculation.
  !> <ol>
  !> <li> Run a single-shot groundstate calculation using the already converged
  !> density and potential. </li>
  !> <li> Obtain the KS wavefunctions, the density and the hamiltonian at 
  !> \( t=0 \). </li>
  !> <li> Evolve the wavefunctions, the density and the hamiltonian using the 
  !> desired time step. </li>
  !> </ol>
  subroutine coordinate_rttddft_calculation()
    use rttddft_init, only: initialize_rttddft
    use m_getunit, only: getunit
    use rttddft_NumberExcitations, only: Obtain_number_excitations
    use rttddft_screenshot, only: screenshot
    use mod_misc, only: filext

    ! counter for the number of iterations of real-time
    integer                 :: it
    ! prints data every nprint steps of the counter "it"
    integer                 :: nprint
    ! indexes of the first and the last k-points
    integer                 :: first_kpt, last_kpt
    integer                 :: niter_screenshot

    logical                 :: predCorrReachedMaxSteps

    character(50)           :: string
    character(len=100)      :: fieldType
    character(len=100)      :: vectorPotentialSolver

    real(dp)                :: timesave, timeiter
    real(dp)                :: timehml, timerest
    real(dp),allocatable    :: nex(:), ngs(:), nt(:)
    real(dp)                :: aindsave(3),pvecsave(3)
    real(dp)                :: jindsave(3),aextsave(3),atotsave(3)

    ! Variables to store data and print
    real(dp),allocatable    :: timestore(:),aindstore(:,:),atotstore(:,:)
    real(dp),allocatable    :: jindstore(:,:), pvecstore(:,:)
    logical                 :: take_screenshots
    logical, allocatable    :: screenshot_was_taken(:)

    type(TotalEnergy),allocatable &
                            :: etotstore(:)


    call timesec(timesave)

    ! Interface with input parameters
    nprint = input%xs%realTimeTDDFT%printAfterIterations
    fieldType = input%xs%realTimeTDDFT%laser%fieldType
    vectorPotentialSolver = input%xs%realTimeTDDFT%vectorPotentialSolver
    take_screenshots = associated( input%xs%realTimeTDDFT%screenshots )
    if( take_screenshots ) niter_screenshot = input%xs%realTimeTDDFT%screenshots%niter

    call sanity_checks( input, nprint, mpiglobal )
    
    ! Outputs general info to RTTDDFT_INFO.OUT and
    ! opens TIMING_RTTDDFT.OUT (if this is the case)
    if(input%xs%realTimeTDDFT%printTimingGeneral) call timesec (timei)
    if( rank == 0 ) then
      call open_file_info
      call write_file_info_header
    end if

    ! Allocate variables to be stored and printed only after nprint steps
    allocate(timestore(nprint))
    allocate(aindstore(3,nprint))
    allocate(atotstore(3,nprint))
    allocate(jindstore(3,nprint))
    allocate(pvecstore(3,nprint))

    call initialize_rttddft


    ! Print the vector potential, the current density and the polarization vector at t=0
    if( rank == 0 ) then
      call open_files_jpa
      call write_jpa( time, aind, atot, label='avec' )
      call write_jpa( time, pvec, label='pvec' )
      call write_jpa( time, jind, label='jind' )
    end if

    ! Initialize integers that contain the first and last k-point
    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    ! Total energy, number of excitations, and screenshots at t=0
    ! Total energy
    if ( calculateTotalEnergy ) then
      allocate(etotstore(nprint))
      call potcoul
      call potxc
      call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_gnd, etotstore(1) )
      ! Trick: we need an array to call the subroutine print_total_energy
      timestore(1) = time
      if( rank == 0 ) then
        call open_file_etot
        call write_total_energy( .True., 1, timestore(1), etotstore(1) )
      end if
    end if

    ! Number of excitations
    if (calculateNexc) then
      allocate(nex(nprint),ngs(nprint),nt(nprint))
      call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
        & evecfv_time, overlap, nex(1), ngs(1), nt(1) )
      ! Trick: we need an array to call the subroutine print_nexc
      timestore(1) = time
      if( rank == 0 ) then
        call open_file_nexc
        call write_nexc( .True., 1, timestore(1), nex(1), ngs(1), nt(1) )
      end if
    end if

    if ( take_screenshots ) call screenshot( 0, first_kpt, last_kpt, overlap, evecfv_gnd, &
        & evecfv_time, ham_time )

    if( printTimesGeneral ) then
      allocate( timingstore(nprint) )
      allocate( screenshot_was_taken(nprint), source=.False. )
      call timesec( timef )
      if( rank == 0 ) then 
        call open_file_timing
        call write_timing( timef-timesave ) !write time for initialization
      end if
    end if

    iprint = 1
    timeiter = timef
    ! This is the most important loop (performed for each time step \(\Delta t\)
    do it = 1, nsteps
      ! Variable to store the timing of each iteration
      timei = timeiter

      ! The "real time" t of our evolution
      time = time + tstep

      ! WAVEFUNCTION
      if ( predictorCorrector ) evecfv_save(:,:,:) = evecfv_time(:,:,:)
      call UpdateWavefunction( .False. )
      if ( printTimesGeneral ) call timesecRTTDDFT(timei,timef,timingstore(iprint)%t_wvf)

      ! Update the paramagnetic component of the induced current density
      call UpdateCurrentDensity( first_kpt, last_kpt, evecfv_time(:,:,:),jparanext(:) )
      if ( input%xs%realTimeTDDFT%subtractJ0 ) jparanext(:) = jparanext(:)-jparaspurious(:)
      if ( printTimesGeneral ) call timesecRTTDDFT(timei,timef,timingstore(iprint)%t_curr)

      ! DENSITY
      call uprho(it, .True.)

      ! KS-POTENTIAL
      call uppot(.True.)

      ! VECTOR POTENTIAL
      if(printTimesGeneral) call timesec(timei)
      ! Check if we need to save aind, pvec, atot and aext
      if( fieldType == 'external' .and. predictorCorrector .and. (vectorPotentialSolver /= 'euler') ) then
        aindsave(:) = aind(:)
        pvecsave(:) = pvec(:)
        atotsave(:) = atot(:)
        aextsave(:) = aext(:)
      end if
      if( fieldType == 'total' ) then
        call update_vector_potential( time, atot )
      else 
        call update_vector_potential( time, atot, aind, aext )
      end if
      if(printTimesGeneral) call timesecRTTDDFT(timei,timef,timingstore(iprint)%t_obtaina)

      ! INDUCED CURRENT
      ! Update the diamagnetic component of the induced current density
      jdia(:) = - atot(:) * chgval / c / omega
      ! Update the paramagnetic component of the induced current density
      jparaold(:) = jpara(:)
      jpara(:) = jparanext(:)
      ! Update the total induced current
      if ( predictorCorrector .and. (vectorPotentialSolver /= 'euler') ) then
        jindsave(:) = jind(:)
      end if
      jind(:) = jpara(:) + jdia(:)

      ! HAMILTONIAN
      call UpdateHam( .False., .False., &
          & printTimesGeneral, printTimesDetailed,&
          & timei, timef, timehml, timerest )
      if( printTimesGeneral ) then
        timingstore(iprint)%t_upham = timef - timei
        if ( printTimesDetailed ) then
          timingstore(iprint)%t_hmlint = timehml
          timingstore(iprint)%t_ham = timerest
        end if
        timei = timef
      end if

      ! Remark: it makes no sense to employ the predictor-corrector method with SE or EH!
      if ( predictorCorrector .and. (method /= 'SE') .and. (method /= 'EH') ) then
        call loopPredictorCorrector( it, maxstepsPredictorCorrector, &
          (fieldType == 'external') .and. (vectorPotentialSolver /= 'euler'), &
          first_kpt, last_kpt, aindsave, atotsave, aextsave, pvecsave, jindsave, &
          jparaold, predCorrReachedMaxSteps )
        
        if ( predCorrReachedMaxSteps .and. rank == 0 ) &
          write(*,*) 'Problems with convergence (PredCorr), time: ', time
        if (printTimesGeneral) call timesecRTTDDFT(timei,timef,timingstore(iprint)%t_predcorr)
      end if !predictor-corrector

      ! Obtain the total energy, if requested
      if( calculateTotalEnergy ) then
        call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_time, etotstore(iprint) )
        if ( printTimesDetailed ) call timesecRTTDDFT( timei, timef, timingstore(iprint)%t_toten )
      end if

      ! Obtain the number of excited electrons, if requested
      if( calculateNexc ) then
        call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
          & evecfv_time, overlap, nex(iprint), ngs(iprint), nt(iprint))
        if( printTimesDetailed ) call timesecRTTDDFT( timei, timef, timingstore(iprint)%t_nexc )
      end if

      ! Check if a screenshot has been requestedprint_jpa
      if ( take_screenshots ) then
        if ( mod( it, niter_screenshot ) == 0 ) then
          if( printTimesGeneral ) screenshot_was_taken(iprint) = .True.
          if( printTimesGeneral ) call timesec(timei)
          call screenshot( it, first_kpt, last_kpt, overlap, evecfv_gnd, &
            & evecfv_time, ham_time )
          if( printTimesGeneral ) call timesecRTTDDFT( timei, timef, timingstore(iprint)%t_screenshot )
        else 
          if( printTimesGeneral ) screenshot_was_taken(iprint) = .False.
        end if
      end if

      ! Store relevant information from this iteration
      timestore(iprint) = time
      aindstore(:, iprint) = aind(:)
      atotstore(:, iprint) = atot(:)
      pvecstore(:, iprint) = pvec(:)
      jindstore(:, iprint) = jind(:)

      ! Print relevant information, every 'nprint' steps
      if ( iprint == nprint ) then
        ! Update the counter
        iprint = 1
        if( rank == 0 ) then
          call write_jpa( timestore, aindstore, atotstore, label='avec' )
          call write_jpa( timestore, pvecstore, label='pvec' )
          call write_jpa( timestore, jindstore, label='jind' )
          if ( calculateTotalEnergy ) call write_total_energy( .False., nprint, &
            timestore(:), etotstore(:) )
          if ( calculateNexc ) call write_nexc( .False., nprint, timestore(:), &
            nex(:), ngs(:), nt(:) )
          if( printTimesGeneral ) then
            call timesecRTTDDFT( timeiter, timef, timingstore(nprint)%t_iteration )
            call write_timing( it, nprint, timingstore, screenshot_was_taken )
          end if
        end if
      else ! if ( iprint == nprint )
        if(printTimesGeneral) call timesecRTTDDFT(timeiter,timef,timingstore(iprint)%t_iteration)
        iprint = iprint + 1
      end if ! if ( iprint == nprint ) 
      ! Make all the processes wait here: the master alone has been writing the files above
      call barrier( mpi_env_k )
    end do ! do it = 1, nsteps

    if ( rank == 0 ) then
      call close_files_jpa
      call write_file_info( 'Real-time TDDFT calculation finished' )
      call close_file_info
      if( calculateTotalEnergy ) call close_file_etot
      if( calculateNexc ) call close_file_nexc
      if( printTimesGeneral ) call close_file_timing
    end if

    ! write wavefunction, and potential and density with _RTTDDFT.OUT as suffix
    string = filext
    filext = '_RTTDDFT'//trim(filext)
    call write_wavefunction( first_kpt, evecfv_time )
    if ( rank == 0 ) call writestate
    filext = string

    call deallocate_global_arrays

  end subroutine coordinate_rttddft_calculation

  !> This is just an interface to call the subroutine `[[UpdateDensity]]`, which
  !> updates the charge density
  subroutine uprho( iteration_counter, computeTiming )
    use rttddft_Density, only: UpdateDensity
    !> Tells how many time steps have already been executed
    integer, intent(in) :: iteration_counter
    !> Tells if we want to measure/store timings
    logical, intent(in) :: computeTiming


    if ( computeTiming .and. printTimesGeneral ) then
      if ( printTimesDetailed ) then
        call UpdateDensity( iteration_counter, timei, timef, timingstore(iprint)%t_dens_rho, &
          & timingstore(iprint)%t_dens_symrf, timingstore(iprint)%t_dens_rfmtctof, &
          & timingstore(iprint)%t_dens_addrhocr, timingstore(iprint)%t_dens_charge, &
          & timingstore(iprint)%t_dens_rhonorm )
      else
        call UpdateDensity( iteration_counter, timei, timef )
      end if
    else
      call UpdateDensity( iteration_counter )
    endif
    if ( computeTiming .and. printTimesGeneral ) then
      timingstore(iprint)%t_dens = timef-timei
      timei = timef
    endif
  end subroutine uprho


  !> This is just an interface to call the subroutines that updates the KS potential
  subroutine uppot(computeTiming)
    !> Tells if we want to measure/store timings
    logical, intent(in) :: computeTiming

    real(8)             :: timeisave

    timeisave = timei

    ! Compute the effective potential (with the updated density)
    call poteff
    if( computeTiming .and. printTimesDetailed ) call timesecRTTDDFT(timei,timef,&
      & timingstore(iprint)%t_poteff)
    ! Fourier transform effective potential to G-space
    call genveffig
    if( computeTiming .and. printTimesDetailed ) call timesecRTTDDFT(timei,timef,&
      timingstore(iprint)%t_genveffig)
    call genmeffig
    if(computeTiming .and. printTimesGeneral ) then
      call timesecRTTDDFT(timeisave,timef,timingstore(iprint)%t_uppot)
      if(printTimesDetailed) timingstore(iprint)%t_genmeffig = timef-timei
      timei = timef
    end if
  end subroutine uppot

  !> (private subroutine) Check if variables given in the input file make sense
  subroutine sanity_checks( inp, nprint, mpi_env  )
    !> type with the variables given in the input file
    type(input_type):: inp
    !> Variable that tells after how many iterations the output should be written
    integer(i32), intent(in) :: nprint
    !> Variable that encapsulates the MPI communicator that will abort if a sanity test fails
    type(mpiinfo) :: mpi_env
    
    call terminate_if_false( mpi_env, nprint > 0, &
      & 'Invalid non-positive number as printAfterIterations' )

    call terminate_if_false( mpi_env, .not. inp%groundstate%solver%packedmatrixstorage, &
      & 'Error: RT-TDDFT does not work with matrices stored in a packed form.' )

    ! Consistency check: check if no spin polarized calculations are requested.
    call terminate_if_false( mpi_env, .not. inp%groundstate%tevecsv, &
      & 'Error: only spin unpolarised calculations are possible with RT-TDDFT now.' )

    ! Consistency check: laser has been defined?
    call terminate_if_false( mpi_env, associated(inp%xs%realTimeTDDFT%laser), &
      & 'Element <laser> in <realTimeTDDFT> not found')

  end subroutine

  subroutine update_vector_potential( t, A_tot, A_ind, A_ext )
    !> time \(t\)
    real(dp), intent(in)              :: t
    !> total vector potential: \(A_{tot} = A_{ind} + A_{ext}\)
    real(dp), intent(inout)           :: A_tot(3)
    !> induced vector potential
    real(dp), intent(inout), optional :: A_ind(3)
    !> external vector potential
    real(dp), intent(inout), optional :: A_ext(3)

    logical :: all_fields_present

    all_fields_present = present(A_ind)
    if( all_fields_present ) call assert( present(A_ext), 'If A_ind is present, then A_ext must also be' )

    !TODO(Ronaldo): Refactor `Evolve_A_ind` to avoid globals
    call Evolve_A_ind( .not. all_fields_present )
    if( all_fields_present ) then
      call Calculate_Vector_Potential( t, A_ext(:) )
      A_tot(:) = A_ind(:) + A_ext(:)
    else
      call Calculate_Vector_Potential( t, A_tot(:) )
    end if
  end subroutine

  subroutine loopPredictorCorrector( it, maxSteps, evolveA, first_kpt, last_kpt, &
    aindsave, atotsave, aextsave, pvecsave, jindsave, jparasave, maxStepsReached )
    !> current iteration number in the RT-TDDFT loop
    integer(i32), intent(in)       :: it
    !> maximum steps in the predictor corrector scheme
    integer(i32), intent(in)       :: maxSteps
    !> If `.True.`, the vector potential is evolved in each step
    logical, intent(in)            :: evolveA
    !> index of the first `k-point` to be considered in the sum
    integer(i32),intent(in)        :: first_kpt
    !> index of the last `k-point` considered
    integer(i32),intent(in)        :: last_kpt
    !> Backup of `aind`
    real(dp), intent(in)           :: aindsave(3)
    !> Backup of `atot`
    real(dp), intent(in)           :: atotsave(3)
    !> Backup of `aext`
    real(dp), intent(in)           :: aextsave(3)
    !> Backup of `pvec`
    real(dp), intent(in)           :: pvecsave(3)
    !> Backup of `jind`
    real(dp), intent(in)           :: jindsave(3)
    !> Backup of `jpara`
    real(dp), intent(in)           :: jparasave(3)
    !> When `.True.`, it informs that the maximum steps have been reached
    logical, intent(out)           :: maxStepsReached

    integer(i32) :: i
    real(dp)     :: err

    do i = 1, maxSteps
      ! WAVEFUNCTION
      evecfv_time = evecfv_save
      call UpdateWavefunction( .True. )

      ! Update the paramagnetic component of the induced current density
      call UpdateCurrentDensity( first_kpt, last_kpt, evecfv_time(:,:,:),jparanext(:))
      if ( input%xs%realTimeTDDFT%subtractJ0 ) jparanext(:) = jparanext(:)-jparaspurious(:)

      ! DENSITY
      call uprho( it, .False. )

      ! KS-POTENTIAL
      call uppot(.False.)

      ! We only should evolve A if fieldType == external
      if( evolveA ) then
        jpara(:) = jparasave(:) !attention: jparaold saves the value of jpara(t-deltat)
        aind(:) = aindsave(:)
        atot(:) = atotsave(:)
        aext(:) = aextsave(:)
        pvec(:) = pvecsave(:)
        jind(:) = jindsave(:)
        call update_vector_potential( time, atot, aind, aext ) 
        call Evolve_A_ind( .False. )
      end if

      ! INDUCED CURRENT
      ! Update the paramagnetic component of the induced current density
      jdia(:) = -atot(:)*chgval/c/omega
      ! Update the paramagnetic component of the induced current density
      jpara(:) = jparanext(:)
      jind(:) = jpara(:)+jdia(:)

      ! HAMILTONIAN
      ham_predcorr(:,:,:) = ham_time(:,:,:)
      call UpdateHam( predcorr=.True., calculateOverlap=.False. )

      ! Check the difference between the two hamiltonians
      err = maxval(abs(ham_predcorr(:,:,:)-ham_time(:,:,:)))
      if ( err <= tolPredCorr ) exit

    end do
    maxStepsReached = (i>maxSteps)
  end subroutine 

  subroutine deallocate_global_arrays

    deallocate( apwalm, evecfv_gnd, evecsv, evecfv_time )
    deallocate( overlap, ham_time, ham_past, pmat )
    if ( predictorCorrector ) deallocate( ham_predcorr, evecfv_save )
    if ( associated(input%xs%realTimeTDDFT%laser) ) then
      if ( nkicks >= 1 ) then
        deallocate( wkick, dirkick, amplkick, t0kick )
      end if
      if ( ntrapcos >= 1 ) then
        deallocate( dirtrapcos, ampltrapcos, omegatrapcos, phasetrapcos )
        deallocate( t0trapcos, trtrapcos, wtrapcos )
      end if
      if ( nsinsq >= 1 ) then
        deallocate( dirsinsq, amplsinsq, omegasinsq )
        deallocate( phasesinsq, t0sinsq, tpulsesinsq )
      end if
    end if
    
    if ( allocated( timingstore ) ) deallocate( timingstore )
  end subroutine

end module rttddft_main
