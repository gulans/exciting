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
  use rttddft_GlobalVariables
  use precision, only: dp
  use rttddft_Energy, only: TotalEnergy, obtain_energy_rttddft
  use rttddft_io, only: open_file_info, write_file_info, write_file_info_header, &
    close_file_info, write_wavefunction
  use mod_kpoint, only: nkpt
#ifdef MPI
  use mpi, only: MPI_BARRIER
#endif
  use modmpi, only: rank, procs, ierr, mpi_env_k, distribute_loop

  implicit none

  private
  integer                 :: iprint
  real(dp)                :: timei, timef

  type(TimingRTTDDFT),&
              allocatable :: timingstore(:)

  private :: print_total_energy, print_nexc, print_jpa, uprho, uppot, printTiming
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
    use rttddft_HamiltonianOverlap, only: UpdateHam
    use rttddft_CurrentDensity, only: UpdateCurrentDensity
    use rttddft_Wavefunction, only: UpdateWavefunction
    use m_getunit, only: getunit
    use rttddft_NumberExcitations, only: Obtain_number_excitations
    use rttddft_screenshot, only: screenshot
    use rttddft_VectorPotential, only: Calculate_Vector_Potential, Evolve_A_ind => Solve_ODE_Vector_Potential
    use mod_lattice, only: omega
    use mod_charge_and_moment, only: chgval
    use mod_misc, only: filext
    use modinput, only: input
    use physical_constants, only: c
    use errors_warnings, only: terminate_if_false
    use modmpi, only: mpiglobal

    ! counter for the number of iterations of real-time
    integer                 :: it
    ! prints data every nprint steps of the counter "it"
    integer                 :: nprint
    ! counter for the predictor-corrector loop
    integer                 :: ipredcorr
    ! indexes of the first and the last k-points
    integer                 :: first_kpt, last_kpt

    character(50)           :: string

    real(dp)                :: timesave, timeiter
    real(dp)                :: timehml, timerest
    real(dp)                :: err
    real(dp),allocatable    :: nex(:), ngs(:), nt(:)
    real(dp)                :: aindsave(3),aindbackup(3),pvecsave(3)
    real(dp)                :: jindsave(3),aextsave(3),atotsave(3)

    ! Variables to store data and print
    real(dp),allocatable    :: timestore(:),aindstore(:,:),atotstore(:,:)
    real(dp),allocatable    :: jindstore(:,:), pvecstore(:,:)

    type(TotalEnergy),allocatable &
                            :: etotstore(:)


    call timesec(timesave)

    ! Variable that tells after how many iterations the output should be printed
    nprint = input%xs%realTimeTDDFT%printAfterIterations
    ! Consistency check: check if nprint > 0. If not, stop the execution.
    call terminate_if_false( mpiglobal, nprint > 0, &
      & 'Invalid non-positive number as printAfterIterations' )

    ! Consistency check: check if matrices are not stored in a packed form.
    ! If they are, stop the execution.
    call terminate_if_false( mpiglobal, .not. input%groundstate%solver%packedmatrixstorage, &
      & 'Error: RT-TDDFT does not work with matrices stored in a packed form.' )

    ! Consistency check: check if no spin polarized calculations are requested.
    ! In the spin-polarized case, stop the execution
    call terminate_if_false( mpiglobal, .not. input%groundstate%tevecsv, &
      & 'Error: only spin unpolarised calculations are possible with RT-TDDFT now.' )

    ! Outputs general info to RTTDDFT_INFO.OUT and
    ! opens TIMING_RTTDDFT.OUT (if this is the case)
    if( rank == 0 ) then
      if(input%xs%realTimeTDDFT%printTimingGeneral) then
        call timesec (timei)
        call getunit(filetime)
        open(filetime,file='TIMING_RTTDDFT'//trim(filext),status='replace')
      end if
      call open_file_info
      call write_file_info_header
    end if

    ! Allocate variables to be stored and printed only after nprint steps
    allocate(timestore(nprint))
    allocate(aindstore(3,nprint))
    allocate(atotstore(3,nprint))
    allocate(jindstore(3,nprint))
    allocate(pvecstore(3,nprint))
    if( input%xs%realTimeTDDFT%printTimingGeneral )  allocate(timingstore(nprint))
    if( input%xs%realTimeTDDFT%calculateTotalEnergy ) allocate(etotstore(nprint))
    if( input%xs%realTimeTDDFT%calculateNExcitedElectrons ) allocate(nex(nprint),ngs(nprint),nt(nprint))

    call initialize_rttddft


    ! Print the vector potential, the current density and the polarization vector
    ! at t=0
    timestore(1) = time ! Trick: we need an array to call the subroutine print_jpa
    ! vector potential
    call print_jpa( fileavec, 1, timestore(1), aind(:), &
      & atot(:)  )
    ! current density
    call print_jpa( filejind, 1, timestore(1), jind(:) )
    ! polarization vector
    call print_jpa( filepvec, 1, timestore(1), pvec(:) )

    ! Initialize integers that contain the first and last k-point
    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    ! Total energy, number of excitations, and screenshots at t=0
    ! Total energy
    if ( calculateTotalEnergy ) then
      call potcoul
      call potxc
      call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_gnd, etotstore(1) )
      ! Trick: we need an array to call the subroutine print_total_energy
      timestore(1) = time
      call print_total_energy( fileetot, .True., 1, timestore(1), etotstore(1) )
    end if

    ! Number of excitations
    if (calculateNexc) then
      call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
        & evecfv_time, overlap, nex(1), ngs(1), nt(1) )
      ! Trick: we need an array to call the subroutine print_nexc
      timestore(1) = time
      call print_nexc( filenexc, .True., 1, timestore(1), nex(1), ngs(1), nt(1) )
    end if

    ! Screenshot
    if ( associated(input%xs%realTimeTDDFT%screenshots) ) then
      call screenshot( 0, first_kpt, last_kpt, overlap, evecfv_gnd, &
        & evecfv_time, ham_time )
    end if


    if( printTimesGeneral ) then
      call timesec( timef )
      if( rank == 0 ) write( filetime, formatTime ) 'Initialization (sec):',timef-timesave
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
      ! Update the induced part of the vector potential
      if( input%xs%realTimeTDDFT%laser%fieldType == 'external' ) then
        ! Check if we need to save aind, pvec, atot and aext
        if( predictorCorrector .and. (input%xs%realTimeTDDFT%vectorPotentialSolver /= 'euler') ) then
          aindsave(:) = aind(:)
          pvecsave(:) = pvec(:)
          atotsave(:) = atot(:)
          aextsave(:) = aext(:)
        end if
        call Evolve_A_ind( .False. )
        ! Update aext
        call Calculate_Vector_Potential( time, aext(:) )
        ! Update the (total) vector potential
        atot(:) = aind(:) + aext(:)
      else
        call Evolve_A_ind( .True. )
        call Calculate_Vector_Potential( time-tstep, aext(:) ) ! This is a trick: aext acts as an auxiliary variable
        call Calculate_Vector_Potential( time+tstep, atot(:) )
        call Calculate_Vector_Potential( time, atot(:) )
      end if

      ! INDUCED CURRENT
      ! Update the diamagnetic component of the induced current density
      jdia(:) = - atot(:) * chgval / c / omega
      ! Update the paramagnetic component of the induced current density
      jparaold(:) = jpara(:)
      jpara(:) = jparanext(:)
      ! Update the total induced current
      if ( predictorCorrector .and. (input%xs%realTimeTDDFT%vectorPotentialSolver /= 'euler') ) then
        jindsave(:) = jind(:)
      end if
      jind(:) = jpara(:) + jdia(:)

      if(printTimesGeneral) call timesecRTTDDFT(timei,timef,timingstore(iprint)%t_obtaina)

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

      if ( predictorCorrector .and. (method /= 'SE') .and. (method /= 'EH') ) then
        ! Remark: it makes no sense to employ the predictor-corrector method with SE or EH!
        do ipredcorr = 1, maxStepsPredictorCorrector
          ! WAVEFUNCTION
          evecfv_time(:,:,:) = evecfv_save(:,:,:)
          call UpdateWavefunction( .True. )

          ! Update the paramagnetic component of the induced current density
          call UpdateCurrentDensity( first_kpt, last_kpt, evecfv_time(:,:,:),jparanext(:))
          if ( input%xs%realTimeTDDFT%subtractJ0 ) jparanext(:) = jparanext(:)-jparaspurious(:)

          ! DENSITY
          call uprho( it, .False. )

          ! KS-POTENTIAL
          call uppot(.False.)

          ! VECTOR POTENTIAL
          ! Update the induced part of the vector potential
          if( (input%xs%realTimeTDDFT%laser%fieldType .eq. 'external') .and. (input%xs%realTimeTDDFT%vectorPotentialSolver .ne. 'euler') ) then
            jpara(:) = jparaold(:) !attention: jparaold saves the value of jpara(t-deltat)
            aindbackup(:) = aind(:)
            aind(:) = aindsave(:)
            atot(:) = atotsave(:)
            aext(:) = aextsave(:)
            pvec(:) = pvecsave(:)
            jind(:) = jindsave(:)
            call Evolve_A_ind( .False. )
            call Calculate_Vector_Potential( time, aext(:) )
            ! Update the (total) vector potential
            atot(:) = aind(:) + aext(:)
          end if

          ! INDUCED CURRENT
          ! Update the paramagnetic component of the induced current density
          jdia(:) = -atot(:)*chgval/c/omega
          ! Update the paramagnetic component of the induced current density
          jpara(:) = jparanext(:)
          jind(:) = jpara(:)+jdia(:)

          ! HAMILTONIAN
          ham_predcorr(:,:,:) = ham_time(:,:,:)
          call UpdateHam(.True.,.False.)

          ! Check the difference between the two hamiltonians
          err = maxval(abs(ham_predcorr(:,:,:)-ham_time(:,:,:)))
          if ( err .le. tolPredCorr ) exit

        end do ! do ipredcorr = 1, maxstepsPredictorCorrector
        if ( ipredcorr >= maxstepsPredictorCorrector ) then
          if( rank == 0 ) write(*,*) 'Problems with convergence (PredCorr), time: ', time
        end if
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

      ! Check if a screenshot has been requested
      if ( associated(input%xs%realTimeTDDFT%screenshots) ) then
        if ( mod( it, input%xs%realTimeTDDFT%screenshots%niter ) == 0 ) then
          call screenshot( it, first_kpt, last_kpt, overlap, evecfv_gnd, &
            & evecfv_time, ham_time )
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
        ! Print the vector potential
        call print_jpa( fileavec, nprint, timestore(:), aindstore(:,:), &
          & atotstore(:,:)  )
        ! Print the current density
        call print_jpa( filejind, nprint, timestore(:), jindstore(:,:) )
        ! Print the polarization vector
        call print_jpa( filepvec, nprint, timestore(:), pvecstore(:,:) )
        ! Print Total Energy - if requested
        if ( calculateTotalEnergy ) then
          call print_total_energy( fileetot, .False., nprint, timestore(:), &
            & etotstore(:) )
        end if
        ! Print number of excitations - if requested
        if ( calculateNexc ) then
          call print_nexc( filenexc, .False., nprint, timestore(:), nex(:), &
            & ngs(:), nt(:) )
        end if
        ! Update the counter
        iprint = 1
        if( printTimesGeneral ) then
          call timesecRTTDDFT( timeiter, timef, timingstore(nprint)%t_iteration )
          call printtiming( it, nprint, timingstore )
        end if
      else ! if ( iprint .eq. nprint ) then
        if(printTimesGeneral) call timesecRTTDDFT(timeiter,timef,timingstore(iprint)%t_iteration)
        iprint = iprint + 1
      end if ! if ( iprint .eq. nprint ) then
#ifdef MPI
      ! Make all the processes wait here
      ! The master alone has been writing the files above
      call MPI_BARRIER(mpi_env_k%comm,ierr)
#endif
    end do ! do it = 1, nsteps

    if ( rank == 0 ) then
      close(fileavec)
      close(filejind)
      close(filepvec)
      call write_file_info( 'Real-time TDDFT calculation finished' )
      call close_file_info
      if( calculateTotalEnergy ) close(fileetot)
      if( calculateNexc ) close(filenexc)
      if( printTimesGeneral ) close(fileTime)
    end if

    ! write wavefunction, and potential and density with _RTTDDFT.OUT as suffix
    string = filext
    filext = '_RTTDDFT'//trim(filext)
    call write_wavefunction( first_kpt, evecfv_time )
    if ( rank == 0 ) call writestate
    filext = string

    !----
    ! deallocate session
    deallocate( apwalm, evecfv_gnd, evecsv, evecfv_time )
    deallocate( overlap, ham_time, ham_past, pmat )
    deallocate( timestore, aindstore, atotstore, jindstore, pvecstore )
    if ( predictorCorrector ) then
      deallocate( ham_predcorr, evecfv_save )
    end if
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
    if( calculateTotalEnergy ) deallocate( etotstore )
    if( calculateNexc ) deallocate( nex, ngs, nt )

  end subroutine coordinate_rttddft_calculation

  !> Prints the total energy
  subroutine print_total_energy( fileUnitNumber, printHeader, nArrayElements, &
    & timeArray, etotArray )

    implicit none

    !> The unit-number of the output file
    integer, intent(in) :: fileUnitNumber
    !> Number of lines to printed = number of elements of the arrays:
    !> `timeArray` and `etotArray`.
    integer, intent(in) :: nArrayElements
    !> Tells if a header must be printed (useful when the file is opened for 
    !> the 1st time)
    logical, intent(in) :: printHeader
    !> Array with the values of time \( t \)
    real(8), intent(in) :: timeArray(nArrayElements)
    !> Array with the energies (total energy, XC, Madelung, etc.)
    type(TotalEnergy), intent(in) :: etotArray(nArrayElements)

    integer :: i

    if ( rank == 0 ) then
      if ( printHeader ) then
        write(fileUnitNumber,'(A9,8A20)') 'Time','ETOT','Madelung','Eigenvalues-Core',&
          & 'Eigenvalues-Valence','Exchange','Correlation','XC-potential',&
          & 'Coulomb pot. energy'
      end if
      do i = 1, nArrayElements
        write(fileUnitNumber,'(F9.3,8F20.10)') timeArray(i), &
          & etotArray(i)%total_energy, etotArray(i)%madelung, &
          & etotArray(i)%eigenvalues_core, etotArray(i)%hamiltonian,&
          & etotArray(i)%exchange, etotArray(i)%correlation, &
          & etotArray(i)%integral_vxc_times_density, etotArray(i)%Coulomb
      end do
    end if
  end subroutine print_total_energy

  !> Prints the number of excitations
  subroutine print_nexc( fileUnitNumber, printHeader, nArrayElements, &
    & timeArray, nex, ngs, ntot )

    implicit none

    !> The unit-number of the output file
    integer, intent(in) :: fileUnitNumber
    !> Number of lines to printed = number of elements of the arrays:
    !> `nex`, `ngs` and `ntot`.
    integer, intent(in)   :: nArrayElements
    !> printHeader: If we need to print a header (useful when we open the file for the 1st time)
    logical, intent(in)   :: printHeader
    !> timeArray         array with the values of time \( t \)
    real(dp), intent(in)  :: timeArray(nArrayElements)
    !> nex               Number of electrons which were excited
    real(dp), intent(in)  :: nex(nArrayElements)
    !> ngs               Number of electrons on the groundstate
    real(dp), intent(in)  :: ngs(nArrayElements)
    !> ntot              Sum of ngs and nex
    real(dp), intent(in)  :: ntot(nArrayElements)

    integer :: i

    if ( rank == 0 ) then
      if ( printHeader ) then
        write( fileUnitNumber,'(A9,3A20)' ) 'Time','N.Elec.GS', &
          & 'N.XS', 'Sum'
      end if
      do i = 1, nArrayElements
        write( fileUnitNumber, '(F9.3,3F20.10)' ) timeArray(i), &
          & ngs(i), nex(i), ntot(i)
      end do
    end if
  end subroutine print_nexc


  !> Prints the current density \(\mathbf{J}\), or the polarization 
  !> \(\mathbf{P}\), or the vector potential \(\mathbf{A}\)
  subroutine print_jpa( fileUnitNumber, nArrayElements, timeArray, &
    & firstArray, secondArray )

    implicit none
   
    !> The unit-number of the output file
    integer, intent(in) :: fileUnitNumber
    !> Number of lines to printed = number of elements of the arrays:
    !> `timeArray`, `firstArray`, and eventually the second dimension of `secondArray`.
    integer, intent(in)   :: nArrayElements
    !> Array with the values of time \( t \)
    real(8), intent(in) :: timeArray(nArrayElements)
    !> Array with the \( x, y, z \) components of \(\mathbf{J}\) , 
    !> \(\mathbf{P}\) or \(\mathbf{A}\) for each time \( t \)
    real(8), intent(in) :: firstArray(3, nArrayElements)
    !> Same as before, but for the second array - usually \(\mathbf{A}\)
    real(8), intent(in), optional :: secondArray(3, nArrayElements)

    integer :: i
    logical :: twoArraysPresent

    twoArraysPresent = present( secondArray )

    if ( rank == 0 ) then
      do i = 1, nArrayElements
        if( twoArraysPresent ) then
          write( fileUnitNumber, '(F9.3,6F20.12)' ) timeArray(i), &
            & firstArray(1, i), secondArray(1, i), firstArray(2, i), &
            & secondArray(2, i), firstArray(3, i), secondArray(3, i)
        else
          write( fileUnitNumber, '(F9.3,3F20.12)' ) timeArray(i), &
            & firstArray(1, i), firstArray(2, i), firstArray(3, i)
        end if
      end do
    end if
  end subroutine print_jpa


  !> This is just an interface to call the subroutine `[[UpdateDensity]]`, which
  !> updates the charge density
  subroutine uprho( iteration_counter, computeTiming )
    use rttddft_Density, only: UpdateDensity

    implicit none

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

    implicit none
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


  !> Subroutine to output the timings into `TIMING_RTTDDFT.OUT`
  subroutine printTiming( itNumber, n, timing )
    use modinput, only: input

    implicit none

    !> itNumber: The actual number of the counter that tells how many time steps 
    !> have already been executed
    integer, intent(in)             :: itNumber
    !> n:  Number of elements of the timing array
    integer, intent(in)             :: n
    !> timing: Array of timings. Each elements contains information
    !>   about how many seconds (timings) were spent in different parts of code
    type(TimingRTTDDFT), intent(in) :: timing(n)

    integer                         :: ip, shift

    shift = itNumber - n

    if ( rank == 0 ) then
      do ip = 1, n
        write(filetime,'(A30,I10)')'Time (sec) spent in iteration:',ip+shift
        write(filetime,formatTime) 'updatewvf:',timing(ip)%t_wvf
        write(filetime,formatTime) 'updatedens:',timing(ip)%t_dens
        if ( printTimesDetailed ) then
          write(filetime,formatTime) '-- rhovalk and genrhoir:',timing(ip)%t_dens_rho
          write(filetime,formatTime) '-- symrf:',timing(ip)%t_dens_symrf
          write(filetime,formatTime) '-- rfmtctof:',timing(ip)%t_dens_rfmtctof
          write(filetime,formatTime) '-- addrhocr:',timing(ip)%t_dens_addrhocr
          write(filetime,formatTime) '-- charge:',timing(ip)%t_dens_charge
          write(filetime,formatTime) '-- rhonorm:',timing(ip)%t_dens_rhonorm
        end if
        write(filetime,formatTime) 'updatepot:',timing(ip)%t_uppot
        if ( printTimesDetailed ) then
          write(filetime,formatTime) '-- poteff:',timing(ip)%t_poteff
          write(filetime,formatTime) '-- genveffig:',timing(ip)%t_genveffig
          write(filetime,formatTime) '-- genmeffig:',timing(ip)%t_genmeffig
        end if
        write(filetime,formatTime) 'UpdateCurrentDensity:',timing(ip)%t_curr
        write(filetime,formatTime) 'ObtainA:',timing(ip)%t_obtaina
        write(filetime,formatTime) 'updatehamiltonian:',timing(ip)%t_upham
        if ( printTimesDetailed ) then
          write(filetime,formatTime) '-- hmlint:',timing(ip)%t_hmlint
          write(filetime,formatTime) '-- other subs:',timing(ip)%t_ham
        end if

        ! Predictor-corrector
        if ( predictorCorrector )  &
          & write(filetime,formatTime) &
          & 'All cycles of predcorr:',timing(ip)%t_predcorr

        ! Total Energy
        if ( calculateTotalEnergy .and. printTimesDetailed ) write(filetime,formatTime) 'Total Energy:',timing(ip)%t_toten
        if ( calculateNexc .and. printTimesDetailed ) write(filetime,formatTime)'nexc:',timing(ip)%t_nexc
        if ( associated(input%xs%realTimeTDDFT%screenshots) ) write(filetime,formatTime) 'Screenshots:',timing(ip)%t_screenshot
        write(filetime,formatTime) 'time per iteration:',timing(ip)%t_iteration
      end do
    end if
  end subroutine printtiming

end module rttddft_main
