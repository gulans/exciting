!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.



!> Performs basic consistency checks as well as allocating and initialising
!> global variables not dependent on the $k$-point set.
!>
!> NOTE: No one should extend this routine with a bare implementation.
!> New code should be in initialisation routines, and called here.
!>
!> !REVISION HISTORY:
!>   Created January 2004 (JKD)
!>   Started clean-up. 2022 (ABuccheri)
Subroutine init0
      use modinteg
      use modbess
      Use modinput
      Use modmain
      Use autormt, only: optimal_rmt
      use cmd_line_args, only: cmd_line_args_type
      Use modxcifc
      Use modmpi, only: mpiglobal, mpi_env_k, mpi_env_band, find_2d_grid
#ifdef MPI
      use modmpi, only: MPI_COMM_SELF
#endif
      Use errors_warnings, only: terminate_if_false
      Use vx_enums, only: HYB_PBE0, HYB_HSE
      Use APW_basis_size, only: determine_rgkmax, determine_APWprecision
      ! TODO(ALEX) Once everyone has done their refactor
      ! initialisation routines should be moved from tmp_mod_init0
      ! and tmp_mod_init0 should be deleted.
      Use tmp_mod_init0
#ifdef XS
      Use modxs
#endif
      use sirius_init, only: sirius_options
      use sirius_api, only: setup_sirius, get_mpi_comm_sirius, gengvec_sirius, warn_array_sizes_sirius,&
                            set_periodic_function_ptr_sirius

      Implicit None

      Integer :: is, js, ia, ias
      Integer :: ist, l, m, lm, iv (3)
      Real (8) :: ts0, ts1, tv3 (3)

      integer :: rows_per_kpt, cols_per_kpt
      integer :: mpi_grid(2)
      integer :: ierr

      !> Bare mpi communicator for k-points
      integer :: comm_k
      !> Bare mpi communicator for bands.
      !> This will only be set when using sirius
      integer :: comm_band
      real(8) :: mb, ylmg_mb, sfacg_mb
      !> Command line arguments
      type(cmd_line_args_type) :: args

      call stopwatch("exciting:init0", 1)

      ! Zero self-consistent loop number
      iscl = 0
      tlast = .False.

      call initialise_groundstate_timings()
      call timesec(ts0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
      ntpll = 770
      lmmaxhf = (2*input%groundstate%lmaxvr+1) ** 2
      lmmaxvr = (input%groundstate%lmaxvr+1) ** 2
      lmmaxapw = (input%groundstate%lmaxapw+1) ** 2
      lmmaxmat = (input%groundstate%lmaxmat+1) ** 2
      lmmaxinr = (input%groundstate%lmaxinr+1) ** 2
      If (input%groundstate%lmaxvr .Gt. input%groundstate%lmaxapw) Then
         Write (*,*)
         Write (*, '("Error(init0): lmaxvr > lmaxapw : ", 2I8)') &
        & input%groundstate%lmaxvr, input%groundstate%lmaxapw
         Write (*,*)
         Stop
      End If
      If (input%groundstate%lmaxmat .Gt. input%groundstate%lmaxapw) &
     & Then
         Write (*,*)
         Write (*, '("Error(init0): lmaxmat > lmaxapw : ", 2I8)') &
        & input%groundstate%lmaxmat, input%groundstate%lmaxapw
         Write (*,*)
         Stop
      End If
! index to (l,m) pairs
      If (allocated(idxlm)) deallocate (idxlm)
      Allocate (idxlm(0:input%groundstate%lmaxapw,-&
     & input%groundstate%lmaxapw:input%groundstate%lmaxapw))
      lm = 0
      Do l = 0, input%groundstate%lmaxapw
         Do m = - l, l
            lm = lm + 1
            idxlm (l, m) = lm
         End Do
      End Do
! array of i**l values
      If (allocated(zil)) deallocate (zil)
      Allocate (zil(0:input%groundstate%lmaxapw))
      Do l = 0, input%groundstate%lmaxapw
         zil (l) = zi ** l
      End Do
!
!------------------------------------!
!     index to atoms and species     !
!------------------------------------!

! find primitive cell if required
      If (input%structure%primcell) Call findprim
!
      call map_atoms_per_species_to_atomic_index (nspecies, natoms, idxas, natmmax, natmtot)

!------------------------!
!     spin variables     !
!------------------------!
      If (isspinspiral()) Then
         Select Case (task)
         Case (2, 3, 7, 15, 51, 52, 53, 61, 62, 63, 120, 121)
            Write (*,*)
            Write (*, '("Error(init0): spin-spirals do not work with ta&
           &sk ", I4)') task
            Write (*,*)
            Stop
         End Select
          If (associated(input%groundstate%OEP)) Then
            Write (*,*)
            Write (*, '("Error(init0): spin-spirals do not work with th&
           &e OEP method")')
            Write (*,*)
            Stop
         End If
      End If
! spin-orbit coupling or fixed spin moment implies spin-polarised calculation
!
! number of spinor components and maximum allowed occupancy
      If (associated(input%groundstate%spin)) Then
         nspinor = 2
         occmax = 1.d0
      Else
         nspinor = 1
         occmax = 2.d0
      End If
! number of spin-dependent first-variational functions per state
      If (isspinspiral()) Then
         nspnfv = 2
      Else
         nspnfv = 1
      End If
! spin-polarised calculations require second-variational eigenvectors
      If (associated(input%groundstate%spin)) input%groundstate%tevecsv &
     & = .True.
! Hartree-Fock/RDMFT requires second-variational eigenvectors
      If ((task .Eq. 5) .Or. (task .Eq. 6) .Or. (task .Eq. 300)) &
     & input%groundstate%tevecsv = .True.

     call initialise_xc_mixing_coefficients(input%groundstate, xctype, xcdescr, xcspin, xcgrad, ex_coef, ec_coef)
      
! reset input%groundstate%Hybrid%excoeff to ex_coef
! in case of libxc: overwritten by ex_coef as defined by libxc
      If (associated(input%groundstate%Hybrid)) input%groundstate%Hybrid%excoeff = ex_coef

      If ((associated(input%groundstate%spin)) .And. (xcspin .Eq. 0)) Then
         Write (*,*)
         Write (*, '("Error(init0): requested spin-polarised run with spin-unpolarised")')
         Write (*, '(" exchange-correlation functional")')
         Write (*,*)
         Stop
      End If

! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
      If (associated(input%groundstate%spin)) Then
         ndmag = 1
         If ((Abs(input%groundstate%spin%bfieldc(1)) .Gt. &
        & input%structure%epslat) .Or. &
        & (Abs(input%groundstate%spin%bfieldc(2)) .Gt. &
        & input%structure%epslat)) ndmag = 3
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               If ((Abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(1)) .Gt. input%structure%epslat) .Or. &
              & (Abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(2)) .Gt. input%structure%epslat)) ndmag = 3
            End Do
         End Do
! source-free fields and spin-spirals are non-collinear in general
         If ((input%groundstate%nosource) .Or. (isspinspiral())) ndmag &
        & = 3
! spin-orbit coupling is non-collinear in general
         If (isspinorb()) ndmag = 3
      Else
         ndmag = 0
      End If

! set the non-collinear flag
      If (ndmag .Eq. 3) Then
         ncmag = .True.
      Else
         ncmag = .False.
      End If
      If ((ncmag) .And. (xcgrad .Gt. 0)) Then
         call warning('Warning(init0):')
         call warning(' GGA inconsistent with non-collinear magnetism')
      End If
! set fixed spin moment effective field to zero
      bfsmc (:) = 0.d0
! set muffin-tin FSM fields to zero
      bfsmcmt (:, :, :) = 0.d0
!
!-------------------------------------!
!     lattice and symmetry set up     !
!-------------------------------------!
! for convenience - create a copy
      avec(:,1) = input%structure%crystal%basevect(:,1)
      avec(:,2) = input%structure%crystal%basevect(:,2)
      avec(:,3) = input%structure%crystal%basevect(:,3)

! generate the reciprocal lattice vectors and unit cell volume
      Call reciplat
! compute the inverse of the lattice vector matrix
      Call r3minv (input%structure%crystal%basevect, ainv)
! compute the inverse of the reciprocal vector matrix
      Call r3minv (bvec, binv)
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
! map atomic lattice coordinates to [0,1) if not in molecule mode
             If ( .Not. input%structure%cartesian) Call r3frac (input%structure%epslat, &
           & input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), iv)
! determine atomic Cartesian coordinates
            Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
           & atposc(:, ia, is))
! lattice coordinates of the muffin-tin magnetic fields
            Call r3mv (ainv, input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:), bflmt(:, ia, is))
         End Do
      End Do
! lattice coordinates of the global magnetic field
      If (associated(input%groundstate%spin)) Then
         tv3 = input%groundstate%spin%bfieldc
      Else
         tv3 = 0
      End If
      Call r3mv (ainv, tv3, bfieldl)
! Cartesian coordinates of the spin-spiral vector
      If (associated(input%groundstate%spin)) Then
         tv3 = input%groundstate%spin%vqlss
      Else
         tv3 = 0
      End If
      Call r3mv (bvec, tv3, vqcss)
! find Bravais lattice symmetries
      Call findsymlat
! use only the identity if required
      If (input%groundstate%nosym) nsymlat = 1
! find the crystal symmetries and shift atomic positions if required
      Call findsymcrys
! find the site symmetries
      Call findsymsite
#ifdef XS
! determine inverse symmery elements
      Call findsymi (input%structure%epslat, maxsymcrys, nsymcrys, &
     & symlat, lsplsymc, vtlsymc, isymlat, scimap)
! generate symmetrization array for rank 2 tensors
      Call gensymt2 (maxsymcrys, nsymcrys, symlatc, lsplsymc, symt2)
! calculate advanced information on symmetry group
      Call setupsym
#endif

! automatically determine the muffin-tin radii if required
      If (input%structure%autormt .and. (idx_species_fixed_rmt .gt. 0)) then 
            Call optimal_rmt(rmt, spzn, input%structure%crystal%basevect, atposc,&
                              &input%structure%autormtscaling, natoms, nspecies, 1, idx_species_fixed_rmt)         
      else if (input%structure%autormt) then 
            Call optimal_rmt(rmt, spzn, input%structure%crystal%basevect, atposc,&
                              &input%structure%autormtscaling, natoms, nspecies, 1)         
      end if 

! check for overlapping muffin-tins
      Call checkmt
!
!-----------------------!
!     radial meshes     !
!-----------------------!
      nrmtmax = 1
      nrcmtmax = 1
      js = 1
      Do is = 1, nspecies
! make the muffin-tin mesh commensurate with lradstp
         nrmt (is) = nrmt (is) - Mod (nrmt(is)-1, &
        & input%groundstate%lradstep)
         nrmtmax = Max (nrmtmax, nrmt(is))
! number of coarse radial mesh points
         nrcmt (is) = (nrmt(is)-1) / input%groundstate%lradstep + 1
         nrcmtmax = Max (nrcmtmax, nrcmt(is))
! smallest muffin-tin radius
         If (rmt(is) .Lt. rmt(js)) js = is
      End Do
      If ((input%groundstate%isgkmax .Lt. 1) .Or. &
     & (input%groundstate%isgkmax .Gt. nspecies)) &
     & input%groundstate%isgkmax = js
! set up atomic and muffin-tin radial meshes
      Call genrmesh

      If (input%groundstate%useAPWprecision) Then
            input%groundstate%rgkmax = determine_rgkmax(input%groundstate%APWprecision, spzn(input%groundstate%isgkmax))
      Else
            input%groundstate%APWprecision = determine_APWprecision(input%groundstate%rgkmax, spzn(input%groundstate%isgkmax))
      End If


!-----------------------!
! initialize modinteg   !
!-----------------------!


if (allocated(mt_integw%fintw)) then
  call dealloc_icoef()
  call gen_icoef(nspecies,spnrmax,nrmt,spnr,spr)

  else
  call gen_icoef(nspecies,spnrmax,nrmt,spnr,spr)
  endif

!-----------------------!
! initialize modbess    !
!-----------------------!

 call init_bess(nrmtmax,nspecies,nrmt,spr(1:nrmtmax,:))


!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
      chgzn = 0.d0
      chgcr = 0.d0
      chgval = 0.d0
      spnstmax = 0
      Do is = 1, nspecies
! nuclear charge
         chgzn = chgzn + spzn (is) * dble (natoms(is))
! find the maximum number of atomic states
         spnstmax = Max (spnstmax, spnst(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
         spze (is) = 0.d0
         Do ist = 1, spnst (is)
            spze (is) = spze (is) + spocc (ist, is)
            If (spcore(ist, is)) Then
               chgcr = chgcr + dble (natoms(is)) * spocc (ist, is)
            Else
               chgval = chgval + dble (natoms(is)) * spocc (ist, is)
            End If
         End Do
      End Do
! add excess charge
      chgval = chgval + input%groundstate%chgexs
! total charge
      chgtot = chgcr + chgval
      If (chgtot .Lt. 1.d-8) Then
         Write (*,*)
         Write (*, '("Error(init0): zero total charge")')
         Write (*,*)
         Stop
      End If
! effective Wigner radius
      rwigner = (3.d0/(fourpi*(chgtot/omega))) ** (1.d0/3.d0)

! find the G-vector grid sizes
      Call gridsize

      ! Setup SIRIUS and MPI communicators
      if ( associated(input%groundstate%sirius) ) then
        call sirius_options%initialise(input)
        call args%parse(mpiglobal)
        call find_2d_grid(mpiglobal%procs, args%kptgroups, rows_per_kpt, &
                          cols_per_kpt, "k-points")
        mpi_grid = [rows_per_kpt, cols_per_kpt]
        call setup_sirius(input, mpiglobal%comm, mpi_grid)
        call get_mpi_comm_sirius(comm_k, comm_band)
      else
        ! Leave simple layout for pure Exciting, with no band distribution
        comm_k = mpiglobal%comm
        comm_band = 0
#ifdef MPI
        comm_band = MPI_COMM_SELF
#endif
      end if

      call mpi_env_band%init(comm_band)
      call mpi_env_k%init(comm_k)

!-------------------------!
!     G-vector arrays     !
!-------------------------!
! generate the G-vectors
      if (associated(input%groundstate%sirius)) then
        call gengvec_sirius(ivg, ivgig, igfft, vgc, gc, ngvec, intgv)
      else
        call gengvec
      end if

      ! Check for ylmg and sfacg, which will become an issue
      ! if running sirius with many MPI instances per node,
      ! for large systems
      call warn_array_sizes_sirius(lmmaxvr, ngvec)

      call genylmg

! allocate structure factor array for G-vectors
      If (allocated(sfacg)) deallocate (sfacg)
      Allocate (sfacg(ngvec, natmtot))
! generate structure factors for G-vectors
      Call gensfacgp (ngvec, vgc, ngvec, sfacg)
! generate the characteristic function
      Call gencfun
!
!-------------------------!
!     atoms and cores     !
!-------------------------!
#ifdef XS
      If (init0symonly) Go To 10
#endif
! solve the Kohn-Sham-Dirac equations for all atoms
      Call allatoms(1)
! allocate core state eigenvalue array and set to default
      If (allocated(evalcr)) deallocate (evalcr)
      Allocate (evalcr(spnstmax, natmtot))
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
               evalcr (ist, ias) = speval (ist, is)
            End Do
         End Do
      End Do
! allocate core state radial wavefunction array
      If (allocated(rwfcr)) deallocate (rwfcr)
      Allocate (rwfcr(spnrmax, 2, spnstmax, natmtot))
! allocate core state charge density array
      If (allocated(rhocr)) deallocate (rhocr)
      Allocate (rhocr(spnrmax, natmtot))
#ifdef XS
10    Continue
#endif
!
!---------------------------------------!
!     charge density and potentials     !
!---------------------------------------!
! allocate charge density arrays
      If (allocated(rhomt)) deallocate (rhomt)
      Allocate (rhomt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(rhoir)) deallocate (rhoir)
      Allocate (rhoir(ngrtot))
! allocate magnetisation arrays
      If (allocated(magmt)) deallocate (magmt)
      If (allocated(magir)) deallocate (magir)
      If (associated(input%groundstate%spin)) Then
         Allocate (magmt(lmmaxvr, nrmtmax, natmtot, ndmag))
         Allocate (magir(ngrtot, ndmag))
      End If

call allocate_coulomb_potentials(lmmaxvr, nrmtmax, natmtot, ngrtot, vclmt, vclir, vmad)

! exchange-correlation potential
      If (allocated(vxcmt)) deallocate (vxcmt)
      Allocate (vxcmt(lmmaxvr, nrmtmax, natmtot))
      !if (allocated(vxmt)) deallocate(vxmt)
      !allocate(vxmt(lmmaxvr,nrmtmax,natmtot))
      !if (allocated(vcmt)) deallocate(vcmt)
      !allocate(vcmt(lmmaxvr,nrmtmax,natmtot))
      If (allocated(vxcir)) deallocate (vxcir)
      Allocate (vxcir(ngrtot))
! exchange-correlation magnetic field
      If (allocated(bxcmt)) deallocate (bxcmt)
      If (allocated(bxcir)) deallocate (bxcir)
      If (associated(input%groundstate%spin)) Then
         Allocate (bxcmt(lmmaxvr, nrmtmax, natmtot, ndmag))
         Allocate (bxcir(ngrtot, ndmag))
      End If
! exchange energy density
      If (allocated(exmt)) deallocate (exmt)
      Allocate (exmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(exir)) deallocate (exir)
      Allocate (exir(ngrtot))
! correlation energy density
      If (allocated(ecmt)) deallocate (ecmt)
      Allocate (ecmt(lmmaxvr, nrmtmax, natmtot))
      If (allocated(ecir)) deallocate (ecir)
      Allocate (ecir(ngrtot))
! effective potential
      If (allocated(veffmt)) deallocate (veffmt)
      Allocate (veffmt(lmmaxvr, nrmtmax, natmtot))
!      If (allocated(vrefmt)) deallocate (vrefmt)
!      Allocate (vrefmt(lmmaxvr, nrmtmax, natmtot))

      If (allocated(veffir)) deallocate (veffir)
      Allocate (veffir(ngrtot))
      If (allocated(veffig)) deallocate (veffig)
      Allocate (veffig(ngvec))
!      If (allocated(vrefig)) deallocate (vrefig)
!      Allocate (vrefig(ngvec))

      if ( associated(input%groundstate%sirius) ) then
        call set_periodic_function_ptr_sirius(rhomt, veffmt, magmt, lmmaxvr, nrmtmax, natmtot)
      end if

! allocate muffin-tin charge and moment arrays
      If (allocated(chgmt)) deallocate (chgmt)
      Allocate (chgmt(natmtot))
      If (allocated(mommt)) deallocate (mommt)
      Allocate (mommt(3, natmtot))
!
!--------------------------------------------!
!     forces and structural optimisation     !
!--------------------------------------------!
      If (allocated(forcehf)) deallocate (forcehf)
      Allocate (forcehf(3, natmtot))
      If (allocated(forcecr)) deallocate (forcecr)
      Allocate (forcecr(3, natmtot))
      If (allocated(forceibs)) deallocate (forceibs)
      Allocate (forceibs(3, natmtot))
!      If (allocated(forcetot)) deallocate (forcetot)
!      Allocate (forcetot(3, natmtot))
!      If (allocated(forcetp)) deallocate (forcetp)
!      Allocate (forcetp(3, natmtot))
!      If (allocated(tauatm)) deallocate (tauatm)
!      Allocate (tauatm(natmtot))
!      If (allocated(tauxyz)) deallocate (tauxyz)
!      Allocate (tauxyz(3, natmtot))
! initialise the previous force
!      forcetp (:, :) = 0.d0
! initial step sizes
!      If (associated(input%relax)) Then
!         tauatm (:) = input%relax%taunewton
!         tauxyz (:, :) = input%relax%taunewton
!      Else
!         tauatm (:) = 0
!         tauxyz (:, :) = 0
!      End If
!
!-------------------------!
!     LDA+U variables     !
!-------------------------!
      If ((ldapu .Ne. 0) .Or. (task .Eq. 17)) Then
! LDA+U requires second-variational eigenvectors
         input%groundstate%tevecsv = .True.
! density matrices
         If (allocated(dmatlu)) deallocate (dmatlu)
         Allocate (dmatlu(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
! potential matrix elements
         If (allocated(vmatlu)) deallocate (vmatlu)
         Allocate (vmatlu(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
! zero the potential
         vmatlu (:, :, :, :, :) = 0.d0
! energy for each atom
         If (allocated(engyalu)) deallocate (engyalu)
         Allocate (engyalu(natmtot))
! interpolation constants (alpha)
         If (allocated(alphalu)) deallocate (alphalu)
         Allocate (alphalu(natmtot))
      End If
!
!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine the nuclear-nuclear energy
      Call energynn
! get smearing function data
      Call getsdata (input%groundstate%stypenumber, sdescr)
! generate the spherical harmonic transform (SHT) matrices
      Call genshtmat
      Call genshtmat3
      ! Call genshtmat2
!
! allocate 1D plotting arrays
      If (allocated(dvp1d)) deallocate (dvp1d)
      Allocate (dvp1d(nvp1d))
      If (allocated(vplp1d)) deallocate (vplp1d)
      Allocate (vplp1d(3, npp1d))
      If (allocated(dpp1d)) deallocate (dpp1d)
      Allocate (dpp1d(npp1d))
!

! initialisation for the Davidson solver
      nsingular=-1
      mine0=input%groundstate%solver%minenergy

      Call timesec (ts1)
!!      timeinit = timeinit + ts1 - ts0
!
      call stopwatch("exciting:init0", 0)
      Return
End Subroutine
!EOC
