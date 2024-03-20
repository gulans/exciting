!> This module contains routines to generate all files needed 
!> by the cc4s code (https://manuals.cc4s.org/user-manual/index.html)
module cc4s_interface

  use precision, only: sp, dp
  use cc4s_interface_output, only: coulomb_potential, coulomb_potential_spencer_alavi

  implicit none

contains

  !> Launches the calculation of the discretized Coulomb Vertex
  !> \( \tilde{\Gamma}^q_{rG'} \). 
  !> See routine 'calculate_coulomb_vertex_g_basis' for documentation.
  !> When referring to equation number, we mean the publication 
  !> Hummel et al., J. Chem. Phys. 146, 124105 (2017).
  subroutine cc4s_interface_main(input)
    ! !USES:
    use modmpi
    use modinput, only: input_type
    use mod_APW_LO, only: lolmax
    use mod_kpoint, only: nkpt
    use mod_qpoint, only: nqpt
    use mod_lattice, only: bvec
    use modxs, only: xsgnt, nwdf, qpari, ngqmax,&
                     & qparf, ngq, unitout, totalqlmt, qvkloff,&
                     & gqdirname, eps0dirname, scrdirname, &
                     & ikmapikq, nkpt0, vkl0, usefilext0, filext0, filexteps,&
                     & iqmt0, iqmt1, evalsv0, istocc0, istunocc0, isto0, isto,&
                     & istu0, istu, nst1, nst2, ngq,  vgql, vgqc
    use mod_xsgrids
    use m_genfilname
    use m_filedel
    use m_writegqpts
    use m_xsgauntgen
    use m_findgntn0
    use mod_Gkvector, only: gkmax
    use m_ematqk
    use mod_eigenvalue_occupancy, only: nstsv, evalsv, efermi
    use sorting, only: sortidx
    use asserts, only: assert

    use cc4s_interface_output

    type(input_type), intent(in) :: input
    !> Filename of Coulomb Vertex (in G-vector or auxiliary field basis)
    character(len=200) :: coulomb_vertex_filename_bin = "CoulombVertex.elements"
    !> File handle for coulomb vertex
    integer :: coulomb_vertex_filehandle
    !> Running indices q -, (k+q) -, k - points 
    integer :: iq, i_kq_counter, ik_plus_q, ik
    !> Indices of Gamma-point in sets of k-/ q-vectors. exciting sets them always to 1
    integer, parameter :: ik_gamma = 1, iq_gamma = 1 ! TODO (Max):  should be made more robust
    !> File extension
    character(256) :: filex
    !> System command
    character(256) :: syscommand
    !> Routine name for INFOXS.OUT
    character(*), parameter :: thisname = 'cc4s_interface_main'
    !> Indices of k -, q -, (k+q)-points
    integer, allocatable :: momentum_triple_list(:, :)
    !> Coulomb vertex for all G-vectors, all states, and one (k,q) combination.
    complex(dp), allocatable :: coulomb_vertex_g_basis(:, :, :)
    !> Map from ik, iq to combined index
    integer, allocatable :: curr_kq_pair(:, :)
    !> Number of spins used in calculation
    integer, parameter :: CC4S_n_spin = 1
    !> Fraction of number auxiliary field vectors to number of G-vectors  <= 1 
    real(dp), parameter :: fraction_gvecs_to_auxfield = .7_dp
    !> Number of auxiliary field vectors used for construction of Coulomb Vertex
    integer :: n_auxfield
    !> True if Coulomb Vertex should be transformed to opt. aux. field basis 
    logical :: do_optimized_vertex
    !> True if a Gamma-point-only (k=q=k+q=0) calculation should be performed
    logical :: do_gamma_only
    !> True if tests should be performed
    logical :: do_tests
    !> Number of momentum grid points (= number of G-vectors or number of
    !> auxiliary field vectors) used for construction of the Coulomb vertex 
    integer :: n_basis_coulomb_vertex
    !> Number of G-vectors used for construction of Coulomb Vertex 
    !> (should be equal for all q-vectors)
    integer :: n_gvecs
    !> First and last index of (k,q)-pairs the current process computes 
    integer :: kq_start, kq_end
    !> Index of (k,q)-combination
    integer :: ikq
    !> Eigenvalues of "squared" Coulomb vertex in G-vector basis 
    real(dp),  allocatable :: eigvals_squared_coulomb(:)
    !> Left eigenvectors of "squared" Coulomb vertex in G-vector basis 
    complex(dp), allocatable :: eigvecs_squared_coulomb_l(:, :)
    !> Coulomb vertex in auxiliary field basis for all states and one (k,q) combination.
    complex(dp), allocatable :: coulomb_vertex_auxfield_basis(:, :, :)
    ! DFT eigenvalues for all k-points and all spins
    real(dp), allocatable :: scf_eigenvalues(:,:,:)
    !> Flattened (= 1 dimensional) version of DFT eigenvalues for all k-points and all spins
    real(dp), allocatable :: scf_eigenvalues_flattened(:)
    !> Indices of sorted (in ascending order) DFT eigenvalues
    integer, allocatable :: sorted_ids_scf_eigvals(:)
    !> Number of spins in calculation
    integer, parameter :: n_spin = 1
    !> CC4S version number - could be changed?
    character(3), parameter :: cc4s_version = '100' 

    do_optimized_vertex = input%xs%cc4s%dooptimizedvertex
    do_gamma_only = input%xs%cc4s%dogammaonly
    do_tests = input%xs%cc4s%dotests


    ! Initialise universal variables
    call init0
    ! Setting up k and G+k variables
    call init1
    ! Save k and G+k grid variables to
    ! modxs (vkl0, ngk0, ...)
    call xssave0
    ! Set q-points and (G+q)-points
    call init2

    ! Set file extension to '_SCR', needed to read EIGVAL_SCR.OUT and EFERMI_SCR.OUT
    call genfilname(dotext='_SCR.OUT', setfilext=.true.)

    ! Read fermy energy from EFERMI_SCR.OUT 
    call readfermi

    allocate(scf_eigenvalues_flattened(nstsv))
    allocate(scf_eigenvalues(nstsv, n_spin, nkpt))

    ! Read DFT eigenvalues (generated by scrgeneigvec) for all k-points
    do ik = 1, nkpt
      call getevalsv(vkl(:,ik), scf_eigenvalues(:,1,ik))
    end do

    ! Flatten DFT eigenvalues array, needed for cc4s output
    scf_eigenvalues_flattened = pack( scf_eigenvalues, .true.)

    ! Get indices of sorted DFT eigenvalues
    allocate(sorted_ids_scf_eigvals(nkpt*nstsv*n_spin))
    call sortidx(nkpt*nstsv*n_spin, scf_eigenvalues_flattened, sorted_ids_scf_eigvals)

    ! Not sure if this is needed... Clear modxs:ecalsv0
    if (allocated(evalsv0)) deallocate (evalsv0) 
    
    ! Number of G-vectors (should be equivalent for all q-vectors)
    ! by setting gmaxtype = |G| in input.xml
    n_gvecs = ngq(1) 

    call assert(ngqmax == n_gvecs, message='Number of G-vectors varies over q-vectors')

    if(do_optimized_vertex) then
        n_auxfield = int(ngq(1) * fraction_gvecs_to_auxfield)
        n_basis_coulomb_vertex = n_auxfield
    else
        n_basis_coulomb_vertex = n_gvecs
        n_auxfield = n_gvecs
    end if

    ! stores indices in order [ik,iq,i(k+q)]
    allocate (momentum_triple_list(3, nkpt*nqpt))

    ! Generate map of (k,q)-pairs to combined index
    allocate (curr_kq_pair(nqpt, nkpt))
    
    i_kq_counter = 0
    do iq = 1, nqpt
      do ik = 1, nkpt
        ik_plus_q = ikmapikq(ik, iq)
        i_kq_counter = i_kq_counter + 1
        momentum_triple_list(:, i_kq_counter) = (/ik, iq, ik_plus_q/)
        curr_kq_pair(iq, ik) = i_kq_counter
      end do
    end do


    call barrier(callername=trim(thisname))

    !------------------------------------------------------------!
    ! Preparation for plane wave matrix elements                 !
    !------------------------------------------------------------!
    ! Generate gaunt coefficients, and store them in modxs:xsgnt
    call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
      & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
    ! Find indices for non-zero gaunt coefficients, and store
    ! relevant maps in the module m_findgntn0
    call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
      & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
    !------------------------------------------------------------!

    ! Read Fermi energy from file EFERMI
    ! Use EFERMI_SCR.OUT (corresponding to the scr groundstate run
    ! for the reference k grid)
    call genfilname(scrtype='', setfilext=.true.)
    call readfermi

    !! Set parameters for the plane-wave matrix element calculation,
    !! e.g. which files to use for bra and ket states.
    call findocclims(iq_gamma, ikmapikq(:, iq_gamma), istocc0, istunocc0, isto0, isto, istu0, istu)
    call ematbdcmbs(input%xs%cc4s%bandcombinations)
    
    ! TODO (Max): Check if setting extensions is needed or was already done before
    ! Set *_SCR.OUT as bra state file
    usefilext0 = .true.
    iqmt0 = 1
    call genfilname(scrtype='', fileext=filext0)
    
    ! Set *_SCR.OUT as ket state file
    iqmt1 = 1
    call genfilname(scrtype='', setfilext=.true.)
    
    ! Use <mk|e^{-i(G+q)r}|nk'> for q=k'-k in dfq
    emat_ccket = .false.

    ! Set *.OUT as file extension for screening files
    call genfilname(fileext=filexteps)

    ! Write all interesting information to file
    if(mpiglobal%is_root) then

      ! QPOINTS_SCR.OUT
      call genfilname(scrtype='', setfilext=.true.)
      call writeqpts

      ! INFOXS.OUT output and command line 
      call write_infoxs_output(n_gvecs,n_auxfield,do_optimized_vertex,do_gamma_only,unitout, thisname, nqpt, istunocc0,istocc0, nst1, nst2)
      call write_cmdline_output(do_gamma_only,do_optimized_vertex,ngq(iq_gamma),n_auxfield, 1,1,istunocc0,istocc0, nst1, nst2)

      ! Making folder for GQPOINTS info files
      gqdirname = 'GQPOINTS'
      if (rank == 0) then
        syscommand = 'test ! -e '//trim(adjustl(gqdirname))//' && mkdir '//trim(adjustl(gqdirname))
        call system(trim(adjustl(syscommand)))
      end if

      ! Write (G+q)-vectors to files in directory GQPOINTS/
      do iq = 1, nqpt
        ! Write q-point number to fileext, filext = "_SCR_QXYZ.OUT"
        call genfilname(scrtype='', iq=iq, fileext=filex)

        ! Write out G+q vectors to file "GQPOINTS_SCR_QXYZ.OUT"
        call writegqpts(iq, filex, dirname=gqdirname)
      end do

      ! Files needed by cc4s
      call CC4S_out_eigenenergies_yaml(do_gamma_only, cc4s_version, 1._dp, efermi, scf_eigenvalues_flattened(sorted_ids_scf_eigvals),1,1)
      ! TODO (Max): Check if lowest state is really always 1
      call CC4S_output_scf_energies(scf_eigenvalues,nstsv,n_spin,nkpt,1, .false.)
      call CC4S_out_orbital_properties_yaml(sorted_ids_scf_eigvals,nkpt,nstsv,n_spin)
      call CC4S_out_coulvertex_yaml(do_gamma_only, cc4s_version, nst1, 1, n_basis_coulomb_vertex, momentum_triple_list, 1._dp, .true., .false.)
      call CC4S_out_spin_yaml(cc4s_version, nstsv,n_spin,nkpt,1._dp)

      ! These ones are probably not needed at the moment
      call CC4S_out_momentum_properties_yaml(ngq)
      call CC4S_out_grid_vectors_yaml(bvec, vgql)
      call CC4S_out_grid_vectors_element(vgql, do_gamma_only)
      call CC4S_out_coulomb_potential_yaml(vgql, do_gamma_only)
      call CC4S_out_coulomb_potential_element(vgqc)
    end if

    ! Open file for Coulomb Vertex on all participating processes
#ifdef MPI
    call MPI_File_open(mpiglobal%comm, coulomb_vertex_filename_bin, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                &MPI_INFO_NULL, coulomb_vertex_filehandle, ierr)
#endif

    allocate(coulomb_vertex_g_basis(n_gvecs, nst1, nst2))
    allocate(coulomb_vertex_auxfield_basis(n_auxfield,nst1, nst2))


    if (do_gamma_only) then

      if(mpiglobal%is_root) then

        call calculate_coulomb_vertex_g_basis(iq_gamma, ik_gamma, coulomb_vertex_g_basis)
        
        if(do_tests) then
          write(*, *) 
          write(*, *) '------Running tests------'
          call calc_coulomb_integral_1s(coulomb_vertex_g_basis)

          call calc_exact_exchange_energy(istocc0, coulomb_vertex_g_basis)
          write(*, *) '------Finished tests------'
          write(*, *) 

        end if

        if (do_optimized_vertex) then
          call calculate_coulomb_vertex_auxfield_basis(coulomb_vertex_g_basis, n_auxfield, eigvals_squared_coulomb, eigvecs_squared_coulomb_l, coulomb_vertex_auxfield_basis)
        end if

      ! Write in parallel to file
#ifdef MPI
        call CC4S_parallel_writing_cmplx(coulomb_vertex_g_basis, coulomb_vertex_filehandle,  curr_kq_pair(iq_gamma, ik_gamma),&
        &CC4S_n_spin, nst1,0,1)
        write (*, *) 'Vertex written for iq, ik: ', iq_gamma, ik_gamma
#endif
 
      end if


  else

    ! TODO: Need to distribute
    if(mpiglobal%is_root) then
      ! Calculate Coulomb-Vertex
      do iq = 1, nqpt
        do ik = 1, nkpt

        ! Loop over k-points
          ! Compute plane-wave matrix element for given q-point
          call calculate_coulomb_vertex_g_basis(iq, ik, coulomb_vertex_g_basis)

        ! Write in parallel to file
#ifdef MPI
        call CC4S_parallel_writing_cmplx(coulomb_vertex_g_basis, coulomb_vertex_filehandle, curr_kq_pair(iq, ik), CC4S_n_spin,&
        &nst1,0,1)
        write (*, *) 'Vertex written for iq, ik: ', iq, ik
#endif

        end do
      end do
    end if
   
  end if
  

#ifdef MPI
    call MPI_File_close(coulomb_vertex_filehandle, ierr)
#endif
 
  if(mpiglobal%is_root) then 
      write(*,*) 
      write (*,*) '--------- cc4s-Interface: Calculations finished ---------'
      write(*,*) 
  end if
    ! Synchronize
    call barrier(callername=trim(thisname))
    !------------------------------------------------------------!

    ! Delete gaunt maps
    call findgntn0_clear

  end subroutine cc4s_interface_main

  !> Calculates  the discretized Coulomb Vertex
  !> \( \tilde{\Gamma}^q_{rG'} \) , defined as
  !>
  !> \[ 
  !>    \tilde{\Gamma}^q_{rG'} = \sqrt{w_{G'}}\Gamma^q_r(\mathbf{G'}_{G'})
  !>                                                                \]
  !>
  !> with the "full" Coulomb Vertex, \(  \Gamma^q_r(\mathbf{G'}) \)
  !>
  !> \[
  !>     \Gamma^q_r(\mathbf{G'}) = 
  !>        \sqrt{\frac{4\pi}{|\mathbf{G'}|^2}}
  !>        \int d\mathbf{x}e^{-i\mathbf{G'\cdot r}}
  !>        \psi^{*q}(\mathbf{x})\psi_r(\mathbf{x}).
  !>                                                  \]
  !>
  !> The integration weights \( w_{G'}  \) are given by the inverse
  !> number of momentum grid points according to \(  w_{G'} = 1 / N_{G'}   \)
  !> (compare J. Chem. Phys. 146, 124105 (2017), Eq. 6) .
  !>
  !> In contrast to the reference, we denote the general reciprocal 
  !> momentum by \( \mathbf{G'}  \), which in practice is expressed by 
  !> \( \mathbf{G'} = \mathbf{q+G}  \), 
  !> where \( \mathbf{q}  \) is a vector from the first Brilllouin zone and 
  !> \( \mathbf{G}  \) is a reciprocal lattice vector. The momentum vector 
  !> \( \mathbf{q}  \) corresponds to the momentum transfer between 
  !> the two electronic states (see following paragraph.)
  !>
  !> The electron orbitals are labeled by \( q = m\mathbf{k}\sigma_q \)
  !> and \( r = n\mathbf{k+q}\sigma_r \) with the momentum transfer
  !> \( \mathbf{q}  \) being the BZ fraction of grid point \( \mathbf{G'} \).
  !> 
  !> The "full" Coulomb vertex is constructed using the matrix elements
  !> of the planewave (which are calculated by existing exciting routines),
  !>
  !> \[
  !>    M^{\mathbf{G}}_{mn}(\mathbf{k,q}) =
  !>     \langle m\mathbf{k}  
  !>              | e^{-i(\mathbf{G+q})\cdot\mathbf{r}} |
  !>                               n\mathbf{k+q} \rangle
  !>                                                          \]
  !> according to
  !>  
  !> \[ 
  !>      \Gamma^{m\mathbf{k}}_{n\mathbf{k+q}}(\mathbf{G+q}) = 
  !>        \sqrt{\frac{4\pi}{|\mathbf{G+q}|^2}}
  !>                       M^{\mathbf{G}}_{mn}(\mathbf{k,q})
  !>                          
  !>                                                        \]
  subroutine calculate_coulomb_vertex_g_basis(iq, ik, coulomb_vertex_g_basis)
! !USES:
    use mod_misc, only: filext
    use modinput, only: input
    use constants, only: zzero, zone, zi
    use mod_kpoint, only: nkpt, wkpt
    use mod_qpoint, only: nqpt, vql, vqc
    use mod_lattice, only: omega
    use modxs, only:ngq,  istocc0, istocc,&
                  & istunocc0, istunocc, isto0, isto,&
                  & istu0, istu, unitout, nst1,&
                  & nst2, nst3, nst4, istl1,&
                  & istu1, istl2, istu2, istl3,&
                  & istu3, istl4, istu4, &
                  &ikmapikq, &
                  & bcbs, qvkloff,fnemat, &
                  & eps0dirname, scrdirname, vgqc
#ifdef TETRA
    use modxs, only: fnwtet
    use modtetra
#endif
    use m_genwgrid
    use m_getpemat
    use m_dftim
    use m_gettetcw
    use m_putx0
    use m_getunit
    use m_writevars
    use m_filedel
    use m_genfilname
    use m_ematqk
    use m_writecmplxparts
    use putgeteps0, only: puteps0_finite_q, puteps0_zero_q
    use mod_variation, only: ematqk_sv
  
    implicit none

    ! Arguments
    integer, intent(in) :: iq
    integer, intent(in) :: ik

    ! Local variables
    character(*), parameter :: thisname = 'calculate_coulomb_vertex_g_basis'
    character(256) :: fnscreen, fneps0, filex
    real(8), parameter :: epstetra = 1.d-8

    integer(4) :: n, j, i1, i2, ikq, igq, iw, wi, wf, ist1, ist2, nwdfp
    integer(4) :: un
    logical :: tq0
    logical ::fintraband
    type(bcbs) :: band_limits
    complex(dp), allocatable :: pw_matrix_elements(:, :, :)
    complex(dp), intent(out) :: coulomb_vertex_g_basis(:, :, :)
    character(256) :: file_coul_vertex

    ! External functions
    logical, external :: tqgamma

    ! Matrix size for response function
    ! Get number of G+q vectors for current q
    n = ngq(iq)

    ! Set whether intraband should be used
    fintraband = input%xs%tddft%intraband .or. input%xs%screening%intraband
    ! File extension for q-point (not in 'screen')

    ! Calculate k+q and G+k+q related variables
    ! by setting the offset generated from vkloff and the q point
    ! and then calling init1
    call init1offs(qvkloff(1:3, iq))

    ! Find limits for band combinations
    call findocclims(iq, ikmapikq(:, iq), istocc0, istunocc0, isto0, isto, istu0, istu)
    istunocc = istunocc0
    istocc = istocc0

    tq0 = tqgamma(iq)

    call ematrad(iq)
    call ematqalloc
    call genfilname(basename='EMAT', iq=iq, filnam=fnemat)
    call genfilname(basename='coulomb_vertex', iq=iq, filnam=file_coul_vertex)

    ikq = ikmapikq(ik, iq)

    ! The plane wave elements for ou and uo transitions are
    ! calculated and stored in xiou and xiuo
    ! Set 12=ou 34=uo
    call ematbdcmbs(input%xs%cc4s%bandcombinations)

    ! Get ou
    if (allocated(pw_matrix_elements)) deallocate (pw_matrix_elements)
    allocate (pw_matrix_elements(nst1, nst2, n))


    band_limits%n1 = nst1
    band_limits%n2 = nst2
    band_limits%il1 = istl1
    band_limits%il2 = istl2
    band_limits%iu1 = istu1
    band_limits%iu2 = istu2
    ikmapikq_ptr => ikmapikq
    call setptr01
    if (.not. input%groundstate%tevecsv) then
      call ematqk(iq, ik, pw_matrix_elements, band_limits)
    else
      call ematqk_sv(iq, ik, pw_matrix_elements, band_limits)
    end if

    call make_coulomb_vertex_g_basis(omega, pw_matrix_elements, vgqc(:,:,iq), coulomb_vertex_g_basis)

    ! Write to file
    call write_planewave_elements(iq, ik, vql(:, iq), vkl(:, ik), &
                                  trim(fnemat), band_limits, &
                                  pw_matrix_elements)

  end subroutine

  !> Transforms the discretized Coulomb Vertex in the G-vector basis 
  !> \( \tilde{\Gamma}^q_{rG'} \) to the  
  !> optimized auxiliary field basis. i.e.  \( \tilde{\Gamma}^q_{rF} \)
  !> where \( F \) labels the auxiliary field momenta. 
  !> Note that this implementation is currently limited 
  !> to the Gamma-only case,  i.e.,
  !>  \( q = m,\mathbf{k=0}, r = n,\mathbf{k=0}  \).
  !> 
  !> The transformation is computed by 
  !>  ( with the compound index \(  I = (q,r) \) )
  !>
  !> \[
  !>    \Gamma_{FI} = U^*_{FG}\tilde{\Gamma}_{GI}, 
  !>                                            \]
  !>
  !> where the matrix stems from the singular value decomposition
  !>
  !> \[
  !>  \tilde{\Gamma}_{GI} = U_{GF}\Sigma_{FJ}W^*_{JI}
  !>                                                    ]
  !> To avoid a brute-force SVD, which for large marices becomes unfeasible, 
  !> we proceed in the following way:
  !> First, the "squared" Coulomb vertex \( E_{GG'} \) is computed according to
  !> Eq. 10:
  !>
  !> \[
  !>    E_{GG'} = \tilde{\Gamma}_{GI}\tilde{\Gamma}^*_{IG'}  
  !>                                                          \]
  !> 
  !> Second, \( E \) is diagonalised. The left eigenvectors associated 
  !>  to the largest eigenvalues are also the left singular vectors of 
  !> \( \tilde{\Gamma}_{GI} \), i.e. the matrix \(  U_{GF}  \) and the 
  !> eigenvalues are the squares of the singular values of 
  !> \( \tilde{\Gamma}_{GI} \), i.e. the matrix \(  \Sigma_{FJ}  \).
  !> 
  !> By keeping only the elements corresponding to the largest \( N_F \)
  !> (this number is given as input) singular values, the size of the
  !> Coulomb vertex can be significantly reduced.
  !> 
  !> ONLY WORKING FOR GAMMA-ONLY CALCS!!!
  subroutine calculate_coulomb_vertex_auxfield_basis(coulomb_vertex_g_basis, n_auxfield, eigvals_squared_coulomb, eigvecs_squared_coulomb_l,coulomb_vertex_auxfield_basis)
    use tensor_contractions, only: complex_tensor_contraction_dp
    use m_diagfull, only: diagfull
    use math_utils, only: all_close, identity_complex_dp

     complex(dp), intent(in) :: coulomb_vertex_g_basis(:, :, :)
     complex(dp), intent(out) :: coulomb_vertex_auxfield_basis(:, :, :)

     complex(dp), allocatable :: left_singular_vectors(:, :)
    integer, intent(in) :: n_auxfield
    
    complex(dp), allocatable :: squared_coulomb_vertex_g_basis(:,:)
     integer :: n_gvecs, ig, n_states
     real(dp), intent(out), allocatable :: eigvals_squared_coulomb(:)
     complex(dp),intent(out), allocatable :: eigvecs_squared_coulomb_l(:, :)
     complex(dp), allocatable :: eigvecs_squared_coulomb_r(:, :)

     n_states = size(coulomb_vertex_g_basis,dim=2)

     n_gvecs = size(coulomb_vertex_g_basis, dim = 1)

     allocate(squared_coulomb_vertex_g_basis(n_gvecs, n_gvecs))

     call complex_tensor_contraction_dp(coulomb_vertex_g_basis,&
                                      & shape(coulomb_vertex_g_basis),&
                                      coulomb_vertex_g_basis,&
                                      & shape(coulomb_vertex_g_basis),&
                                      & squared_coulomb_vertex_g_basis,&
                                      shape(squared_coulomb_vertex_g_basis),&
                                      trans_A = 'n', trans_B ='c')

      allocate(eigvals_squared_coulomb(n_gvecs))
      allocate(eigvecs_squared_coulomb_l(n_gvecs, n_gvecs))
      allocate(eigvecs_squared_coulomb_r(n_gvecs, n_gvecs))

      call diagfull(n_gvecs,squared_coulomb_vertex_g_basis,eigvals_squared_coulomb,&
                    evecl=eigvecs_squared_coulomb_l, fsort=.true.)

      allocate(left_singular_vectors(n_gvecs, n_auxfield))

      ! Keep only larges N_F auxiliary field vectors
      left_singular_vectors = eigvecs_squared_coulomb_l(:, n_gvecs - n_auxfield + 1:n_gvecs)

      call complex_tensor_contraction_dp(left_singular_vectors,&
                                        shape(left_singular_vectors),&
                                        coulomb_vertex_g_basis,&
                                        shape(coulomb_vertex_g_basis),&
                                        coulomb_vertex_auxfield_basis,&
                                        shape(coulomb_vertex_auxfield_basis),&
                                        trans_A='c')

  end subroutine calculate_coulomb_vertex_auxfield_basis

  subroutine make_coulomb_vertex_g_basis(unit_cell_volume,pw_matrix_elements, vgqc, coulomb_vertex_g_basis)
    use constants, only: fourpi, twopi
    use asserts, only: assert

    real(dp), intent(in) :: vgqc(:,:)
    complex(dp), intent(in) :: pw_matrix_elements(:, :, :)
    complex(dp), intent(out) :: coulomb_vertex_g_basis(:, :, :)
    real(dp), intent(in) :: unit_cell_volume
    real(dp) :: brillouin_zone_volume

    !> number of G-vectors
    integer :: n_gvecs
    integer :: ig
    real(dp) :: pref
    real(dp) :: square_root_ngvecs
    real(dp) :: integration_weight

    !> Shift that makes G finite for G = 0
    real(dp), parameter :: g_shift = 1e-6_dp

    n_gvecs = size(vgqc, dim=2)
    square_root_ngvecs = sqrt(dble(n_gvecs))

    brillouin_zone_volume = twopi**3._dp / unit_cell_volume

    integration_weight = brillouin_zone_volume / (twopi**3._dp)
    write(*, *) 'Integration weight:', integration_weight

    call assert(size(coulomb_vertex_g_basis,dim=1) == n_gvecs, message='Size of Coulomb vertex not correct')


    do ig = 1, n_gvecs

      ! pref =  sqrt(coulomb_potential(vgqc(:,ig)))
      pref =  sqrt(coulomb_potential_spencer_alavi(unit_cell_volume, vgqc(:,ig)))
      coulomb_vertex_g_basis(ig, :, :) = pref*pw_matrix_elements(:, :, ig)
    end do

    coulomb_vertex_g_basis = sqrt(integration_weight) * coulomb_vertex_g_basis

    write(*, *) "Convergence calculated:", coulomb_potential_spencer_alavi(unit_cell_volume,vgqc(:,1))

  end subroutine

  !> Writes the matrix elements of the planewave, defined by
  !>
  !> \[
  !>   M_{mn}^{\mathbf G}(\mathbf{k,q}) =
  !>        \langle m\mathbf{k}|e^{i\mathbf{(G+q)r)}}|n\mathbf{k+q}\rangle
  !>                                                        \],
  !> for a given (q,k)-combination  from file.
  subroutine write_planewave_elements(iq, ik, q_vec, k_vec, &
                                      filename, band_limits, pw_matrix_elements)
    Use modmpi
    Use m_getunit
    use modxs, only: bcbs

    !> q-vector index
    Integer, intent(in)  :: iq
    !> k-vector index
    Integer, intent(in)  ::  ik
    !> q-vector in lattice coordinates
    real(dp), intent(in) :: q_vec(3)
    !> k-vector in lattice coordinates
    real(dp), intent(in) :: k_vec(3)
    !> Number of (G+q)-vectors
    integer:: n_gqvecs
    !> Filename
    Character(*), intent(in)  :: filename
    !> Band intervals
    class(bcbs), intent(in)  :: band_limits
    !> Matrix elements of the plane for one (k,q)-combination
    complex(dp), intent(in) :: pw_matrix_elements(:, :, :)
    !> Fileunit
    integer :: un
    !> Record length
    integer ::  recl
    !> Lowest band for index \( m \)
    integer :: lower_limit_m
    !> Lowest band for index \( n \)
    integer :: lower_limit_n
    !> Highest band for index \( m \)
    integer :: upper_limit_m
    !> Highest band for index \( n \)
    integer :: upper_limit_n

    lower_limit_m = band_limits%il1
    lower_limit_n = band_limits%il2
    upper_limit_m = band_limits%iu1
    upper_limit_n = band_limits%iu2

    Call getunit(un)

    n_gqvecs = size(pw_matrix_elements, dim=3)

    ! I/O record length
    Inquire (IoLength=Recl) q_vec, k_vec, &
    &  n_gqvecs, lower_limit_m, upper_limit_m, &
    lower_limit_n, upper_limit_n, pw_matrix_elements

    Open (Unit=un, File=trim(filename), Form='unformatted', &
    & Action='write', Access='direct', Recl=Recl)
    Write (un, Rec=ik) q_vec, k_vec, &
    &  n_gqvecs, lower_limit_m, upper_limit_m, &
    lower_limit_n, upper_limit_n, pw_matrix_elements

    close (un)
  end Subroutine write_planewave_elements

  !> Reads the matrix elements of the planewave, defined by
  !>
  !> \[
  !>   M_{mn}^{\mathbf G}(\mathbf{k,q}) =
  !>        \langle m\mathbf{k}|e^{i\mathbf{(G+q)r)}}|n\mathbf{k+q}\rangle
  !>                                                        \],
  !> for a given (q,k)-combination from file.
  subroutine read_planewave_elements(iq, ik, q_vec, k_vec, &
                                     n_second_variation, filename, band_limits, pw_matrix_elements)
    Use modmpi
    Use m_getunit
    use modxs, only: bcbs

    !> q-vector index
    Integer, intent(in)  :: iq
    !> k-vector index
    Integer, intent(in)  ::  ik
    !> q-vector in lattice coordinates
    real(dp), intent(in) :: q_vec(3)
    !> k-vector in lattice coordinates
    real(dp), intent(in) :: k_vec(3)
    !> Number of (G+q)-vectors
    integer :: n_gqvecs
    !> Number of (G+q)-vectors
    integer, intent(in) :: n_second_variation
    !> Filename
    Character(*), intent(in)  :: filename
    !> Band intervals
    class(bcbs), intent(in)  :: band_limits
    !> Matrix elements of the plane for one (k,q)-combination
    complex(dp), intent(out) :: pw_matrix_elements(:, :, :)
    !> Fileunit
    integer :: un
    !> Record length
    integer ::  recl
    !> Lowest band for index \( m \)
    integer :: lower_limit_m
    !> Lowest band for index \( n \)
    integer :: lower_limit_n
    !> Highest band for index \( m \)
    integer :: upper_limit_m
    !> Highest band for index \( n \)
    integer :: upper_limit_n

    ! Variables for sanity checks

    !> q-vector index stored in file
    Integer  :: iq_read
    !> k-vector index stored in file
    Integer  ::  ik_read
    !> q-vector in lattice coordinates stored in file
    real(dp) :: q_vec_read(3)
    !> k-vector in lattice coordinates stored in file
    real(dp) :: k_vec_read(3)
    !> Number of (G+q)-vectors stored in file
    integer :: n_gqvecs_read
    !> Number of (G+q)-vectors stored in file
    integer :: n_second_variation_read
    !> Lowest band for index \( m \) stored in file
    integer :: lower_limit_m_read
    !> Lowest band for index \( n \) stored in file
    integer :: lower_limit_n_read
    !> Highest band for index \( m \) stored in file
    integer :: upper_limit_m_read
    !> Highest band for index \( n \) stored in file
    integer :: upper_limit_n_read

    lower_limit_m = band_limits%il1
    lower_limit_n = band_limits%il2
    upper_limit_m = band_limits%iu1
    upper_limit_n = band_limits%iu2

    n_gqvecs = size(pw_matrix_elements, dim=3)

    Call getunit(un)

    ! I/O record length
    Inquire (IoLength=Recl) q_vec_read, k_vec_read, &
    & n_second_variation_read, n_gqvecs_read, lower_limit_m_read, upper_limit_m_read, &
    lower_limit_n_read, upper_limit_n_read, pw_matrix_elements

    Open (Unit=un, File=trim(filename), Form='unformatted', &
    & Action='read', Access='direct', Recl=Recl)
    read (un, Rec=ik) q_vec_read, k_vec_read, &
    & n_second_variation_read, n_gqvecs_read, lower_limit_m_read, upper_limit_m_read, &
    lower_limit_n_read, upper_limit_n_read, pw_matrix_elements
    close (un)

  end Subroutine read_planewave_elements

  ! Writes relevant information to INFOXS.OUT
  subroutine write_infoxs_output(n_gvecs,n_auxfield, do_optimized_vertex,do_gamma_only,unitout, thisname, nqpt, istunocc0,istocc0, nst1, nst2)

    integer, intent(in) :: unitout
    integer, intent(in) :: nqpt
    integer, intent(in) :: istunocc0
    integer, intent(in) :: istocc0
    integer, intent(in) :: nst1
    integer, intent(in) :: nst2
    !> True if Coulomb Vertex should be transformed to opt. aux. field basis 
    logical, intent(in) :: do_optimized_vertex
    !> True if a Gamma-point-only (k=q=k+q=0) calculation should be performed
    logical,intent(in) :: do_gamma_only
    integer, intent(in) :: n_gvecs
    integer, intent(in) :: n_auxfield
  
    character(*), intent(in) :: thisname

    call printline(unitout, "+")
    write (unitout, '(a, i8)') 'Info('//thisname//'):&
      & Calculating Coulomb vertex for unshifted q-grid.'
    write (unitout, '(a, i8)') 'Info('//thisname//'):&
      & Number of q points:', nqpt
    call printline(unitout, "+")
    write (unitout, *)

    if(do_gamma_only) then
      write(unitout,'(a, a8)') 'Gamma-only calculation:'
    else 
      write(unitout,'(a, a8)') 'All k-point calculation'
    end if

    if(do_optimized_vertex) then
      write(unitout,'(a, a8)') 'Coulomb Vertex in optimized auxiliary field basis'
    else
      write(unitout,'(a, a8)') 'Coulomb Vertex in G-vector basis'
    end if

       
    write (unitout, '(a, 4i6)') 'Info('//thisname//'):&
      & lowest (partially)  unoccupied state: ', istunocc0
    write (unitout, '(a, 4i6)') 'Info('//thisname//'):&
      & highest (partially) occupied state  : ', istocc0
    write (unitout, '(a, i5)') 'Info('//thisname//'):&
      & number of states in first dimension:', nst1
    write (unitout, '(a, i5)') 'Info('//thisname//'):&
      & number of states in second dimension:', nst2
      call printline(unitout, "-")
    write (unitout, '(a, 4i6)') 'Info('//thisname//'):&
      & lowest (partially)  unoccupied state: ', istunocc0
    write (unitout, '(a, i6)') 'Info('//thisname//'):&
      & Number of G-vectors: ', n_gvecs
    write (unitout, '(a, i6)') 'Info('//thisname//'):&
      & Number of auxiliary field vectors: ', n_auxfield
  end subroutine



  ! Writes relevant information to as command line output
  subroutine write_cmdline_output(do_gamma_only,do_optimized_vertex, n_gvecs, n_auxfield, nqpt,nkpt, istunocc0,istocc0, nst1, nst2)

    !> True if Coulomb Vertex should be transformed to opt. aux. field basis 
    logical :: do_optimized_vertex
    !> True if a Gamma-point-only (k=q=k+q=0) calculation should be performed
    logical :: do_gamma_only

    integer, intent(in) :: nqpt
    integer, intent(in) :: nkpt
    integer, intent(in) :: n_gvecs


    integer, intent(in) :: istunocc0
    integer, intent(in) :: istocc0
    integer, intent(in) :: nst1
    integer, intent(in) :: nst2
    integer, intent(in) :: n_auxfield
  
    write(*,*) 
    write (*,*) '--------- cc4s-Interface: Calculations started ---------'
    write(*,*) 

    write(*,*) 'Gamma-only calculation:', do_gamma_only
    write(*,*) 'Optimized auxiliary field calculation:', do_optimized_vertex

    write(*,*)  'Number of q points:', nqpt
    write(*,*)  'Number of k points:', nkpt
    write(*,*)  'Number of G vectors:', n_gvecs
    if(do_optimized_vertex) then
      write(*,*)  'Number of auxiliary-field vectors:', n_auxfield
    else 
    write(*,*)  'Number of auxiliary-field vectors:  None'
    end if
  
    write(*,*)  'lowest (partially)  unoccupied state: ', istunocc0
    write(*,*)  'highest (partially) occupied state  : ', istocc0
    write(*,*)  'number of states in first dimension:', nst1
    write(*,*)  'number of states in second dimension:', nst2


  end subroutine

  
  
  !> Computes the Coulomb integral <1s1s|V|1s1s> for the
  !> Helium atom. Expected result is around 1 Ha.
  subroutine calc_coulomb_integral_1s(coulomb_vertex)
    use constants, only: zzero
    use unit_conversion, only: hartree_to_ev

    complex(dp), intent(in) :: coulomb_vertex(:,:,:)
    integer :: ig, n_gvecs
    complex(dp), allocatable :: coulomb_integral_1s
    integer, parameter :: id_1s = 1
    integer, parameter :: print_steps = 500

    n_gvecs = size(coulomb_vertex, dim=1)
    
    
    coulomb_integral_1s = zzero

    do ig = 1, n_gvecs
      
      coulomb_integral_1s = coulomb_integral_1s + abs(coulomb_vertex(ig, id_1s, id_1s))**2
      
      if (mod(ig, print_steps) == 0) then
        write(*, *) 'G-vector index:', ig
        write(*, "(x,A40, F10.5)") 'Calculated value for <1s,1s|V|1s,1s> in Ha:', coulomb_integral_1s      
      end if
    end do
    

    write(*, *) 
    write(*, "(x,A50, F10.5)") 'Final: Calculated value for <1s,1s|V|1s,1s> in Ha:', real(coulomb_integral_1s)
    write(*, "(x,A50, F10.5)") 'Final: Calculated value for <1s,1s|V|1s,1s> in eV:', real(coulomb_integral_1s)*hartree_to_ev
    write(*, "(x,A50, F10.5)") 'Reference value VASP for <1s,1s|V|1s,1s> in eV:', 27.597764019596607_dp

  end subroutine
  subroutine calc_exact_exchange_energy(n_occupied,coulomb_vertex)
    use constants, only: zzero
    use unit_conversion, only: hartree_to_ev

    complex(dp), intent(in) :: coulomb_vertex(:,:,:)
    integer, intent(in) :: n_occupied
    integer :: ig, n_gvecs
    complex(dp), allocatable :: exact_exchange_energy
    integer :: i_occ, j_occ
    integer, parameter :: print_steps = 500

    n_gvecs = size(coulomb_vertex, dim=1)
    
    
    exact_exchange_energy = zzero

    do j_occ = 1, n_occupied
      do i_occ = 1, n_occupied
        do ig = 1, n_gvecs

          exact_exchange_energy = exact_exchange_energy + abs(coulomb_vertex(ig, i_occ, j_occ))**2

        end do
      end do
    end do

    exact_exchange_energy = -1._dp / 2._dp * exact_exchange_energy
    
    write(*, *) 
    write(*, *) 'Final: Calculated value for exact exchange energy:', real(exact_exchange_energy)

  end subroutine

end module cc4s_interface
