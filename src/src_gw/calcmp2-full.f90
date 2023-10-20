subroutine calcmp2_full(iq,iomstart,iomend)
!
! Compute the RPA dielectric matrix
!
  use precision, only: dp
#ifdef MPI
    use mpi
#endif MPI   
    use modinput
    use modmain, only : zone, zzero, zi
    use modgw
    use mod_mpi_gw, only : myrank, mycomm
    use modxs,      only : symt2
    use m_getunit
    use mp2_energy
    use sorting, only: sortidx
    use cc4s_interface_output

    implicit none
    ! input/output
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: iomstart, iomend
    ! local
    integer(4) :: ie1, ie2,ie3,ie4
    integer(4) :: iom
    integer(4) :: ik, jk, ispn
    integer(4) :: im, iop, jop
    integer(4) :: ndim, mdim, nmdim,nocc,nunocc, n_states
    integer(4) :: nblk, iblk, mstart, mend
    integer(8) :: recl
    real(8)    :: tstart, tend
    real(8)    :: wto, wlo
    complex(8) :: head(3,3), f, w
    complex(8), allocatable :: minm(:,:,:),Viabj(:,:,:,:),minm2(:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    ! Fermi energy
    real(8) :: fermi_energy 
    ! DFT eigenvalues for all k-points and all spins
    real(8), allocatable :: scf_eigenvalues(:,:,:)
    !> Flattened (= 1 dimensional) version of DFT eigenvalues for all k-points and all spins
    real(dp), allocatable :: scf_eigenvalues_flattened(:)
    !> Indices of sorted (in ascending order) DFT eigenvalues
    integer, allocatable :: sorted_ids_scf_eigvals(:)
    !> Number of spins in calculation
    integer, parameter :: n_spin = 1, n_kpt = 1
    !> CC4S version number - could be changed?
    character(3), parameter :: cc4s_version = '100'
    integer, allocatable :: momentum_triple_list(:, :)
    logical, parameter :: do_gamma_only = .true. ! Change when implementing no only gamma
    
    !> Filename of Coulomb Vertex (in G-vector or auxiliary field basis)
    character(len=200) :: coulomb_vertex_filename_bin = "CoulombVertex.elements"
    !> File handle for coulomb vertex
    integer :: coulomb_vertex_filehandle
    
    integer :: ierr


    call timesec(tstart)

    

    !=============================
    ! Initialization
    !=============================


    ndim = nstdf !nomax
    mdim = nstdf !nstdf-numin+1
    nmdim = ndim*mdim

    nocc=nomax
    nunocc=nstdf-nocc

    write(*,*) 'nocc,nunocc',nocc,nunocc

    ! block size
    if (mblksiz >= mdim) then
      nblk = 1
    else
      nblk = mdim / mblksiz
      if ( mod(mdim,mblksiz) /= 0 ) nblk = nblk+1
    end if

    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))

    !==================================================
    ! Calculate the q-dependent BZ integration weights
    !==================================================
    select case (trim(input%gw%qdepw))
    case('sum')
        call qdepwsum(iq, iomstart, iomend, ndim)
    case('tet')
        call qdepwtet(iq, iomstart, iomend, ndim)
    case default
        stop "Error(calcepsilon): Unknown qdepw method!"
    end select
    
    ! Open file for Coulomb Vertex on all participating processes
#ifdef MPI
    call MPI_File_open(mycomm, coulomb_vertex_filename_bin, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                &MPI_INFO_NULL, coulomb_vertex_filehandle, ierr)
#endif


    !=================
    ! BZ integration
    !=================
    ! write(*,*)
    do ik = 1, kqset%nkpt
        ! k-q point
        jk = kqset%kqid(ik, iq)

!        if (Gamma) then
            ! read the momentum matrix elements
!            call getpmatkgw(ik)
            ! and compute the head of the dielectric function
!            call calchead(ik, iomstart, iomend, ndim)
!        end if

        ! get KS eigenvectors
        allocate(evecfv(nmatmax,nstfv))
        call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
        eveckp = conjg(evecfv)
        call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
        eveck = evecfv
        deallocate(evecfv)

        ! compute products \sum_G C_{k}n * A_{lm}
        call expand_evec(ik,'t')
        call expand_evec(jk,'c')

        !=================================================
        ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
        !=================================================
!        do iblk = 1, nblk

            ! call timesec(ta)

            mstart = 1 !numin + (iblk-1)*mblksiz
            mend   = nstdf !min(nstdf, mstart+mblksiz-1)
            nmdim  = ndim * (mend-mstart+1)
            ! print*, iblk, nblk, mstart, mend
!
!  Coulomb vertex dimensions 
!  mixed basis size  X  number of all bands  X  number of all bands 
            allocate(minmmat(mbsiz,nstdf,nstdf))
            msize = sizeof(minmmat)*b2mb

            ! compute M^i_{nm}+M^i_{cm}
write(*,*) 'expand_products',nstdf
            call expand_products(ik, iq, 1, nstdf, nstdf, 1, nstdf, -1, minmmat)
write(*,*) 'expand_products done'
!            call expand_products(ik, iq, 1, mend, nstdf, 1, mend, -1, minmmat)
write(*,*)sum(minmmat(:, 1, 1), 1)
!            if (Gamma) then
                ! wings of the dielectric matrix
!                call calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend)
!            end if

            ! Body



 
if (.true.) then
  !> gamma_ph = gamma(:, N_occ + 1:, :N_occ)
  !> gamma_hp = gamma(:, :N_occ, N_occ + 1:)
            allocate(minm(mbsiz,1:nunocc,1:nocc))
            allocate(minm2(mbsiz,1:nocc,1:nunocc))

            do ie1=1,nocc
              do ie2=1,nunocc
                minm(:,ie2,ie1)=minmmat(:,ie2+nocc,ie1)
              enddo
            enddo
            do ie1=1,nunocc
              do ie2=1,nocc
                minm2(:,ie2,ie1)=minmmat(:,ie2,ie1+nocc)
              enddo
            enddo

            !write(*,*) 'E_MP2',calculate_mp2_energy(minm, minm2, evalfv(1:nocc,ik), evalfv(nocc+1:nocc+nunocc,ik))
             call calculate_mp2_energy(minm, minm2, evalfv(1:nocc,ik), evalfv(nocc+1:nocc+nunocc,ik))
    endif
      ! Flatten DFT eigenvalues array, needed for cc4s output
      n_states = nocc + nunocc
      allocate(scf_eigenvalues(n_states, n_spin, n_kpt))
      scf_eigenvalues(:, 1, 1) = evalfv(:, 1)
      scf_eigenvalues_flattened = pack( scf_eigenvalues, .true.)
      fermi_energy = evalfv(nocc, 1) + (evalfv(nocc+1, 1) - evalfv(nocc, 1)) / 2 ! Take the middle in the band gap. Fine for cc4s 

      ! Get indices of sorted DFT eigenvalues
      allocate(sorted_ids_scf_eigvals(size(scf_eigenvalues)))
      call sortidx(size(scf_eigenvalues), scf_eigenvalues_flattened, sorted_ids_scf_eigvals)

      allocate(momentum_triple_list(3, n_kpt ** 2)) ! only needed for more than one k-point (not implemented yet)

      call CC4S_out_eigenenergies_yaml(do_gamma_only, cc4s_version, 1._dp, fermi_energy, scf_eigenvalues_flattened(sorted_ids_scf_eigvals), 1, 1)
      ! TODO (Max): Check if lowest state is really always 1
      call CC4S_output_scf_energies(scf_eigenvalues, n_states, n_spin, n_kpt, 1, .false.)
      call CC4S_out_orbital_properties_yaml(sorted_ids_scf_eigvals, n_kpt, n_states, n_spin)
      call CC4S_out_coulvertex_yaml(do_gamma_only, cc4s_version, n_states, 1, mbsiz, momentum_triple_list, 1._dp, .true., .false.)

             ! Write in parallel to file
#ifdef MPI
     call CC4S_parallel_writing_cmplx(minmmat, coulomb_vertex_filehandle, 1, n_spin, n_states)
#endif

end do

#ifdef MPI
    call MPI_File_close(coulomb_vertex_filehandle, ierr)
#endif

end subroutine
