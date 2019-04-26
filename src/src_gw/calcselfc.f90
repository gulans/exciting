subroutine calcselfc(iq)
    use modinput
    use modmain,    only : nstfv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                      zzero, nmatmax
    use modgw
    use mod_mpi_gw, only : myrank
    use m_getunit
    implicit none
    ! input/output
    integer(4), intent(in) :: iq
    ! local
    integer(4) :: ik, ikp, jk, ispn
    integer(4) :: mdim, iblk, nblk, mstart, mend
    integer(4) :: fid
    character(120) :: fname_mwm
    real(8) :: tstart, tend, t0, t1
    complex(8), allocatable :: evecfv(:,:)
    character(len=10), external :: int2str

    call timesec(tstart)

    !------------------------
    ! total number of states
    !------------------------
    if (input%gw%coreflag=='all') then
      mdim = nstse+ncg
    else
      mdim = nstse
    end if

    !-------------------------------------------------------
    ! determine the number of blocks used in minm operation
    !-------------------------------------------------------
    if (mblksiz >= mdim) then
      nblk = 1
    else
      nblk = mdim / mblksiz
      if (mod(mdim,mblksiz) /= 0) nblk = nblk+1
    end if

    !-------------------------------------------
    ! products M*W^c*M
    !-------------------------------------------
    allocate(mwm(ibgw:nbgw,1:mdim,1:freq%nomeg))
    ! msize = sizeof(mwm)*b2mb
    ! write(*,'(" calcselfc: size(mwm) (Mb):",f12.2)') msize

    !----------------------------
    ! q-dependent M*W*M products
    !----------------------------
    if (input%gw%taskname == 'gw0') then
      fname_mwm = 'MWM'//'-q'//trim(int2str(iq))//'.OUT'
      call getunit(fid)
      open(fid,File=fname_mwm,Action='Write',Form='Unformatted')
    end if

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))

    !================================
    ! loop over irreducible k-points
    !================================
    write(*,*)
    do ikp = 1, kset%nkpt
      write(*,*) 'calcselfc: rank, (iq, ikp):', myrank, iq, ikp

      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector
      jk = kqset%kqid(ik,iq)

      ! get KS eigenvectors
      allocate(evecfv(nmatmax,nstfv))
      call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
      eveckp = conjg(evecfv)
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
      eveck = evecfv
      deallocate(evecfv)

      call expand_evec(ik, 't')
      call expand_evec(jk, 'c')

      !=================================
      ! Loop over m-blocks in M^i_{nm}
      !=================================
      do iblk = 1, nblk

        mstart = 1 + (iblk-1)*mblksiz
        mend = min(mdim, mstart+mblksiz-1)

        ! m-block M^i_{nm}
        allocate(minmmat(mbsiz,ibgw:nbgw,mstart:mend))
        msize = sizeof(minmmat)*b2mb
        ! write(*,'(" calcselfc: rank, size(minmmat) (Mb):",3i4,f12.2)') myrank, mstart, mend, msize

        call expand_products(ik, iq, ibgw, nbgw, -1, mstart, mend, nstse, minmmat)

        !================================================================
        ! Calculate weight(q)*Sum_ij{M^i*W^c_{ij}(k,q;\omega)*conjg(M^j)}
        !================================================================
        call calcmwm(ibgw, nbgw, mstart, mend, minmmat)

        deallocate(minmmat)

      end do ! iblk

      !=======================================
      ! Calculate the correlation self-energy
      !=======================================
      if (input%gw%taskname=='cohsex') then
        call calcselfc_cohsex(ikp, iq, mdim)
      else
        if (input%gw%selfenergy%method == 'cd') then
          ! Contour deformation technique
          call calcselfc_freqconv_cd(ikp, iq, mdim)
        else if (input%gw%selfenergy%method == 'ac') then
          ! Imaginary frequency formalism
          call calcselfc_freqconv_ac(ikp, iq, mdim)
        end if
      end if

      if (input%gw%taskname=='gw0') then
        ! store M*W*M in files
        call timesec(t0)
        write(fid) mwm
        call timesec(t1)
        time_io = time_io+t1-t0
      end if

    end do ! ikp

    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    ! delete MWM
    deallocate(mwm)
    ! and close the file
    if (input%gw%taskname.eq.'gw0') close(fid)

    ! timing
    call timesec(tend)
    time_selfc = time_selfc+tend-tstart

    return
end subroutine
