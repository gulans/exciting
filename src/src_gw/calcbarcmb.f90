!BOP
!
!!ROUTINE: calcbarcmb
!
!!INTERFACE:
!
subroutine calcbarcmb(iq)
!
!!DESCRIPTION:
!
!This subroutine calculates the matrix of the bare coulomb potential
!
!!USES:
    use modinput
    use modgw
    use mod_coulomb_potential
    use mod_mpi_gw, only: myrank
    use modmain
!    use iso_fortran_env, only: int64, real64
!!INPUT PARAMETERS:
    implicit none
    integer, intent(in) :: iq ! index of the q-point

!!LOCAL VARIABLES:
    integer :: imix, jmix, igq, jgq
    real(8) :: tstart, tend, t0, t1!, t2, t3
    character(len=256) :: filename
    complex(8), allocatable :: barc_lr(:,:)
    ! for diagonalization subroutine
    real(8) :: vl, vu, abstol
    integer :: il, iu, neval, info, i, j
    integer :: liwork, lrwork, lwork, matsiz11
    complex(8), allocatable :: work(:)
!#ifdef zheevrr
    complex(8), allocatable :: v_copy(:,:)
!#endif
    real(8),    allocatable :: rwork(:)
    integer,    allocatable :: iwork(:), ifail(:), isuppz(:)
    real(8), external :: dlamch

!!REVISION HISTORY:
!
! Created Jan 2014 by DIN
!
!EOP
!BOC
    call timesec(tstart)

!===============================================================================
! Setup the bare Coulomb potential matrix in MB representation
!===============================================================================

    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,matsiz))
    barc(:,:) = 0.d0

    select case (trim(input%gw%barecoul%basis))

    case('pw')

      call calcmpwmix(iq)
      call calcbarcmb_pw(iq)

    case('mb')

      if (Gamma) then
        !------------------------------------------------
        ! Matrix elements for the singular q=0, L=0 case
        !------------------------------------------------
        call barcq0
      end if

      !-----------------------------------------------------------
      ! Matrix elements between MT and MT mixed product functions
      !-----------------------------------------------------------
      call calcbarcmb_mt_mt(iq)

      !-----------------------------------------------------------
      ! Matrix elements between an atomic mixed function and an IPW
      !-----------------------------------------------------------
      call calcbarcmb_ipw_mt(iq)

      !-----------------------------------------------------------
      ! Matrix elements between two IPW's
      !-----------------------------------------------------------
      call calcbarcmb_ipw_ipw(iq)

      !----------------------------------------
      ! Only in case of HSE self-consitent run
      !----------------------------------------
      if ((task == 7) .and. (xctype(1)==408)) then
        ! Compute the short-range Coulomb potential
        allocate(barc_lr(matsiz,matsiz))
        call calcbarcmb_lr(iq, barc_lr)
        barc(:,:) = barc(:,:) - barc_lr(:,:)
        deallocate(barc_lr)
      endif

    case default

      write(*,*) 'ERROR(calcbarcmb): Unknown basis type!'
      stop

    end select


!===============================================================================
! Diagonalize the bare coulomb matrix
!===============================================================================
write(*,*)"start diagonalization, calcbarcmb"
    if (allocated(vmat)) deallocate(vmat)
    allocate(vmat(matsiz,matsiz))
    vmat(:,:) = barc(:,:)
    deallocate(barc)

    if (allocated(barcev)) deallocate(barcev)
    allocate(barcev(matsiz))

#ifdef zheevrr
    lrwork = -1
    liwork = -1
    lwork = -1
    iu = matsiz
    matsiz11=matsiz
    !abstol = 2.d0*dlamch('S')
    abstol=1d-12
    allocate(work(1),rwork(1),iwork(1), isuppz(2*matsiz))
    write(*,*)"before first"
    allocate(v_copy(matsiz, matsiz))
    v_copy=0.d0
    barcev=0d0
    vl=0
    vu=1d10
    il=1
    iu=matsiz11
    !call zheevd('V', 'U', matsiz, vmat, matsiz, barcev, work, lwork, rwork, lrwork, iwork, liwork, info)
    !call zheevr('V', 'A', 'U', matsiz11, vmat, matsiz11, vl, vu, il, iu, abstol, matsiz11, barcev, v_copy, matsiz11, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
!    call errmsg(info.ne.0, 'CALCBARCMB', "Fail to diag. barc by zheevd !!!")
write(*,*)"rwork", rwork(1)
    lrwork= 230880 !int(rwork(1))*5
    liwork= 96200  !int(iwork(1))*5
    lwork=  625300  !int(work(1))*5
    write(*,*) lrwork,liwork,lwork
    deallocate(work,rwork,iwork)
write(*,*)"previous error", info
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    write(*,*)"lwork", lwork, "mmatsiz", matsiz
    write(*,*) lrwork,liwork,lwork
 
    !do i = 1,10
    !    do j=  1,10
    !      write(*,*)i,j,"------    ", vmat(i,j)
     !   end do
    !end do
    !call timesec(t2)
    !call zheevd('V', 'U', matsiz, vmat, matsiz, barcev, work, lwork, rwork, lrwork, iwork, liwork, info)
    call zheevr('V', 'A', 'U', matsiz11, vmat, matsiz11, vl, vu, il, iu, abstol, matsiz11, barcev, v_copy, matsiz11, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

    !call timesec(t3)
    !write(*,*)"Diagonalization time :", t3-t2
    vmat=v_copy
    deallocate(v_copy, isuppz)
    write(*,*)"after zheevrrrrrr"
    call errmsg(info.ne.0, 'CALCBARCMB', "Fail to diag. barc by zheevd !!!")
    deallocate(work,rwork,iwork)
#endif
#ifdef bsesol
write(*,*)"BSE subroutine"
allocate(v_copy(matsiz, matsiz))
v_copy=0.d0
matsiz11=matsiz
call bsesoldiag(matsiz11, matsiz11, vmat, barcev, v_copy)
vmat=v_copy
deallocate(v_copy)
#endif

#ifdef zheevdd
 
    lrwork = -1
    liwork = -1
    lwork = -1
    !iu = matsiz
    !abstol = 2.d0*dlamch('S')
    !abstol=0.d0
    allocate(work(1),rwork(1),iwork(1))!, isuppz(2*matsiz))
    write(*,*)"before first"
    !allocate(v_copy(matsiz, matsiz))
    !v_copy=0.d0
    call zheevd('V', 'U', matsiz, vmat, matsiz, barcev, work, lwork, rwork, lrwork, iwork, liwork, info)
    !call zheevr('V', 'A', 'U', matsiz, vmat, matsiz, vl, vu, il, iu, abstol, matsiz,&
    !       & barcev, v_copy, matsiz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    call errmsg(info.ne.0, 'CALCBARCMB', "Fail to diag. barc by zheevd !!!")
write(*,*)"rwork", rwork(1)
    lrwork=int(rwork(1))
    liwork=int(iwork(1))
    lwork=int(work(1))
    write(*,*) lrwork,liwork,lwork
    deallocate(work,rwork,iwork)
write(*,*)"previous error", info
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    write(*,*)"lwork", lwork, "mmatsiz", matsiz
    write(*,*) lrwork,liwork,lwork
 
    !do i = 1,10
    !    do j=  1,10
    !      write(*,*)i,j,"------    ", vmat(i,j)
     !   end do
    !end do
    !call timesec(t2)
    call zheevd('V', 'U', matsiz, vmat, matsiz, barcev, work, lwork, rwork, lrwork, iwork, liwork, info)
    !call zheevr('V', 'A', 'U', matsiz, vmat, matsiz, vl, vu, il, iu, abstol, matsiz,&
    !       & barcev, v_copy, matsiz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    !call timesec(t3)
    !write(*,*)"Diagonalization time :", t3-t2
    !vmat=v_copy
    !deallocate(v_copy, isuppz)
    write(*,*)"after zheevd"
    call errmsg(info.ne.0, 'CALCBARCMB', "Fail to diag. barc by zheevd !!!")
    deallocate(work,rwork,iwork)
#endif
    write(*,*)"all done calcbarcmb",info

    !stop
!----------------------
! debug info
!----------------------
!      do imix = 1, matsiz
!        write(fdebug,'(i5,e16.6)') imix, barcev(imix)
!      end do
!stop

    if (input%gw%debug) then
      msize = sizeof(barcev)*b2mb+sizeof(vmat)*b2mb
      write(fdebug,'("calcbarcmb: rank, size(Coulomb potential) (Mb):",i4,f12.2)') myrank, msize
      write(fdebug,*) "### barcev ###"
      do imix = 1, matsiz
        write(fdebug,'(i5,e16.6)') imix, barcev(imix)
      end do
      write(fdebug,*) "### vmat ###"
      do imix = 1, matsiz, matsiz/10
        do jmix = 1, matsiz, matsiz/10
          write(fdebug,'(2i5,2e16.6)') imix, jmix, vmat(imix,jmix)
        end do
      end do
    endif !debug

    call timesec(tend)
    time_barcmb = time_barcmb+tend-tstart

end subroutine
!EOC
