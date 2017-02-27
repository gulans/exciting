
subroutine getevalqp(nkp2,kvecs2,eqp2)

  use modinput
  use modmain
  use modgw,    only: kset, ibgw, nbgw, nkp1, kvecs1, eks1, eqp1, eferqp
  use mod_wannier

  implicit none
      
  integer, intent(in) :: nkp2
  real(8), intent(in) :: kvecs2(3,nkp2)
  real(8), intent(out):: eqp2(nstsv,nkp2)

  logical       :: exist
  integer(4)    :: ik, ib, nb, nk, nqp
  integer(4)    :: recl
  real(8)       :: eferks, kvecs2_( 3, nkp2)
  character(30) :: fname
  integer(4), allocatable :: idx(:)
  real(8),    allocatable :: eqp(:), eqpwan(:,:), eqpwanint(:,:)
  complex(8), allocatable :: de1(:,:), de2(:,:)

  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      

  fname='EVALQP.OUT'
  inquire(File=fname, Exist=exist)
  if (.not.exist) then
    write(*,*)'ERROR(getevalqp): File EVALQP.OUT does not exist!'
    stop
  end if
      
  inquire(IoLength=recl) nkp1, ibgw, nbgw
  open(70, File=fname, Action='READ', Form='UNFORMATTED', &
  &    Access='DIRECT', Recl=recl)
  read(70, Rec=1) nkp1, ibgw, nbgw
  close(70)
      
  allocate(kvecs1(1:3,nkp1))
  allocate(eqp1(ibgw:nbgw,nkp1))
  allocate(eks1(ibgw:nbgw,nkp1))
  
  ! old format (gwmod-boron)    
  ! inquire(IoLength=recl) nkp1, ibgw, nbgw, kvecs1(1:3,1), &
  ! &       eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1)      

  ! new format (carbon)
  inquire(IoLength=recl) nkp1, ibgw, nbgw, kvecs1(1:3,1), &
  &       eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1), &
  &       eferqp, eferks
  
  open(70, File=fname, Action='READ', Form='UNFORMATTED', &
  &    Access='DIRECT', Recl=recl)
  
  nqp = nbgw-ibgw+1
  allocate(idx(nqp))
  allocate(eqp(nqp))

  do ik = 1, nkp1

    ! old format (gwmod-boron) 
    ! read(70, Rec=ik) nk, ib, nb, kvecs1(:,ik), &
    ! &    eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik)

    ! new format (carbon)
    read(70, Rec=ik) nk, ib, nb, kvecs1(:,ik), &
    &    eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik), &
    &    eferqp, eferks

    !write(fgw,*) '# ik    kvecs1    ibgw,    nbgw'
    !write(fgw,*) ik, kvecs1(:,ik), ib, nb
    !write(fgw,*) '# ib    eqp1    eks1'
    !do ib = ibgw, nbgw
    !  write(fgw,*) ib, eqp1(ib,ik), eks1(ib,ik)
    !end do          
  end do ! ik
  close(70)
  deallocate(idx)

  !------------------------------
  ! Data-set consistency check
  !------------------------------
  if ((ibgw /= 1).or.(nbgw < nstsv)) then
    write(*,*)
    write(*,*)'WARNING(getevalqp):'
    write(*,*)'  Quasiparticle energies has been calculated for the interval of bands:'
    write(*,*)'  [ibgw,nbgw]=[', ibgw, nbgw,']'
    write(*,*)'  Check whether it is that you need ...'
    write(*,*)
  end if

  !----------------------------------------------
  ! Special case of only one k-point (molecules)
  !----------------------------------------------
  if (nkp1==1) then
    if (nkp2==1) then
      do ib = ibgw, min(nbgw,nstsv)
        eqp2(ib,1) = eqp1(ib,1)
      end do
      deallocate(kvecs1,eqp1,eks1)
      return
    else
      write(*,*) 'ERROR(getevalqp):' 
      write(*,*) '  Interpolation is not possible!'
      write(*,*) '  EVALQP.OUT file contains data only for a single k-point.'
      stop
    end if
  end if 

  !-----------------------------------------------------------------------------
  ! Interpolate the energies
  !-----------------------------------------------------------------------------      
  allocate(de1( nkp1, ibgw:nbgw))

  allocate(de2( nkp2, ibgw:nbgw))
  de2(:,:) = zzero

  if( input%gw%taskname .eq. "wannier") then
    call readfermi
    kvecs2_ = kvecs2
    allocate( eqpwan( ibgw:nbgw, wf_nkpt), eqpwanint( ibgw:nbgw, nkp2))
    if( allocated( vkl)) deallocate( vkl)
    allocate( vkl( 3, nkp1))
    vkl = kvecs1
    do ik = 1, wf_nkpt
      call findkpt( wf_vkl( :, ik), nb, ib)
      !write(*,*) ik, ib
      !write(*,*) wf_vkl( :, ik)
      !write(*,*) vkl( :, ib)
      !write(*,*) kvecs1( :, ib)
      eqpwan( :, ik) = eqp1( ibgw:nbgw, ib)
      !write(*,'(100F13.6)') eqpwan( :, ik)
    end do
    deallocate( vkl)
    allocate( vkl( 3, nkp2))
    vkl = kvecs2_
    nkpt = nkp2
      
    !write(*,*) shape( eqpwan), shape( wf_vkl), wf_nkpt
    !write(*,*) shape( eqpwanint), shape( kvecs2), nkp2
    !write(*,*) shape( eqp2)
    call wannier_interpolate_eval( eqpwan, wf_nkpt, wf_vkl, eqpwanint, nkp2, kvecs2_, ibgw, nbgw)
    do ib = ibgw, min(nbgw,nstsv)
       do ik = 1, nkp2
          eqp2( ib, ik) = eqpwanint( ib, ik) - eferqp - efermi
       end do 
     end do
  else
    do ik = 1, nkp1
      de1(ik,:) = cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
    end do
    call fourintp(de1, nkp1, kvecs1, de2, nkp2, vkl, nbgw-ibgw+1)
    do ib = ibgw, min(nbgw,nstsv)
       do ik = 1, nkpt
          eqp2(ib,ik) = eqp2(ib,ik)+dble(de2(ik,ib))-eferqp
       end do 
     end do
  end if

  deallocate(de1,de2)
  deallocate(kvecs1,eqp1,eks1)

  return
end subroutine
