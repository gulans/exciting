module modsurf

  use mod_atoms, only: natmtot
implicit none

integer              :: naxismax
integer,allocatable  :: naxis(:) !ias
integer,allocatable  :: axisy(:,:),axisz(:,:) !iax,ias
real(8),allocatable  :: axisx_sph(:,:) !iax,ias
complex(8),allocatable :: ylm_mat(:,:,:) !lm,iax,ias
complex(8),allocatable :: ylm_tmat(:,:,:) !lm,iax,ias
contains

subroutine generate_surf_grid(lmax)
  use invert, only: zinvert_lapack
  use m_linalg, only: zlsp

  use modinput
  use mod_Gvector, only: ngrtot, ngrid, ngvec, cfunir
  use mod_atoms, only: nspecies, natoms, atposc, idxas,natmtot
  use mod_muffin_tin, only: rmt
  use constants, only: zzero
  implicit none
  
  real(8) :: a1(3),a2(3),a3(3),r0(3),rv(3),rv1(3)
  real(8) :: xmu
  real(8) :: ratom(3), disc
  integer, intent(In) :: lmax
  integer :: lmmax
  integer :: is,ia,ias,yi,zi,p1,p2,p3
  integer :: nax,iax
  real(8),allocatable  :: axisx_sph_tmp(:,:), tp_tmp(:,:,:) ! 2,iax,ias
  integer,allocatable  :: axisy_tmp(:,:),axisz_tmp(:,:)
  complex(8),allocatable  :: A(:,:),Ainv(:,:),AA(:,:)
  
  

  lmmax=(lmax+1)**2
  write(*,*)"Režģa izmērs",ngrid
  write(*,*)"x ass garums",sqrt(sum(input%structure%crystal%basevect(1, :)**2))


  a1=input%structure%crystal%basevect(1, :)/ngrid(1)
  a2=input%structure%crystal%basevect(2, :)/ngrid(2)
  a3=input%structure%crystal%basevect(3, :)/ngrid(3)

  Allocate(axisy_tmp(ngrid(2)*ngrid(3),natmtot),axisz_tmp(ngrid(2)*ngrid(3),natmtot),axisx_sph_tmp(ngrid(2)*ngrid(3),natmtot))
  allocate(naxis(natmtot))
  allocate(tp_tmp(2,ngrid(2)*ngrid(3),natmtot))
  ia=1
  is=1
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      nax=0
      do p1=-1,1
        do p2=-1,1
          do p3=-1,1
            ratom=atposc (:, ia, is) + p1*input%structure%crystal%basevect(1, :) + p2*input%structure%crystal%basevect(2, :) + p3*input%structure%crystal%basevect(3, :)
            do yi=0,ngrid(2)-1
              do zi=0,ngrid(3)-1

                r0 = yi*a2 + zi*a3
                rv=r0-ratom

                disc= sum(a1*rv)**2 - sum(a1**2)*(sum(rv**2)-rmt(is)**2)
              
                if (disc.gt.1d-9) then
                  xmu= ( -sum(a1*rv) - sqrt(disc) )/ sum(a1**2)
                  if ((xmu.ge.0d0).and.(xmu.lt.(ngrid(1)-1))) then
                    nax = nax+1
                    axisx_sph_tmp(nax,ias)=xmu
                    axisy_tmp(nax,ias)=yi
                    axisz_tmp(nax,ias)=zi
                    rv1=r0 + a1*xmu - ratom

                    !write(*,*)"rmt1",sqrt(sum((rv1)**2)) 
                    tp_tmp(1,nax,ias) = atan2(sqrt(rv1(1)**2+rv1(2)**2),rv1(3))
                    tp_tmp(2,nax,ias) = atan2(rv1(2),rv1(1))
 
                  endif

                  xmu= ( -sum(a1*rv) + sqrt(disc) )/ sum(a1**2)
                  if ((xmu.ge.0d0).and.(xmu.lt.(ngrid(1)-1))) then
                    nax = nax+1
                    axisx_sph_tmp(nax,ias)=xmu
                    axisy_tmp(nax,ias)=yi
                    axisz_tmp(nax,ias)=zi
                    rv1=r0 + a1*xmu - ratom

                    !write(*,*)"rmt2",sqrt(sum((rv1)**2)) 
                    tp_tmp(1,nax,ias) = atan2(sqrt(rv1(1)**2+rv1(2)**2),rv1(3))
                    tp_tmp(2,nax,ias) = atan2(rv1(2),rv1(1))
                  endif
                  
                endif
              enddo !zi
            enddo !yi
          enddo !p3
        enddo !p2
      enddo !p1
      naxis(ias)=nax
    enddo !ia
  enddo !is
  naxismax=maxval(naxis(:))
write(*,*)"maxasis",naxismax
write(*,*)naxis

Allocate(axisy(naxismax,natmtot),axisz(naxismax,natmtot),axisx_sph(naxismax,natmtot),ylm_mat(lmmax,naxismax,natmtot),ylm_tmat(lmmax,naxismax,natmtot))
ylm_mat=zzero
axisy=0
axisz=0
axisx_sph=0d0
do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    axisy(1:naxis(ias),ias)=axisy_tmp(1:naxis(ias),ias)
    axisz(1:naxis(ias),ias)=axisz_tmp(1:naxis(ias),ias)
    axisx_sph(1:naxis(ias),ias)=axisx_sph_tmp(1:naxis(ias),ias)
    do iax=1, naxis(ias)
      Call genylm (lmax, tp_tmp(:, iax,ias), ylm_mat(:,iax,ias))
    enddo
  enddo
enddo

! do is=1, nspecies
!   do ia=1, natoms(is)
!     ias=idxas(ia,is)
!     write(*,*)"shapes"
!     write(*,*)shape(ylm_mat(:,:,ias))
!     write(*,*)"ylm"
!     write(*,*)ylm_mat(1,1:3,ias)
!     write(*,*)ylm_mat(2,1:3,ias)
!     write(*,*)ylm_mat(3,1:3,ias)

!     allocate( A(naxis(ias),naxis(ias)), Ainv(naxis(ias),naxis(ias)), AA(naxis(ias),naxis(ias)))

!     write(*,*)"shape1",shape( matmul(conjg(transpose(ylm_mat(:,1:naxis(ias),ias))),ylm_mat(:,1:naxis(ias),ias)) )
!     write(*,*)"shape2", shape(A)
!     A=matmul(conjg(transpose(ylm_mat(:,1:naxis(ias),ias))),ylm_mat(:,1:naxis(ias),ias))
!     call zinvert_lapack(A,Ainv)
!     AA=matmul(Ainv,A)
!     write(*,*)"pirmais A"
!     write(*,*)A(1,1:3)
!     write(*,*)A(2,1:3)
!     write(*,*)A(3,1:3)
!     write(*,*)"pirmais Ainv"
!     write(*,*)Ainv(1,1:3)
!     write(*,*)Ainv(2,1:3)
!     write(*,*)Ainv(3,1:3)
!     write(*,*)"pirmais Ainv*A"
!     write(*,*)AA(1,1:3)
!     write(*,*)AA(2,1:3)
!     write(*,*)AA(3,1:3)

!     A=matmul(conjg(transpose(ylm_mat(:,1:naxis(ias),ias))),ylm_mat(:,1:naxis(ias),ias))
!     call zinver(A,Ainv,naxis(ias))
!     AA=matmul(Ainv,A)
!     write(*,*)"otrs"
!     write(*,*)A(1,1:3)
!     write(*,*)A(2,1:3)
!     write(*,*)A(3,1:3)
!     write(*,*)"otrais Ainv"
!     write(*,*)Ainv(1,1:3)
!     write(*,*)Ainv(2,1:3)
!     write(*,*)Ainv(3,1:3)
!     write(*,*)"otrais Ainv*A"
!     write(*,*)AA(1,1:3)
!     write(*,*)AA(2,1:3)
!     write(*,*)AA(3,1:3)

! stop

!     call zinver(A,Ainv,naxis(ias))
!     write(*,*)shape(matmul(Ainv,transpose(ylm_mat(:,1:naxis(ias),ias))))



!     call zinver(Ainv,AA,naxis(ias))

!     write(*,*)A(1:3,1:3)
!     write(*,*)AA(1:3,1:3)
!     write(*,*)"stop" 
!     stop
!   enddo
! enddo



deallocate(axisy_tmp,axisz_tmp,axisx_sph_tmp,tp_tmp)
open(11,file='surf.dat',status='replace')
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      do iax=1, naxis(ias)
        rv(:)=axisx_sph(iax,ias)*a1 + axisy(iax,ias)*a2 + axisz(iax,ias)*a3
        write(11,*)rv(1),rv(2),rv(3)
      enddo
    enddo
  enddo
  close(11)



end subroutine 

subroutine surf_pot(lmax,zvclir,igfft,vlm)
  use modinput
  use m_linalg, only:zlsp
  use mod_Gvector, only: ngrtot, ngrid, ngvec, cfunir,vgc
  use mod_atoms, only: nspecies, natoms, atposc, idxas,natmtot
  use mod_muffin_tin, only: rmt
  use constants, only: zzero,pi
  implicit none
  integer  , intent(in) :: lmax
  complex(8) ,intent(inout) :: zvclir(:)
  integer, intent(in) :: igfft(:)
  complex(8) ,intent(inout) :: vlm(:,:)!lm,ias
  complex(8) :: zvax(ngrid(1)), zvaxft(ngrid(1))
  !complex(8) :: zvmu((lmax+1)**2)
  integer :: ig1,ig2,iax,ir,is,ia,ias,iy,iz,imu, nn
  integer :: lmmax,lm
  complex(8),allocatable :: vmu(:,:),vmu2(:,:)
  real(8) :: xx, irf,ir2,v1(3),a1(3),a2(3),a3(3),t1
  complex(8) ::zt1, ii
  integer :: ifg,ig
 
  complex(8) :: vlm2((lmax+1)**2,1)
  ii=cmplx(0d0,1d0,8)
  lmmax=(lmax+1)**2
  ias=1

  a1=input%structure%crystal%basevect(1, :)/ngrid(1)
  a2=input%structure%crystal%basevect(2, :)/ngrid(2)
  a3=input%structure%crystal%basevect(3, :)/ngrid(3)


  allocate(vmu(naxismax,natmtot))
  if(.false.)then
    !Tests ar sadalīšanu pa v_lm komponentēm
    ias=1
    vmu(1:naxis(ias),ias) = matmul(transpose(ylm_mat(:,1:naxis(ias),ias)),vlm(:,ias))
    open(11,file='vmu0.dat',status='replace')
    do iax=1,naxis(ias)
      write(11,*)dble(vmu(iax,ias)),imag(vmu(iax,ias))
    enddo
    close(11)
    allocate(vmu2(naxis(ias),1))
    vmu2(:,1)=vmu(1:naxis(ias),ias) 
    call zlsp(transpose(ylm_mat(:,1:naxis(ias),ias)), vmu2, vlm2) 
    open(11,file='lm-set.dat',status='replace')
    do lm=1,lmmax
      write(11,*)dble(vlm(lm,ias)),imag(vlm(lm,ias)),dble(vlm2(lm,1)),imag(vlm2(lm,1))
    enddo
    close(11)
  endif


if(.false.)then
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      do iax=1,naxis(ias)
        iz=axisz(iax,ias)
        iy=axisy(iax,ias)
        xx=axisx_sph(iax,ias)
        v1=xx*a1 + iy*a2 + iz*a3
        zt1=zzero
        Do ig = 1, ngvec
          ifg = igfft (ig)
          t1 = vgc(1,ig)*v1(1) + vgc(2,ig)*v1(2) + vgc(3,ig)*v1(3)
          zt1 = zt1 + zvclir(ifg)*cmplx(Cos(t1), Sin(t1), 8)
        End Do
        vmu(iax,ias)=zt1
      enddo !iax
      open(11,file='vmu3d.dat',status='replace')
      do iax=1,naxis(ias)
        write(11,*)dble(vmu(iax,ias)),imag(vmu(iax,ias))
      enddo
      close(11)
    enddo!ia
  enddo!is
endif





  
  
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      
      do iax=1,naxis(ias)
        iz=axisz(iax,ias)
        iy=axisy(iax,ias)

        ig1 = iz*ngrid(2)*ngrid(1) + iy*ngrid(1) + 1
        ig2 = iz*ngrid(2)*ngrid(1) + iy*ngrid(1) + ngrid(1)
        zvax(:)=zvclir(ig1:ig2)
        !write(*,*)"krustpunkts ar sferu",axisx_sph(iax,ias)
        zvaxft=zvax
        call cfftnd(1,ngrid(1),-1,zvaxft)

    !Fourier interpolation of potential on x=axisx_sph(iax,ias)
        xx=axisx_sph(iax,ias)
        zt1=zzero
        do ir=1, ceiling(dble(ngrid(1))/2d0)
          ir2=ir-1
          !write(*,*)ir,",",ir2
          irf=2d0*pi*dble(ir2)/dble(ngrid(1))
          zt1=zt1 + zvaxft(ir) * exp (ii*irf*xx)
        enddo
        do ir=ngrid(1),ceiling((dble(ngrid(1))-1)/2)+2,-1
          ir2=ngrid(1)-ir+1
          !write(*,*)ir,",",ir2
          irf=2d0*pi*dble(-ir2)/dble(ngrid(1))
          zt1=zt1 + zvaxft(ir) * exp (ii*irf*xx)
        enddo
        if (mod(ngrid(1),2).eq.0)then
          ir=ngrid(1)/2+1
          ir2=ngrid(1)/2
        ! write(*,*)"vidus",ir,ir2
          irf=2d0*pi*dble(ir2)/dble(ngrid(1))
          zt1=zt1 + zvaxft(ir) * cos (irf*xx)
        endif
        vmu(iax,ias)=zt1
      enddo !iax
      ! open(11,file='vmu1.dat',status='replace')
      ! do iax=1,naxis(ias)
      !   write(11,*)dble(vmu(iax,ias)),imag(vmu(iax,ias))
      ! enddo
      ! close(11)

    !Create a set of v_lm from v_mu
      allocate(vmu2(naxis(ias),1))
      vmu2(:,1)=vmu(1:naxis(ias),ias)
      call zlsp(transpose(ylm_mat(:,1:naxis(ias),ias)), vmu2, vlm2) 
      
      ! open(11,file='lm-set.dat',status='replace')
      !  do lm=1,lmmax
      !    write(11,*)dble(vlm2(lm,1)),imag(vlm2(lm,1))
      !  enddo
      !  close(11)
      ! stop
        
      vlm(:,ias)=vlm2(:,1)
      deallocate(vmu2)
    
    enddo !ia
  enddo !is
deallocate(vmu)
end subroutine




end module
