module modsurf

  use mod_atoms, only: natmtot
implicit none

integer              :: naxismax
integer,allocatable  :: naxis(:) !ias
integer,allocatable  :: axisy(:,:),axisz(:,:) !iax,ias
real(8),allocatable  :: axisx_sph(:,:) !iax,ias
complex(8),allocatable :: ylm_mat(:,:,:) !lm,iax,ias
complex(8),allocatable :: ylm_tmat(:,:,:) !lm,iax,ias
real(8),allocatable   :: raxis(:,:,:)! 3 ,iax,ias
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
  real(8),allocatable ::  raxis_tmp(:,:,:) !3,iax,ias
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
  allocate(raxis_tmp(3,ngrid(2)*ngrid(3),natmtot))
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
                  if ((xmu.ge.0d0).and.(xmu.lt.(ngrid(1)))) then
                    nax = nax+1
                    axisx_sph_tmp(nax,ias)=xmu
                    axisy_tmp(nax,ias)=yi
                    axisz_tmp(nax,ias)=zi
                    rv1=r0 + a1*xmu - ratom
                    raxis_tmp(:,nax,ias)=rv1
                    !write(*,*)"rmt1",sqrt(sum((rv1)**2)) 
                    tp_tmp(1,nax,ias) = atan2(sqrt(rv1(1)**2+rv1(2)**2),rv1(3))
                    tp_tmp(2,nax,ias) = atan2(rv1(2),rv1(1))
 
                  endif

                  xmu= ( -sum(a1*rv) + sqrt(disc) )/ sum(a1**2)
                  if ((xmu.ge.0d0).and.(xmu.lt.(ngrid(1)))) then
                    nax = nax+1
                    axisx_sph_tmp(nax,ias)=xmu
                    axisy_tmp(nax,ias)=yi
                    axisz_tmp(nax,ias)=zi
                    rv1=r0 + a1*xmu - ratom
                    raxis_tmp(:,nax,ias)=rv1
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
Allocate(raxis(3,naxismax,natmtot))
ylm_mat=zzero
axisy=0
axisz=0
axisx_sph=0d0
do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    ratom=atposc (:, ia, is) 
    axisy(1:naxis(ias),ias)=axisy_tmp(1:naxis(ias),ias)
    axisz(1:naxis(ias),ias)=axisz_tmp(1:naxis(ias),ias)
    axisx_sph(1:naxis(ias),ias)=axisx_sph_tmp(1:naxis(ias),ias)
    
    do iax=1, naxis(ias)
      raxis(:,iax,ias)=raxis_tmp(:,iax,ias)+ratom
      Call genylm (lmax, tp_tmp(:, iax,ias), ylm_mat(:,iax,ias))
    enddo
  enddo
enddo

do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia,is)
    call inv_svd(naxis(ias),lmmax,transpose(ylm_mat(:,1:naxis(ias),ias)),ylm_tmat(:,1:naxis(ias),ias))
  enddo!ia
enddo!is



deallocate(axisy_tmp,axisz_tmp,axisx_sph_tmp,tp_tmp,raxis_tmp)
open(11,file='surf.dat',status='replace')
open(12,file='surf2.dat',status='replace')
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      do iax=1, naxis(ias)
        rv(:)=axisx_sph(iax,ias)*a1 + axisy(iax,ias)*a2 + axisz(iax,ias)*a3
        write(11,*)rv(1),rv(2),rv(3)
        write(12,*)raxis(1,iax,ias),raxis(2,iax,ias),raxis(3,iax,ias)
      enddo
    enddo
  enddo
close(11)
close(12)


end subroutine 

subroutine surf_pot(lmax,zvclir,igfft,qvec,vlm)
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
  real(8),intent(in) :: qvec(3)
  complex(8) ,intent(inout) :: vlm(:,:)!lm,ias
  complex(8) :: zvax(ngrid(1)), zvaxft(ngrid(1))
  !complex(8) :: zvmu((lmax+1)**2)
  integer :: ig1,ig2,iax,ir,is,ia,ias,iy,iz,imu, nn
  integer :: lmmax,lm
  complex(8),allocatable :: vmu(:,:),vmu2(:,:)
  real(8) :: xx, irf,ir2,v1(3),a1(3),a2(3),a3(3),t1,rv(3)
  complex(8) ::zt1, ii
  integer :: ifg,ig
 
  complex(8) :: vlm2((lmax+1)**2,1),phase
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

        rv(:)=raxis(:,iax,ias)
        
        !!!!!!!!!!!!!!!!!!!!!!
        !!! This fixes the result for ZnO
        !v1=rv
        !!! This fixes the result for ZnO
        !!!!!!!!!!!!!!!!!!!!!!
        
        phase=exp(cmplx(0,1,8)*sum(qvec*rv))

        zt1=zzero
        Do ig = 1, ngvec
          ifg = igfft (ig)
          t1 = vgc(1,ig)*v1(1) + vgc(2,ig)*v1(2) + vgc(3,ig)*v1(3)
          zt1 = zt1 + zvclir(ifg)*cmplx(Cos(t1), Sin(t1), 8)
        End Do
        vmu(iax,ias)=zt1*phase
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
        rv(:)=raxis(:,iax,ias)
        phase=exp(cmplx(0,1,8)*sum(qvec*rv))
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
        vmu(iax,ias)=zt1*phase
      enddo !iax
      ! open(11,file='vmu1.dat',status='replace')
      ! do iax=1,naxis(ias)
      !   write(11,*)dble(vmu(iax,ias)),imag(vmu(iax,ias))
      ! enddo
      ! close(11)

    !Create a set of v_lm from v_mu
      ! allocate(vmu2(naxis(ias),1))
      ! vmu2(:,1)=vmu(1:naxis(ias),ias)
      ! call zlsp(transpose(ylm_mat(:,1:naxis(ias),ias)), vmu2, vlm2)
      ! vlm(:,ias)=vlm2(:,1)
      ! deallocate(vmu2)

    !Second way to do it:
      vlm(:,ias)=matmul(ylm_tmat(:,1:naxis(ias),ias),vmu(1:naxis(ias),ias))

      ! open(11,file='lm-set.dat',status='replace')
      !  do lm=1,lmmax
      !    write(11,*)dble(vlm2(lm,1)),imag(vlm2(lm,1))
      !  enddo
      !  close(11)
      ! stop
        

    
    enddo !ia
  enddo !is
deallocate(vmu)
end subroutine

subroutine inv_svd(m,n,A_in,Ainv)
  implicit none
  integer  , intent(in) :: m, n
  complex(8) ,intent(in) :: A_in(m,n)
  complex(8) ,intent(out) :: Ainv(n,m)
  complex(8) :: A(m,n),U(m,m),VT(n,n)
  real(8)    :: S(n),Sinv(n,m)
  real(8)    :: rwork(5*m) !can be decreased

  integer :: info ,work_size, lm
  complex(8),allocatable  :: work(:) 
  real(8),allocatable     :: Smat(:,:)
  complex(8),allocatable  :: E(:,:)
  
  A=A_in
  work_size=5*m
  allocate(work(work_size))
  !call zgesvd('A', 'A', nax, lmmax, A, nax, S, U, nax, VT, lmmax, work, work_size, rwork, info)
  call zgesvd('A', 'A', m, n, A, m, S, U, m, VT, n, work, work_size, rwork, info)

  deallocate(work)
  write(*,*)"info",info
  if (.false.) then !check if decomposition was correct
    write(*,*)"matrix before SVD:"
    write(*,*)A_in(1,1:3)
    write(*,*)A_in(2,1:3)
    write(*,*)A_in(3,1:3)
    allocate(Smat(m,n))
    Smat=0d0
    do lm=1,n
      Smat(lm,lm)=S(lm)
    enddo
    A = matmul( matmul(U ,Smat) , VT )
    write(*,*)"matrix constructed back:"
    write(*,*)A(1,1:3)
    write(*,*)A(2,1:3)
    write(*,*)A(3,1:3)
    deallocate (Smat)
  endif
  Sinv=0d0
  do lm=1,n
    if (abs(S(lm)).gt.1e-10) then
      Sinv(lm,lm)=1d0/S(lm)
    else
      write(*,*)"Matricas dalīšana SVD bija aizdomīga"
      Sinv(lm,lm)=0d0
    endif
  enddo


  Ainv = matmul( matmul(conjg(transpose(VT)),Sinv)  , conjg(transpose(U)))

  if (.false.) then !check if A^(-1)*A=E
    Allocate(E(n,n))
    E=matmul(Ainv,A)
    write(*,*)"mazā Vienības matica:"
    write(*,*)E(1,1:3)
    write(*,*)E(2,1:3)
    write(*,*)E(3,1:3)
    Write(*,*)"Dioganal:"
    do lm=1,20!n
      write(*,*)E(lm,lm)
    enddo
    deallocate(E)
  endif

end subroutine


end module
