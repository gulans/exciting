subroutine LS_iteration(Ngrid,is,ia,hybx_coef,r, vloc,l,shell_l,shell_occ,sp, nmax, shell0,&
                Nshell,Nspin,relativity,lmax,psi_in,&
                psi,vx_psi,eig )
        ! Ngrid
        ! r
        ! vfull - potential (v_n+v_h+v_xc)
        ! l - quantum number l
        ! num - Number Of eigvals and eigfun to solve
        ! nummax - size of eigval array
        ! eigval (OUT) - erigval array
        ! eigfun (OUT) - eigfun array
use modinteg
Use mod_hybrids, only: ex_coef

use modmain, only: killflag

!--input and output variables--
implicit none

integer, intent(in) :: Ngrid
integer, intent(in) :: is,ia
real(8), intent(in) :: r(Ngrid)
real(8), intent(in) :: vloc(Ngrid),hybx_coef
integer, intent(in) :: l
integer, intent(in) :: shell_l(Nshell)
real(8), intent(in) :: shell_occ(Nshell)
integer, intent(in) :: sp
integer, intent(in) :: nmax
integer, intent(in) :: shell0
integer, intent(in) :: Nshell
integer, intent(in) :: Nspin
logical, intent(in) :: relativity
integer, intent(in) :: lmax
real(8), intent(in) :: psi_in(Ngrid,Nshell,Nspin)

real(8), intent(out) :: psi(Ngrid,nmax)
real(8), intent(out) :: vx_psi(Ngrid,nmax)

real(8), intent(inout) :: eig(nmax)


!integer, intent(in) :: Nrsfun,lmax
!complex(8), intent(in) :: rsfunC(Nrsfun+1,2)
!real(8), intent(in) :: hybx_w(3,2)
!integer, intent(in) :: Ngrid, is,nmax, shell0, Nshell,Nspin
!integer, intent(in) :: l,sp, shell_l(Nshell) 
!real(8), intent(in) :: r(Ngrid),vloc(Ngrid),shell_occ(Nshell,Nspin)
!real(8), intent(in) :: psi_in(Ngrid,Nshell,Nspin),psi_in_p(Ngrid,Nshell,Nspin),F_mix

real(8) :: vx_psi_sr(Ngrid,nmax)
!logical, intent(in) :: relativity
real(8) :: v_rel(Ngrid)
!complex(8)::Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
!real(8), intent(inout) :: eig(nmax)
integer :: iner_loop
!real(8), intent(out) :: psi(Ngrid,nmax)

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)
integer :: inn,inp,ish,i,j 
real(8) :: vx_chi(Ngrid,nmax),vx_chi_sr(Ngrid,nmax)
real(8) :: vx_chi2(Ngrid,nmax),vx_chi2_sr(Ngrid,nmax)
real(8) :: f1(Ngrid),f2(Ngrid),f3(Ngrid),f4(Ngrid),f5(Ngrid)
real(8) :: f6(Ngrid),f7(Ngrid)
real(8) :: phi(Ngrid,nmax),norm,eigp(nmax)
integer :: iscl,maxscl,ir,isp
real(8) :: f(Ngrid)
logical :: spin

real(8) :: t1, vx_psi_in(Ngrid,nmax),vx_psi_sr_in(Ngrid,nmax)
integer :: n1





vx_psi_in=vx_psi
vx_psi_sr_in=vx_psi_sr

vx_chi=psi*0d0
vx_chi_sr=psi*0d0

if (Nspin.eq.2)then
        spin=.true.
else
        spin=.false.
endif




maxscl=40
iner_loop=maxscl+1
do iscl=1,maxscl


!!!!!!!!!!!!!!!!!!!!!!!!
!! Convergence check  !!
!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,*)l,eig,iscl!,". eig-eigp: ",eig-eigp
!convergence check
if((maxval(abs((eig-eigp)/(eig-1d0)))).lt.1d-14)then
        iner_loop=iscl-1
!       write(*,*)"Convergence of internal cycle reached, iteration ",iscl
       !write(*,*)"max(eig-eigp) absolute : ",maxval(abs(eig-eigp))," relative: ",maxval(abs(eig-eigp)/(abs(eig-1d0)))
        exit
endif



if (abs(hybx_coef).gt.1d-20) then
  do inn=1,nmax
    ish=inn+shell0
    call get_Fock_ex(Ngrid,r,is,ia,ish,Nshell,shell_l,&
            shell_occ,lmax,psi(:,inn), psi_in(:,:,sp),vx_psi(:,inn))
  enddo
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Construct right side of scr Poisson  and solve !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do inn=1,nmax
  ish=shell0+inn

  if (.not.relativity) then 

          f=-2d0* (vloc(:)*psi(:,inn)+ hybx_coef *vx_psi(:,inn) )
                
     !     f=-2d0*( vloc(:)*psi(:,inn) + hybx_w(2,1)*vx_psi(:,inn)&
     !            + hybx_w(3,1)*vx_psi_sr(:,inn) )
  else !zora

    f1=1d0/(1d0-alpha2*v_rel)
    call deriv_f(Ngrid,is,log(f1),f2,atom_integw)
    call deriv_f(Ngrid,is,psi_in(:,ish,sp),f3,atom_integw)
    f4=-f2*f3
    
    f2=2d0*vloc(:)*psi_in(:,ish,sp)/f1
    f3=2d0*alpha2*v_rel*eig(inn)*psi_in(:,ish,sp)
    f=f4+f2+f3
    f=-f

  endif
  call scrPoisson(Ngrid,is, r,l, f, eig(inn), psi(:,inn))

  call integ_v(Ngrid,is,r**2*psi(:,inn)**2,norm,atom_integw)

  psi(:,inn)=psi(:,inn)/dsqrt(norm)

enddo



if (abs(hybx_coef).gt.1d-20) then
   do inn=1,nmax
   ish=inn+shell0
   call get_Fock_ex(Ngrid,r,is,ia,ish,Nshell,shell_l,shell_occ,lmax,&
           psi(:,inn),psi_in(:,:,sp),vx_chi(:,inn))
   enddo
endif

eigp=eig

call orthonorm_get_eig(Ngrid,is,r,vloc,l,nmax,relativity,v_rel,hybx_coef,&
        vx_chi,vx_chi_sr,&
        psi,eig,&
        vx_psi,vx_psi_sr)


if (abs(hybx_coef).gt.1d-20) then
  do inn=1,nmax
    ish=inn+shell0
    call get_Fock_ex(Ngrid,r,is,ia,ish,Nshell,shell_l, shell_occ,lmax,&
            psi(:,inn), psi_in(:,:,sp),vx_psi(:,inn))
  
  enddo
endif



enddo !self consistent loop
end subroutine

subroutine scrPoisson(Ngrid,is,r,l,f, e, psi)
use modinteg
integer, intent(in) :: Ngrid
integer, intent(in) :: l,is
real(8), intent(in) :: r(Ngrid)
real(8), intent(inout) :: e
real(8), intent(out) :: psi(Ngrid)

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8) :: f(Ngrid),lam,temp1(Ngrid),temp2(Ngrid)
real(8) :: besrezi(0:50), besrezk(0:50)
integer :: ri,i,endpoint
real(8) :: f11(Ngrid),f12(Ngrid),f21(Ngrid),f22(Ngrid),int1(Ngrid),int2(Ngrid)
!result=f11*(integral_(0->r)f12)+f21(integral_(r->Ngrid)f22)
real(8) :: besi,besk


if (e.gt.0) then
!        write(*,*)"scrPoisson Error: positive eigenvalue!"
        e=-2d-2
endif

lam=dsqrt(-2d0*e)


f11=0d0*r
f12=0d0*r
f21=0d0*r
f22=0d0*r
!write(*,*)"argument       besi         besk"
endpoint=Ngrid
do ri=1, Ngrid
if((lam*r(ri)).gt.200d0) then
!        write(*,*)"endpoint",ri
        endpoint=ri
exit
endif
call msbesseli(l,lam*r(ri), besrezi)
call msbesselk(l,lam*r(ri), besrezk)


besi=besrezi(l)
besk=besrezk(l)
!write(*,*)lam*r(ri),besi,besk

f11(ri)=lam*besk
f12(ri)=besi*f(ri)
f21(ri)=lam*besi
f22(ri)=besk*f(ri)
enddo

call integ_f(Ngrid,is,f12*r**2,int1,atom_integw)


!call integ_f_rev(Ngrid,r,is,f22*r**2,int2,atom_integw)
call integ_f_rev(endpoint,r(1:endpoint),is,f22(1:endpoint)*r(1:endpoint)**2,int2(1:endpoint),atom_integw)
int2(endpoint+1:Ngrid)=0d0*r(endpoint+1:Ngrid)


psi=f11*int1+f21*int2


end subroutine


subroutine diasym_small(a,eig,n)
 implicit none
 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))
 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)
end subroutine

subroutine inver(A,Ainv,si)
    implicit none
    integer            :: si
    real(8)            :: A(si,si)
    real(8)            :: Ainv(si,si)
    real(8)            :: work(si*si)            ! work array for LAPACK
    integer         :: n,info,ipiv(si)     ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = si
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) then
      print *, 'Matrix is numerically singular!'
      stop
    endif
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call DGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) then
      print *,  'Matrix inversion failed!'
      stop
    endif
end subroutine
