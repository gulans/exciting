subroutine LS_iteration(Ngrid, r,tools,tools_info,l,sp, nmax,&
                Nspin,vr,psi,eig)
        ! Ngrid (in)
        ! r (in)
        ! vr - full potential (in)
        ! l - quantum number l (in)
        ! sp - spin up or down (in)
	! nmax - number of orbitals for a specific l (in)
	! Nspin - if 1 spin not taken into account, if 2 then is ( for now 2 is not implemented) (in)
        ! psi - wavefunction (inout)
        ! eig - eigenvalues (inout)

implicit none

integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: Ngrid, nmax, Nspin
integer, intent(in) :: l,sp
real(8), intent(in) :: r(Ngrid)
real(8), intent(in) :: vr(Ngrid)
real(8), intent(inout) :: eig(nmax)
real(8), intent(inout) :: psi(Ngrid,nmax,Nspin)

integer :: inn,iscl,maxscl,ir
real(8) :: norm,eigp(nmax)
real(8) :: f(Ngrid), psi_in(Ngrid,nmax,1), eig_in(nmax,1)

real(8) :: elimit, eshift,wf_sum
logical :: new_algorithm
logical :: eig_limiter


new_algorithm=.true.
eig_limiter=.false.
elimit=-0.0002d0

eigp=0d0
maxscl=200 !maximum number of iterations



!________________________________
!check whether it is the first iteration in scf calculation
wf_sum = sum(psi(:,1,1))
if ((wf_sum.lt.1).or.(wf_sum.gt.200000)) then
	call iteration0(Ngrid,r,1d0,1,1,l,nmax,1,eig_in,psi_in)!guess for wavefunction, filled with laguerra polynomials
	write(*,*)"0th iteration"
	psi(:,:,sp) =psi_in(:,:,1)
else
	write(*,*)"non 0th iteration"
end if

do iscl=1,maxscl!self consistent loop

!---------------convergence check
	if((maxval(abs((eig-eigp)/(eig-1d0)))).lt.1d-15)then
		write(*,*)"Convergence of internal cycle reached, iteration ",iscl
       		write(*,*)"max(eig-eigp) absolute : ",maxval(abs(eig-eigp))," relative: ",maxval(abs(eig-eigp)/(abs(eig-1d0)))
		exit
	endif


!--------solve screened Poisson equation for different n, equal l orbitals
	do inn=1,nmax
  		f=-2d0*(vr(:))*psi(:,inn,sp) 
		elimit=-2d0
		if ((eig_limiter).and.(eig(inn).gt.elimit))then
		 	eshift=-eig(inn)+elimit
			eig(inn)=elimit
			
		  	if(new_algorithm)then
				call scrPoisson(Ngrid, r,tools,tools_info,l, f-2*(eshift)*psi(:,inn,sp), eig(inn), psi(:,inn,sp))
			else
				call scrPoisson(Ngrid, r,tools,tools_info,l, f-2*(eshift)*psi(:,inn,sp), eig(inn), psi(:,inn,sp))
		  	endif    
		  	eig(inn)=eig(inn)-eshift
		else
		  	call scrPoisson(Ngrid, r,tools,tools_info,l, f, eig(inn), psi(:,inn,sp))
		endif
  		call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn,sp)**2,norm)
  		psi(:,inn,sp)=psi(:,inn,sp)/dsqrt(norm)
	enddo


eigp=eig

!------------get eigenvalues and wavefunctions drom diagonalization
call diag_get_eig(Ngrid,r,tools,tools_info,l,vr,nmax,psi,eig)
enddo !self consistent loop



end subroutine

subroutine scrPoisson(Ngrid, r,tools,tools_info,l,f, e, psi)

integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: Ngrid
integer, intent(in) :: l
real(8), intent(in) :: r(Ngrid)
real(8), intent(inout) :: e
real(8), intent(out) :: psi(Ngrid)

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8) :: f(Ngrid),lam
real(8) :: besrezi(0:50), besrezk(0:50)
integer :: ri,i,endpoint
real(8) :: f11(Ngrid),f12(Ngrid),f21(Ngrid),f22(Ngrid),int1(Ngrid),int2(Ngrid)
!result=f11*(integral_(0->r)f12)+f21(integral_(r->Ngrid)f22)
real(8) :: besi,besk
real(8), allocatable :: r_e(:), f22_e(:), int2_e(:), tools_e(:,:)

!write(*,*)"e=",e
if (e.gt.0) then
!        write(*,*)"scrPoisson Error: positive eigenvalue!"
        e=-2d-2
endif

lam=dsqrt(-2d0*e)


f11=0d0*r
f12=0d0*r
f21=0d0*r
f22=0d0*r

endpoint=Ngrid
do ri=1, Ngrid
if((lam*r(ri)).gt.200d0) then

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
allocate(r_e(endpoint), f22_e(endpoint), int2_e(endpoint), tools_e(endpoint, tools_info(1)))
call integ_BodesN_fun(Ngrid,r,tools,tools_info,1,f12*r**2,int1)
tools_e = tools
r_e = r(1:endpoint)
tools_e = tools(1:endpoint,:)
f22_e = f22(1:endpoint)

call integ_BodesN_fun(endpoint, r_e, tools_e, tools_info,-1,f22_e*r_e**2,int2_e)
int2(endpoint+1:Ngrid)=0d0*r(endpoint+1:Ngrid)
int2(1:endpoint) = int2_e

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
