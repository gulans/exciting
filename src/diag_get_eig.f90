subroutine diag_get_eig(Ngrid,r,tools,tools_info,l,vr,nmax,&
	psi,eig)
!------------Solves general eigenvalue problem Hx = \lambda Sx
implicit none
integer(8), intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: l,nmax
real(8), intent(in) :: vr(Ngrid)
real(8), intent(inout) :: psi(Ngrid,nmax) 
real(8), intent(out) :: eig(nmax)


integer :: inn,inp,i,j
real(8) :: Snn,Hnn
real(8) :: S(nmax,nmax),H(nmax,nmax),lambda(nmax)
real(8) :: Hevec(nmax,nmax),Hp(nmax,nmax), Sevec(nmax,nmax),Seval(nmax), s12(nmax,nmax)
real(8) :: x(nmax,nmax),W(nmax,nmax),Winv(nmax,nmax)
real(8) :: psi_new(Ngrid,nmax)
real(8) :: f(Ngrid),f1(Ngrid),f2(Ngrid)


!-----calculate overlap matrix S_{ij} = <\psi_i|psi_j>
!-----calculate Hamiltonian matrix H
do inn=1,nmax !matrix is symetric !lower triangular matrix 
do inp=inn,nmax
    call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn)*psi(:,inp),Snn)
    S(inn,inp)=Snn
    S(inp,inn)=Snn

    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inp),f2)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inn),f1)
    call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*f2*f1*r**2,H(inn,inp))


    f=(0.5d0*dble(l)*dble(l+1)/r**2+vr(:))*psi(:,inp)
           
    call integ_BodesN_value(Ngrid,r,tools,tools_info,psi(:,inn)*f*r**2,Hnn)

    H(inn,inp)=H(inn,inp)+Hnn
    H(inp,inn)=H(inn,inp)

    enddo
enddo


!!!!!!!!Caulculate W=S^-0.5 matrix
 Sevec=S
 call diasym_small(Sevec,Seval,nmax)

!!!!!!!! construct s12 matrix, (eigenvalue diagonal-matrix with elements in power -0.5)
 do i=1, nmax
   do j=1, nmax
     if (i==j) Then
       s12(i,j)=Seval(i)**(-0.5d0)
     else
       s12(i,j)=0d0
     endif
   enddo
 enddo
 W=matmul(matmul(Sevec,s12),transpose(Sevec))

!!!!!!!!! Make a substitution x=Wx' un H'=WHW
call inver(W,Winv,nmax)


 Hp=matmul(matmul(W,H),W)
 Hevec=Hp
!!!!!!!! SOLVE H'x'=x'
 call diasym_small(Hevec,lambda,nmax)




 x=matmul(W,Hevec)


!!CONSTRUCT wave functions
do inn=1,nmax
  psi_new(:,inn)=0d0*r
  do inp=1,nmax
    psi_new(:,inn)=psi_new(:,inn)+x(inp,inn)*psi(:,inp)
  enddo 
enddo

!-------output wavefunctions and eigenvalues
do inn=1,nmax
  psi(:,inn)=psi_new(:,inn)
  eig(inn)=lambda(inn)
enddo


end subroutine




