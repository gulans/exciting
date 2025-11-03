subroutine getu(m,is,ia,nodes,nrmt,v,l,e,u0,u,u1,q0,q1)
use modinteg
Use mod_hybrids, only: ex_coef
use modmain, only: spnr,spr
use modinput, only: input        
        implicit none
integer, intent(in) :: nrmt,l,is,ia,m
real(8), intent(in) :: v(nrmt),e,u0(nrmt)
real(8), intent(out) :: u(nrmt)
real(8), intent(out) :: u1(nrmt)
integer, intent(out) :: nodes
real(8), intent(out) :: q0(nrmt)
real(8), intent(out) :: q1(nrmt)
real(8) :: p0p(nrmt)


real(8) :: rmt(nrmt)
real(8) :: uold(nrmt),vx_u(nrmt)
real(8) :: vx_uatom(spnr(is)),uatom(spnr(is))
integer :: iter
integer :: ii,ir
real(8) :: diff,norm
real(8) :: rmfactor
logical :: original

integer :: nratom
real(8) :: ratom(spnr(is))
integer :: maxit
nratom=spnr(is)
ratom(:)=spr(:spnr(is),is)

original=.true.

if (input%groundstate%ValenceRelativity.eq."zora") then
  rmfactor=1d0
else
  rmfactor=0d0
endif

rmt=ratom(1:nrmt)

vx_u=0d0
!write(*,*)"l=",l,"e=",e

if(original)then


  if (m.eq.0) then
    p0p=0d0
    Call rschrodint (m, l, e, nrmt, rmt, v, nodes, rmfactor, p0p, u, u1, q0, q1)
  else
    p0p=-(vx_u/dble(m)-u0)
    Call rschrodint (m, l, e, nrmt, rmt, v, nodes, rmfactor, p0p, u, u1, q0, q1)
  endif


else
write(*,*)"getu.f90 "
stop
      
        if (m.eq.0) then
!    call getu_iter(nrmt,rmt,v,l,nodes,e,u,vx_u)
  else
!    call getu_iter(nrmt,rmt,v,l,nodes,e,u,vx_u-dble(m)*u0)
  endif


endif


if( ex_coef.ne.0d0) then
maxit=2000
do iter=1, maxit
  uold=u
  uatom=0d0
  uatom(1:nrmt)=u
  call getrFock(nratom,ratom,is,ia,l,uatom,vx_uatom)
  vx_u=ex_coef*vx_uatom(1:nrmt)


if(original)then
  if (m.eq.0) then
    p0p=-vx_u
    Call rschrodint (1, l, e, nrmt, rmt, v, nodes, rmfactor, p0p, u, u1, q0, q1)
  else
    p0p=-(vx_u/dble(m)-u0)
    Call rschrodint (m, l, e, nrmt, rmt, v, nodes, rmfactor, p0p, u, u1, q0, q1)
  endif
else
  if (m.eq.0) then
!    call getu_iter(nrmt,rmt,v,l,nodes,e,u,vx_u)
  else
!    call getu_iter(nrmt,rmt,v,l,nodes,e,u,vx_u-dble(m)*u0)
  endif
endif

  diff = sum((u-uold)**2)/sum(u**2)

  if (abs(diff).lt.1e-15) exit
  if (iter.ge.maxit)then
    write(*,'("WARNING!! getu l=",i1,"m=",i1,"e=",F6.2," it:",i4, "sum:", ES10.1E3," diff:", ES10.1E3)')l,m,e,iter,sum(u**2),diff
  endif 
enddo

if (iter.gt.40)then
  write(*,'("getuOK l=",i1,"m=",i1,"e=",F6.2," it:",i4, "sum:", ES10.1E3," diff:", ES10.1E3,"n:",i1)')l,m,e,iter,sum(u**2),diff,nodes
endif


endif !ex_coef.ne.0d0



end subroutine

