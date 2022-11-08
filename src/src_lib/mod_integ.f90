!> Module to evaluate integrals and derivatives using Newton-Cotes and Lagrange polynomial interpolation methods respectively. 
module modinteg
implicit none 
  !> Data type that holds arrays of coeficients that are used for integral and derivative evaluation. 
Type IntegWeightsType
      real(8), allocatable :: fintw(:, :, :) !!array of weigts to evaluate integral as a funciton
      real(8), allocatable :: intw(:, :)     !!array of weigts to evaluate integral value
      real(8), allocatable :: fderivw(:, :, :)
      integer,allocatable :: istart(:,:)
      integer,allocatable :: dstart(:,:)
End Type IntegWeightsType

Type(IntegWeightsType) :: mt_integw, atom_integw


logical :: integrate_0_r1
integer :: d_order
integer :: i_order
integer :: Npoints2
real(8),allocatable :: w2(:)

contains

  !> Generates IntegWeigthType data type variables that consists of arrays for each species for grid of MT size and atom size
subroutine gen_icoef(nspecies,nrmax,nrmt,nratom,r)

implicit none
integer, intent(in) :: nspecies
integer, intent(in) :: nrmax
integer, intent(in) :: nrmt(nspecies)
integer, intent(in) :: nratom(nspecies)
real(8), intent(in) :: r(nrmax,nspecies)

integer :: iNpoints, dNpoints
integer :: isp
integer :: ir,i

!!!!!global module variables!!!!!!

d_order=9
i_order=9
integrate_0_r1=.true.
Npoints2=2 ! 2,3 or 4   
allocate(w2(Npoints2))

iNpoints=i_order+1
dNpoints=d_order+1

allocate(mt_integw%fintw    (nrmax,iNpoints,nspecies))
allocate(mt_integw%intw    (nrmax,nspecies))
allocate(mt_integw%fderivw (nrmax,dNpoints,nspecies))
allocate(mt_integw%istart  (nrmax,nspecies))
allocate(mt_integw%dstart  (nrmax,nspecies))

allocate(atom_integw%fintw    (nrmax,iNpoints,nspecies))
allocate(atom_integw%intw    (nrmax,nspecies))
allocate(atom_integw%fderivw (nrmax,dNpoints,nspecies))
allocate(atom_integw%istart  (nrmax,nspecies))
allocate(atom_integw%dstart  (nrmax,nspecies))


mt_integw%fintw=0d0
mt_integw%intw=0d0
mt_integw%fderivw=0d0
mt_integw%istart=0
mt_integw%dstart=0

atom_integw%fintw=0d0   
atom_integw%intw=0d0   
atom_integw%fderivw=0d0
atom_integw%istart=0
atom_integw%dstart=0




if (Npoints2.eq.2) then !Trapezoidal rule
        w2(1)=1d0
        w2(2)=1d0
        w2=w2/2d0
elseif (Npoints2.eq.3) then !Simpson's rule
        w2(1)=1d0
        w2(2)=4d0
        w2(3)=1d0
        w2=1d0*w2/3d0
elseif (Npoints2.eq.4) then !Simpson's 3/8 rule
        w2(1)=1d0
        w2(2)=3d0
        w2(3)=3d0
        w2(4)=1d0
        w2=3d0*w2/8d0
endif


do isp=1, nspecies
  !write(*,*)"species",isp," MT grid points",nrmt(isp)," atom: ",nratom(isp)
  call gen_icoef_one_grid(nrmt(isp),&
          r(1:nrmt(isp),isp),iNpoints,dNpoints,&
          mt_integw%istart(1:nrmt(isp),isp),mt_integw%dstart(1:nrmt(isp),isp),&
          mt_integw%intw(1:nrmt(isp),isp),&
          mt_integw%fintw(1:nrmt(isp),:,isp),&
          mt_integw%fderivw(1:nrmt(isp),:,isp) )

  call gen_icoef_one_grid(nratom(isp),&
          r(1:nratom(isp),isp),iNpoints,dNpoints,&
          atom_integw%istart(1:nratom(isp),isp),atom_integw%dstart(1:nratom(isp),isp),&
          atom_integw%intw(1:nratom(isp),isp),&
          atom_integw%fintw(1:nratom(isp),:,isp),&
          atom_integw%fderivw(1:nratom(isp),:,isp) )

enddo



end subroutine



subroutine gen_icoef_one_grid(Ngrid,r,iNpoints,dNpoints,istart1,dstart1,icv,icf,dcf)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: iNpoints
integer, intent(in) :: dNpoints
integer, intent(out) :: istart1(Ngrid)
integer, intent(out) :: dstart1(Ngrid)
real(8), intent(out) :: icv(Ngrid)          !icv(Integral Coefficients value) 
real(8), intent(out) :: dcf(Ngrid,dNpoints) !dcf(derivative Coefficients function) 
real(8), intent(out) :: icf(Ngrid,iNpoints) !icf(Integral Coefficients function)


real(8) :: w(iNpoints)
real(8) :: lc(Ngrid,iNpoints,iNpoints-2)
real(8) :: xx1(iNpoints)   !Grid section
real(8) :: xx2(dNpoints)
integer :: ir,i,j

real(8) :: hh(Ngrid)

if (iNpoints.eq.2) then !Trapezoidal rule
        w(1)=1d0
        w(2)=1d0
        w=w/2d0
elseif (iNpoints.eq.3) then !Simpson's rule
        w(1)=1d0
        w(2)=4d0
        w(3)=1d0
        w=1d0*w/3d0
elseif (iNpoints.eq.4) then !Simpson's 3/8 rule
        w(1)=1d0
        w(2)=3d0
        w(3)=3d0
        w(4)=1d0
        w=3d0*w/8d0
elseif (iNpoints.eq.5) then ! Boole's rule
        w(1)=7d0
        w(2)=32d0
        w(3)=12d0
        w(4)=32d0
        w(5)=7d0
        w=2d0*w/45d0
elseif (iNpoints.eq.6) then
        w(1)=19d0
        w(2)=75d0
        w(3)=50d0
        w(4)=50d0
        w(5)=75d0
        w(6)=19d0
        w=5d0*w/288d0
elseif (iNpoints.eq.7) then
        w(1)=41d0
        w(2)=216d0
        w(3)=27d0
        w(4)=272d0
        w(5)=27d0
        w(6)=216d0
        w(7)=41d0
        w=w/140d0
elseif (iNpoints.eq.8) then
        w(1)=751d0
        w(2)=3577d0
        w(3)=1323d0
        w(4)=2989d0
        w(5)=w(4)
        w(6)=w(3)
        w(7)=w(2)
        w(8)=w(1)
        w=7d0*w/17280d0
elseif (iNpoints.eq.9) then
        w(1)=989d0
        w(2)=5888d0
        w(3)=-928d0
        w(4)=10496d0
        w(5)=-4540d0
        w(6)=w(4)
        w(7)=w(3)
        w(8)=w(2)
        w(9)=w(1)
        w=4d0*w/14175d0
elseif (iNpoints.eq.10) then
        w(1)=2857d0
        w(2)=15741d0
        w(3)=1080d0
        w(4)=19344d0
        w(5)=5778d0
        w(6)=w(5)
        w(7)=w(4)
        w(8)=w(3)
        w(9)=w(2)
        w(10)=w(1)
        w=9d0*w/89600d0
elseif (iNpoints.eq.11) then
        w(1)=16067d0
        w(2)=106300d0
        w(3)=-48525d0
        w(4)=272400d0
        w(5)=-260550d0
        w(6)=427368d0
        w(7)=w(5)
        w(8)=w(4)
        w(9)=w(3)
        w(10)=w(2)
        w(11)=w(1)
        w=5d0*w/299376d0
endif


do ir=1,Ngrid-1
  hh(ir)=(r(ir+1)-r(ir))/(iNpoints-1)
enddo

!!!!!!!!!!!!!!!Generate istart array !!!!!!!!!!!!!!!

do ir=1,int(iNpoints/2)
   istart1(ir)=1
enddo
do ir=int(iNpoints/2)+1, Ngrid-int(iNpoints/2)-1
 if (MOD(iNpoints,2) .eq. 0) then
   istart1(ir)=ir-int(iNpoints/2)+1
 else
   istart1(ir)=ir-int(iNpoints/2)
 endif
enddo
do ir=Ngrid-int(iNpoints/2),Ngrid-1
   istart1(ir)=Ngrid-iNpoints+1
enddo


!!!!!!!!!!!!!!!END Generate istart array !!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!Generate dstart array !!!!!!!!!!!!!!!

do ir=1,int(dNpoints/2)
   dstart1(ir)=1
enddo
do ir=int(dNpoints/2)+1, Ngrid-int(dNpoints/2)-1
 if (MOD(dNpoints,2) .eq. 0) then
   dstart1(ir)=ir-int(dNpoints/2)+1
 else
   dstart1(ir)=ir-int(dNpoints/2)
 endif
enddo
do ir=Ngrid-int(dNpoints/2),Ngrid
   dstart1(ir)=Ngrid-dNpoints+1
enddo


!!!!!!!!!!!!!!!END Generate dstart array !!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!Interpolation coeficients !!!!!!!!!!!!!!!!!!!!!
if(.false.)then
xx1(1)=0d0
xx1(2:iNpoints)=r(1:iNpoints-1)
do j=1,iNpoints-2
   call lagr_c(0d0+j*r(1)/(iNpoints-1),iNpoints,xx1,lc(1,:,j))
enddo
else
xx1(1)=0d0
xx1(2:Npoints2)=r(1:Npoints2-1)
do j=1,Npoints2-2
   call lagr_c(0d0+j*r(1)/(Npoints2-1),Npoints2,xx1,lc(1,1:Npoints2,j))
enddo

endif



do ir=1,Ngrid-1
   do j=1,iNpoints-2  !How many points have to be interpolated in between grid points
     xx1=r(istart1(ir):istart1(ir)+iNpoints-1)
     call lagr_c(r(ir)+j*hh(ir),iNpoints,xx1,lc(ir+1,:,j))
   enddo
enddo

!!!!!!!!!!!!!!!!!!!! END Interpolation coeficients!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!! Build array of coficients to get function

do i=1, iNpoints
  icf(:,i)=0d0*r
enddo

if (integrate_0_r1) then
  icf(1,1)=icf(1,1)+w2(1)
  icf(1,2)=icf(1,2)+w2(Npoints2)
  do i=1, Npoints2
    do j=1, Npoints2-2
      icf(1,i)=icf(1,i)+w2(j+1)*lc(1,i,j)
    enddo
  enddo
endif

do ir=1, Ngrid-1
  do i=1, iNpoints
    if (istart1(ir)+i-1.eq.ir) then
      icf(ir+1,i)=icf(ir+1,i)+w(1)
    elseif (istart1(ir)+i-1.eq.ir+1) then
      icf(ir+1,i)=icf(ir+1,i)+w(iNpoints)
    endif
    do j=1, iNpoints-2
      icf(ir+1,i)=icf(ir+1,i)+w(j+1)*lc(ir+1,i,j)
    enddo
  enddo
enddo


do i=1, iNpoints
  icf(1,i)=icf(1,i)*r(1)/(Npoints2-1)
  do ir=1,Ngrid-1
     icf(ir+1,i)=icf(ir+1,i)*hh(ir)
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!END Build array of coficients to get function

!!!!!!!!!!!!!!!!!!!! Build array of coficients to get integral value

icv=r*0d0

do i=2, Npoints2
    icv(i-1)=icv(i-1)+icf(1,i)
enddo

do ir=1, Ngrid-1
  do i=1, iNpoints
     icv(istart1(ir)+i-1)=icv(istart1(ir)+i-1)+icf(ir+1,i)
  enddo

enddo


!!!!!!!!!!!!!!!!!!!!END Build array of coficients to get integral value




!!!!!!!!!!!!!!!!!!!! Build array of coficients to get derivative function

do ir=1,Ngrid
  xx2=r(dstart1(ir):dstart1(ir)+dNpoints-1)
  call lagr_c_d(r(ir),dNpoints,xx2,dcf(ir,:))
enddo

!!!!!!!!!!!!!!!!!!!! Build array of coficients to get derivative function

end subroutine


subroutine deriv_f(Ngrid,isp,fin,fout,integw)
integer, intent(in) :: Ngrid
integer, intent(in) :: isp
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: fout(Ngrid)
Type(IntegWeightsType) , intent(in):: integw


real(8) :: r1
integer :: i,ir
integer :: Npoints

Npoints=d_order+1
do ir=1, Ngrid
  r1=0d0
  do i=1, Npoints
    r1=r1+integw%fderivw(ir,i,isp)*fin(integw%dstart(ir,isp)+i-1)

    !r1=r1+fderivw_mt(ir,i,isp)*fin(dstart_mt(ir,isp)+i-1)
  enddo
  fout(ir)=r1
enddo
end subroutine



subroutine integ_f(Ngrid,isp,fin,fout,integw)
integer, intent(in) :: Ngrid
integer, intent(in) :: isp
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: fout(Ngrid)
Type(IntegWeightsType) , intent(in):: integw

real(8) :: r1
integer :: ir,i,Npoints


Npoints=i_order+1
if (integrate_0_r1) then

  r1=0d0
  r1=r1+integw%fintw(1,1,isp)*0d0

!  r1=r1+fintw_mt(1,1,isp)*0d0
    do i=2, Npoints2
      r1=r1+integw%fintw(1,i,isp)*fin(i-1)

      !r1=r1+fintw_mt(1,i,isp)*fin(i-1)
    enddo
  fout(1)=r1
else
  fout(1)=0d0
endif

do ir=1, Ngrid-1
  r1=0d0
  do i=1, Npoints
    r1=r1+integw%fintw(ir+1,i,isp)*fin(integw%istart(ir,isp)+i-1)

!    r1=r1+fintw_mt(ir+1,i,isp)*fin(istart_mt(ir,isp)+i-1)
  enddo
  fout(ir+1)=fout(ir)+r1
enddo
end subroutine


subroutine integ_v(Ngrid,isp,fin,vout,integw)
integer, intent(in) :: Ngrid
integer, intent(in) :: isp
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: vout
Type(IntegWeightsType) , intent(in):: integw

real(8) :: r1
integer :: ir
r1=0d0

do ir=1, Ngrid
   r1=r1+integw%intw(ir,isp)*fin(ir)

!   r1=r1+intw_mt(ir,isp)*fin(ir)

enddo
vout=r1
end subroutine



subroutine lagr_c(xx,k,x,l)
implicit none
integer, intent(in) :: k
real(8), intent(in) :: x(k),xx
real(8), intent(out):: l(k)

real(8) :: prod
integer :: i,j,m

!Get coeficients
do j=1,k
 prod=1d0
  do m=1,k
    if (m.ne.j) then
      prod=prod*(xx-x(m))/(x(j)-x(m))
    endif
  enddo
  l(j)=prod
enddo
end subroutine

subroutine lagr_c_d(xx,k,x,l)
implicit none
integer, intent(in) :: k
real(8), intent(in) :: x(k),xx
real(8), intent(out):: l(k)

real(8) :: prod
integer :: i,j,m

!Get coeficients
do j=1,k
  l(j)=0d0
  do i=1,k
    if (i.ne.j) then
      prod=1d0
      do m=1,k
         if ((m.ne.i).and.(m.ne.j)) then
         prod=prod*(xx-x(m))/(x(j)-x(m))
         endif
      enddo
      l(j)=l(j)+prod/(x(j)-x(i))
    endif
  enddo
enddo

end subroutine

subroutine dealoc_icoef()
implicit none

deallocate(mt_integw%fintw)
deallocate(mt_integw%intw)
deallocate(mt_integw%fderivw)
deallocate(mt_integw%istart)
deallocate(mt_integw%dstart)

deallocate(atom_integw%fintw)
deallocate(atom_integw%intw)
deallocate(atom_integw%fderivw)
deallocate(atom_integw%istart)
deallocate(atom_integw%dstart)

end subroutine

end module
