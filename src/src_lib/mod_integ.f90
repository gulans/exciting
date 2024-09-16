!> Module to evaluate radial integrals and derivatives using Newton-Cotes and Lagrange polynomial interpolation methods. 
module modinteg
implicit none 
  !> Data type that holds arrays of coefficients that are used for integral and derivative evaluation. 
Type IntegWeightsType
      real(8), allocatable :: fintw(:, :, :)   !!array of weights to evaluate integral as a function
      real(8), allocatable :: intw(:, :)       !!array of weights to evaluate integral value
      real(8), allocatable :: fderivw(:, :, :) !!array of weights to evaluate derivative as a function
      integer,allocatable :: istart(:,:)       !!array of grid point indexes for the beginning of the range that is being used for
                                               !!interpolation during integral calculation 
      integer,allocatable :: dstart(:,:)       !!array of grid point indexes for the beginning of the range that is being used for
                                               !!interpolation during derivative calculation

End Type IntegWeightsType

Type(IntegWeightsType) :: mt_integw    !!all the arrays of coefficients that are used for integral and derivative evaluation
                                       !!for mt region for all the atom species
Type(IntegWeightsType) :: atom_integw  !!all the arrays of coeficients that are used for integral and derivative evaluation
                                       !!for atom grid for all the atom species
integer :: i_order           !! order of the interpolated polynomial (can be set from 1 to 10) for integral calculation.
integer :: d_order           !! order of the interpolated polynomial (can be set from 1 to 10) for derivative calculation.
logical :: integrate_0_r1    !! region between r=0 and r(1) is included in the integrals, if integrate_0_r1=True
integer :: i_order2          !! order of the polynomial used for interpolation between r=0 and r(1) (can be set from 1 to 3)
real(8),allocatable :: w2(:) !! array of Newton-Cotes coefficients for interpolation between r=0 and r(1). 

contains

!> Stores all the needed arrays in IntegWeigthType data type variables for each species for MT grid and atom grid.
subroutine gen_icoef(nspecies,nrmax,nrmt,nratom,r)
implicit none
integer, intent(in) :: nspecies          !! number of atom species
integer, intent(in) :: nrmax             !! the largest grid size   
integer, intent(in) :: nrmt(nspecies)    !! MT grid size for all species 
integer, intent(in) :: nratom(nspecies)  !! atom grid size for all species 
real(8), intent(in) :: r(nrmax,nspecies) !! all the radial grids for all species

integer :: iNpoints, dNpoints, Npoints2
integer :: isp
integer :: ir,i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Variables to adjust module settings !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Optimal parameters

d_order=10  
i_order=8
integrate_0_r1=.true.
i_order2=1 ! 1,2 or 3   

!!!!!!!!!!!!!!!!!!!! 
! End of Variables !
!!!!!!!!!!!!!!!!!!!!

Npoints2=i_order2+1
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
!write(*,*)"species",isp,"mt",nrmt(isp),"atom",nratom(isp)
!write(*,*)"r(mt)",r(nrmt(isp),isp),"r(inf)",r(nratom(isp),isp)

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


  !> Generates and returns all the needed coefficients for one given grid.
subroutine gen_icoef_one_grid(Ngrid,r,iNpoints,dNpoints,istart1,dstart1,icv,icf,dcf)
implicit none
integer, intent(in) :: Ngrid             !! number of grid points
real(8), intent(in) :: r(Ngrid)          !! grid 
integer, intent(in) :: iNpoints          !! number of grid points being used in interpolation during integral calculations
integer, intent(in) :: dNpoints          !! number of grid points being used in interpolation during derivative calculations
integer, intent(out) :: istart1(Ngrid)   !! array of grid point indexes for the beginning of the range that is being used for
                                         !! interpolation during integral calculation
integer, intent(out) :: dstart1(Ngrid)   !! array of grid point indexes for the beginning of the range that is being used for
                                         !! interpolation during derivative calculation
real(8), intent(out) :: icv(Ngrid)          !! array of weights to evaluate integral value
real(8), intent(inout) :: dcf(Ngrid,dNpoints) !! array of weights to evaluate derivative as a function
real(8), intent(inout) :: icf(Ngrid,iNpoints) !! array of weights to evaluate integral as a function


real(8) :: w(iNpoints)
real(8) :: lc(Ngrid,iNpoints,iNpoints-2)
real(8) :: xx1(iNpoints)   !Grid section
real(8) :: xx2(dNpoints)   !Grid section
integer :: Npoints2
integer :: ir,i,j
real(8) :: hh(Ngrid)

Npoints2=i_order2+1

!--------------------------!
! Newton-Cotes coefficents !
!--------------------------!

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

!-----------------------!
! Generate istart array !
!-----------------------!
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

!-----------------------!
! Generate dstart array !
!-----------------------!

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

!------------------------------------!
! Generate interpolation coeffiencts !
!------------------------------------!

xx1(1)=0d0
xx1(2:Npoints2)=r(1:Npoints2-1)
do j=1,Npoints2-2
   call lagr_c(0d0+j*r(1)/(Npoints2-1),Npoints2,xx1,lc(1,1:Npoints2,j))
enddo

do ir=1,Ngrid-1
   do j=1,iNpoints-2  !How many points have to be interpolated in between grid points
     xx1=r(istart1(ir):istart1(ir)+iNpoints-1)
     call lagr_c(r(ir)+j*hh(ir),iNpoints,xx1,lc(ir+1,:,j))
   enddo
enddo

!--------------------------------------------------------------------!
! Build an array of coficients for integral evaluation as a function !
!--------------------------------------------------------------------!

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

!------------------------------------------------------------!
! Build an array of coficients for integral value evaluation !
!------------------------------------------------------------!

icv=r*0d0
        
do i=2, Npoints2
    icv(i-1)=icv(i-1)+icf(1,i)
enddo

do ir=1, Ngrid-1
  do i=1, iNpoints
     icv(istart1(ir)+i-1)=icv(istart1(ir)+i-1)+icf(ir+1,i)
  enddo

enddo

!--------------------------------------------------------!
! Build an array of coficients for derivative evaluation !
!--------------------------------------------------------!

do ir=1,Ngrid
  xx2=r(dstart1(ir):dstart1(ir)+dNpoints-1)
  call lagr_c_d(r(ir),dNpoints,xx2,dcf(ir,:))
enddo

end subroutine

!> Evaluates derivative for the given function
subroutine deriv_f(Ngrid,isp,fin,fout,integw)
integer, intent(in) :: Ngrid !! Number of grid points
integer, intent(in) :: isp   !! species index
real(8), intent(in) :: fin(Ngrid) !! input function
real(8), intent(out) :: fout(Ngrid) !! output function (derivative)
Type(IntegWeightsType) , intent(in):: integw !! coefficient array, that lets the subroutine distinguish if this is a MT or atom grid 

real(8) :: r1
integer :: i,ir
integer :: Npoints

Npoints=d_order+1
do ir=1, Ngrid
  r1=0d0
  do i=1, Npoints
    r1=r1+integw%fderivw(ir,i,isp)*fin(integw%dstart(ir,isp)+i-1)
  enddo
  fout(ir)=r1
enddo
end subroutine

!> Evaluates integral function for the given function
subroutine integ_f(Ngrid,isp,fin,fout,integw)
integer, intent(in) :: Ngrid  !! Number of grid points
integer, intent(in) :: isp    !! species index
real(8), intent(in) :: fin(Ngrid) !! input function
real(8), intent(out) :: fout(Ngrid) !! output function (derivative)
Type(IntegWeightsType) , intent(in):: integw !! coefficient array, that lets the subroutine distinguish if this is a MT or atom grid

real(8) :: r1
integer :: ir,i,Npoints,Npoints2

Npoints2=i_order2+1

Npoints=i_order+1

if (integrate_0_r1) then
  r1=0d0
  r1=r1+integw%fintw(1,1,isp)*0d0
    do i=2, Npoints2
      r1=r1+integw%fintw(1,i,isp)*fin(i-1)
    enddo
  fout(1)=r1
else
  fout(1)=0d0
endif

do ir=1, Ngrid-1
  r1=0d0
  do i=1, Npoints
    r1=r1+integw%fintw(ir+1,i,isp)*fin(integw%istart(ir,isp)+i-1)
  enddo
  fout(ir+1)=fout(ir)+r1
enddo
end subroutine

subroutine integ_cf(Ngrid,isp,fin,rez,integw)
implicit none
integer, intent(in) :: Ngrid
integer, intent(in) :: isp
complex(8), intent(in)  :: fin(Ngrid)
complex(8), intent(out) :: rez(Ngrid)
Type(IntegWeightsType) , intent(in):: integw 
integer :: ir

real(8) :: finRe(Ngrid),finIm(Ngrid),rezRe(Ngrid),rezIm(Ngrid)

do ir=1, Ngrid
  finRe(ir)=dble(fin(ir))
  finIm(ir)=imag(fin(ir))
enddo

call integ_f(Ngrid,isp,finRe,rezRe,integw)
call integ_f(Ngrid,isp,finIm,rezIm,integw)

do ir=1, Ngrid
  rez(ir)=cmplx(rezRe(ir),rezIm(ir),8)
enddo
end subroutine

subroutine integ_cv(Ngrid,isp,fin,rez,integw)
  implicit none
  integer, intent(in) :: Ngrid
  integer, intent(in) :: isp
  complex(8), intent(in)  :: fin(Ngrid)
  complex(8), intent(out) :: rez
  Type(IntegWeightsType) , intent(in):: integw 
  integer :: ir
  
  real(8) :: finRe(Ngrid),finIm(Ngrid),rezRe,rezIm
  
!  do ir=1, Ngrid
!    finRe(ir)=dble(fin(ir))
!    finIm(ir)=imag(fin(ir))
!  enddo
  
!  call integ_v(Ngrid,isp,finRe,rezRe,integw)
!  call integ_v(Ngrid,isp,finIm,rezIm,integw)
  
!  rez=cmplx(rezRe,rezIm,8)

 rez=0d0

!can be done with dot product
 do ir=1, Ngrid
   rez=rez+integw%intw(ir,isp)*fin(ir)
 enddo
 
  end subroutine

subroutine integ_f_rev(Ngrid,r,isp,fin,fout,integw)
integer, intent(in) :: Ngrid  !! Number of grid points
integer, intent(in) :: isp    !! species index
real(8), intent(in) :: r(Ngrid) !! 
real(8), intent(in) :: fin(Ngrid) !! input function
real(8), intent(out) :: fout(Ngrid) !! output function (derivative)
Type(IntegWeightsType) , intent(in):: integw !! coefficient array, that lets the subroutine distinguish if this is a MT or atom grid

real(8) :: r1
integer :: ir,i,Npoints,Npoints2

Npoints2=i_order2+1

Npoints=i_order+1


fout(Ngrid)=0d0

do ir=Ngrid-1, Ngrid-int(Npoints/2),-1
  !interpolation gives incorect results if the fin functioon is small in the current point, but is a few orders higher around the
  !point ,trapezoidal rule is used insead:
  fout(ir)=fout(ir+1)+(r(ir+1)-r(ir))*(fin(ir)+fin(ir+1))/2d0
!write(*,*),ir,fout(ir),fin(ir)
enddo




do ir=Ngrid-int(Npoints/2)-1, 1,-1

  r1=0d0
  do i=1, Npoints
    r1=r1+integw%fintw(ir+1,i,isp)*fin(integw%istart(ir,isp)+i-1)
  enddo
  fout(ir)=fout(ir+1)+r1

!write(*,*),ir,fout(ir),fin(ir)

enddo
end subroutine




  !> Evaluates integral value for the given function
subroutine integ_v(Ngrid,isp,fin,vout,integw)
integer, intent(in) :: Ngrid  !! Number of grid points
integer, intent(in) :: isp    !! species index
real(8), intent(in) :: fin(Ngrid) !! input function
real(8), intent(out) :: vout !! output (result)
Type(IntegWeightsType) , intent(in):: integw  !! coefficient array, that lets the subroutine distinguish if this is a MT or atom grid

real(8) :: r1
integer :: ir
r1=0d0

!can be done with dot product
do ir=1, Ngrid
   r1=r1+integw%intw(ir,isp)*fin(ir)
enddo
vout=r1
end subroutine


   !> Calculates Lagrange interpolation coefficients
subroutine lagr_c(xx,k,x,l)
implicit none

real(8), intent(in) :: xx   !! coordinate on the grid where the interpolation is needed
integer, intent(in) :: k    !! k-1 - order of the polynomial
real(8), intent(in) :: x(k) !! section of the grid region that is being used for interpolation
real(8), intent(out):: l(k) !! result - Lagrange interpolation coefficients

real(8) :: prod
integer :: i,j,m

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

   !> Calculates Lagrange interpolation coefficients for derivative
subroutine lagr_c_d(xx,k,x,l)
implicit none
real(8), intent(in) :: xx   !! coordinate on the grid where the derivative is needed
integer, intent(in) :: k    !! k-1 - order of the polynomial
real(8), intent(in) :: x(k) !! section of the grid region that is being used for interpolation
real(8), intent(out):: l(k) !! result - Lagrange interpolation coefficients for derivative

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

!>Deallocates both wariables mt_integw and atom_integw
subroutine dealloc_icoef()
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
deallocate(w2)
end subroutine

end module
