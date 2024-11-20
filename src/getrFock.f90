subroutine getrFock(Ngrid,r,is,ia,l,u,vx_u)
use modmain, only: mt_dm,apword,nlorb,lorbl,idxlm,idxas,nrmt,lofr,apwfr,spl,spocc
use mod_potential_and_density, only: dm_copy
use modinput!, only: input
use mod_corestate, only: rwfcr,c_list,c_count

use modinteg
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0

integer, intent(in) :: is,ia
integer, intent(in)  :: Ngrid,l
real(8), intent(in)  :: r(Ngrid),u(Ngrid)
real(8), intent(out) :: vx_u(Ngrid)

integer ::shell_l(c_count(is))
real(8) ::u_all(Ngrid,c_count(is)),shell_occ(c_count(is))


!complex(8) ::rez(Ngrid)
!complex(8) :: integc1(Ngrid),integc2(Ngrid)
complex(8) :: t1c,t2c
integer :: k,i,ir,ish,lpri,lpripri
real(8) :: gc
real(8) :: integ(Ngrid)

integer :: l1,io1,io2,m1,m2,if1,if1old,lm1,ias
integer :: if2,if2old,ioo
integer :: ilo1

integer ::   maxnlo, maxaa, wfsize,maxnorb
integer, pointer :: losize(:)
real(8),allocatable :: basis(:,:,:)
integer,allocatable :: matmap(:,:)
integer :: norb(0:input%groundstate%lmaxmat) 
integer :: tempindex(0:input%groundstate%lmaxmat)
real(8) :: vx_u_part(Ngrid),occ



ias=idxas (ia, is)

u_all(:,:)=rwfcr(1:Ngrid,1,c_list(1:c_count(is),is),ias)
shell_l(:)=spl(1:c_count(is),is)
shell_occ(:)=spocc(1:c_count(is),is)




vx_u=0d0


!write(*,*)"Fock sākums"





if ((mt_dm%maxnlo.ne.0) .or. (mt_dm%maxaa.ne.0))then
!if (.false.)then

maxnlo=mt_dm%maxnlo
losize=>mt_dm%losize
maxaa=mt_dm%maxaa
wfsize=maxaa+maxnlo

if (.false.)then !write dm to a file
  write(*,*)"maxaa",maxaa
  write(*,*)"dm size",wfsize
  write(*,*)"lmaxmat",input%groundstate%lmaxmat

  open (11, file = 'dm.dat', status = 'replace')
  do if1=1, wfsize
    do if2=1, wfsize
      write(11, '(E11.4,",")',advance="no")dreal(mt_dm%main%ff(if1,if2,ias))
    enddo
    write(11,*)
  enddo
  write(11,*)

  do if1=1, wfsize
    do if2=1, wfsize
      write(11, '(E11.4,",")',advance="no")dimag(mt_dm%main%ff(if1,if2,ias))
    enddo
    write(11,*)
  enddo
  write(11,*)
  close(11)
endif !exta outputs

if (.false.) then !write diagonal to a file
   open (11, file = 'dm-diago1.dat', status = 'replace')

   do if1=1, wfsize
     write(11,'(E11.4,",",E11.4)')dreal(mt_dm%main%ff(if1,if1,ias)),dimag(mt_dm%main%ff(if1,if1,ias))
   enddo
   close(11)

   open (11, file = 'dm-diago-copy.dat', status = 'replace')

   do if1=1, wfsize
     write(11,'(E11.4,",",E11.4)')dreal(dm_copy(if1,if1,ias)),dimag(dm_copy(if1,if1,ias))
   enddo
   close(11)


endif
!write(*,*)"dm ir sarakstīta failā"
!read(*,*)
! if (.true.) then  !test mt_dm if diagonal is real
!  !write(*,*)"*********apw************"
!   if1=1
!   do l1=0, input%groundstate%lmaxmat
!     do io1=1, apword (l1, is)
!         do m1=-l1, l1
!     !        lm1 = idxlm (l1, m1)
!            ! write(*,'(i3,". l=",i1," m=",i2," lm=",i2," dm=")')if1,l1,m1,lm1
!            ! write(*,*)"                   ",mt_dm%main%ff(if1,if1,ias)
!            ! if (dimag(mt_dm%main%ff(if1,if1,ias)).gt.1e-10) then

!             if (dimag(dm_copy(if1,if1,ias)).gt.1e-10) then

!               write(*,'(i3,". l=",i1," m=",i2," lm=",i2," dm=")')if1,l1,m1,lm1
! !              write(*,*)"                   ",mt_dm%main%ff(if1,if1,ias)
!               write(*,*)"                   ",dm_copy(if1,if1,ias)
!               write(*,*)"get_Fock_ex.f90 - complex element in the diagonal of mt_dm"
!               stop
!             endif
!         if1=if1+1
!         enddo
!      enddo
!   enddo
!   !write(*,*)"*********lo************"
!   if1=1
!   do ilo1 = 1, nlorb (is)
!     l1 = lorbl (ilo1, is)
!     do m1=-l1, l1
!       lm1 = idxlm (l1, m1)
! !      write(*,'(i3,". l=",i1," m=",i2," lm=",i2)')maxaa+if1,l1,m1,lm1
! !      write(*,*)"                   ",mt_dm%main%ff(maxaa+if1,maxaa+if1,ias)

! !            if (dimag(mt_dm%main%ff(if1,if1,ias)).gt.1e-10) then
!             if (dimag(dm_copy(if1,if1,ias)).gt.1e-10) then

!               write(*,'(i3,". l=",i1," m=",i2," lm=",i2," dm=")')if1,l1,m1,lm1
! !              write(*,*)"                   ",mt_dm%main%ff(if1,if1,ias)

!               write(*,*)"                   ",dm_copy(if1,if1,ias)


!               write(*,*)"get_Fock_ex.f90 - complex element in the diagonal of mt_dm"
!               stop
!             endif

!       if1=if1+1
!     enddo
!   enddo
! endif !test mt_dm


        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count apw+lo functions by each l !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do l1=0, input%groundstate%lmaxmat
  norb(l1)=apword (l1, is)
enddo
do ilo1 = 1, nlorb (is)
  l1 = lorbl (ilo1, is)
  norb(l1)=norb(l1)+1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create a density matrix map for locations of m=-l element !
! store basis function the same way as the matmap           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
maxnorb=maxval(norb)
allocate( matmap( maxnorb, 0:input%groundstate%lmaxmat ) )
allocate( basis( nrmt(is), maxnorb, 0:input%groundstate%lmaxmat ) )

tempindex=1
if1=1
do l1=0, input%groundstate%lmaxmat
  do io1=1, apword (l1, is)
    basis(:,tempindex(l1),l1)=apwfr (1:nrmt(is), 1, io1, l1, ias) 
    matmap(tempindex(l1),l1)=if1
    tempindex(l1)=tempindex(l1)+1 
    if1=if1+(2*l1+1)
  enddo
enddo

if1=maxaa+1
do ilo1 = 1, nlorb (is)
  l1 = lorbl (ilo1, is)
  basis(:,tempindex(l1),l1)=lofr (1:nrmt(is), 1, ilo1, ias)
  matmap(tempindex(l1),l1)=if1
  tempindex(l1)=tempindex(l1)+1
  if1=if1+(2*l1+1)
enddo


open (2, file = 'dm_map.dat', status = 'replace')
write(2,*)"***********MATMAP***********"
do l1=0, input%groundstate%lmaxmat
  write(2,*)l1,".",matmap(1:norb(l1),l1)
enddo
close(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! go through the density matrix and and construct vx_psi form valece orbitals !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!basis(ir,orb,l)
!ir    1:nrmt(is)     
!orb   1:norb(l)


integ=0d0
do l1=0, input%groundstate%lmaxmat
 ioo=1
 do io1=1,norb(l1)
    if1=matmap(io1,l1)
    do io2=1,norb(l1)
      if2=matmap(io2,l1)
      t1c=cmplx(0d0,0d0,8)
      do m1=0, 2*l1
        !t2c=mt_dm%main%ff(if1+m1,if2+m1,ias)
        t2c=dm_copy(if1+m1,if2+m1,ias)
        t1c=t1c+t2c
        !write(*,'("l=",I1," m=",I2,"|",I3,"|" ,I3,"| Re=",E11.4," Im=",E11.4)')l1,m1,if1+m1,if2+m1,&
        !        dreal(t2c),dimag(t2c)
      enddo

      !write(*,*)"           value=",t1c
      occ=dreal(t1c)
      vx_u_part=0d0

      call insum(nrmt(is),r(1:nrmt(is)),is,l,l1,occ,&
           basis(:,io1,l1)*r(1:nrmt(is)),&
           basis(:,io2,l1)*r(1:nrmt(is)),&
           u(1:nrmt(is)),&
           vx_u_part(1:nrmt(is)),.true.)

      integ(1:nrmt(is))=integ(1:nrmt(is))+vx_u_part(1:nrmt(is))
      ioo=ioo+1
    enddo !io2
  enddo !io1
enddo !l1
vx_u=vx_u+integ

endif!(mt_dm%maxnlo.ne.0).and.(mt_dm%maxaa.ne.0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! go through the core orbitals and and construct vx_psi !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integ=0d0
do ish=1, c_count(is)
  lpri=shell_l(ish)
  call insum(Ngrid,r,is,l,lpri,shell_occ(ish),u_all(:,ish),u_all(:,ish),u,vx_u_part,.false.)
  integ = integ + vx_u_part
enddo !ish


!open (2, file = 'vx_psi_core_apw.dat', status = 'replace')
! do ir = 1,Ngrid
!      write(2,*) r(ir), integ(ir)/r(ir), vx_psi(ir), integ(ir)/r(ir)+vx_psi(ir),
! end do
!close(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sum in vx_psi contribution from calence and core       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vx_u=vx_u+integ



end subroutine

subroutine  insum(Ngrid,r,is,l,lpri,occ,u1,u2,u,vx_u,mt)
  use modinteg
  use modinput, only:input
  use modbess, only: erfc_fit,nfit,msbesselic,msbesselkc
  use constants, only: zzero
  implicit none
  
  integer, intent(in) :: Ngrid,is, l, lpri
  real(8), intent(in) :: occ,r(Ngrid),u1(Ngrid),u2(Ngrid),u(Ngrid)
  logical, intent(in) :: mt
  real(8), intent(out) ::vx_u(Ngrid)
  
  complex(8) :: zbi(Ngrid,nfit), zbk(Ngrid,nfit),bi,bk,z
  complex(8) :: integc1(Ngrid),integc2(Ngrid),integc3(Ngrid),zf1(Ngrid),zf2(Ngrid)
  integer :: ifit, ir
  real (8)  :: wigner3j
  
  integer :: lpripri
  real(8) :: gc,integ1(Ngrid),integ2(Ngrid),integ3(Ngrid)
  logical :: erfc_kernel


  if (input%groundstate%hybrid%erfcapprox.ne."none") then
    erfc_kernel=.true.
  else
    erfc_kernel=.false.
  endif

  vx_u=0d0
    do lpripri=abs(l-lpri),l+lpri,2
  !    call wigner3j_list(l,lpri,lpripri,gc)
      gc=wigner3j(l,lpri,lpripri,0d0,0d0,0d0)
  
      gc=0.5d0*occ*gc**2
      !    write(*,*)"(l,l',l'') (",l,",",lpri,",",lpripri,")", " Gaunt_coef=",gc
      if (gc.ne.0d0) then
        if (erfc_kernel) then
          !!Generate set of complex bessel functions for current lpripri
          do ifit=1,nfit
            do ir = 1,Ngrid
              z=erfc_fit(ifit,2)*r(ir)
              call msbesselic (lpripri, z, bi)
              call msbesselkc (lpripri, z, bk)
              zbi(ir,ifit) = bi
              zbk(ir,ifit) = bk
            end do
          end do
          !! Bessel functoins ready
          
          integc3=zzero
          do ifit=1, nfit
            if (mt) then
              call integ_cf (Ngrid, is, zbi(:,ifit)*u2*u , zf1, mt_integw)
            else
              call integ_cf (Ngrid, is, zbi(:,ifit)*u2*u , zf1, atom_integw)
            endif
            integc1 = erfc_fit(ifit,2) * zbk(:,ifit) * zf1(:)
            if (mt) then
              call integ_cf (Ngrid, is, zbk(:,ifit)*u2*u, zf2, mt_integw)
            else
              call integ_cf (Ngrid, is, zbk(:,ifit)*u2*u, zf2, atom_integw)
            endif
            integc2= erfc_fit(ifit,2)* zbi(:,ifit) * (zf2(Ngrid)-zf2)
            
            integc3=integc3 + erfc_fit(ifit,1) * (integc1+integc2)
          enddo

          do ifit=2, nfit
            if (mt) then
              call integ_cf (Ngrid, is, conjg(zbi(:,ifit))*u2*u , zf1, mt_integw)
            else
              call integ_cf (Ngrid, is, conjg(zbi(:,ifit))*u2*u , zf1, atom_integw)
            endif
              integc1 = conjg(erfc_fit(ifit,2)) * conjg(zbk(:,ifit)) * zf1(:)
            if (mt) then
              call integ_cf (Ngrid, is, conjg(zbk(:,ifit))*u2*u, zf2, mt_integw)
            else
              call integ_cf (Ngrid, is, conjg(zbk(:,ifit))*u2*u, zf2, atom_integw)
            endif
              integc2= conjg(erfc_fit(ifit,2)) * conjg(zbi(:,ifit)) * (zf2(Ngrid)-zf2)
            
            integc3=integc3 + conjg(erfc_fit(ifit,1)) * (integc1+integc2)
          enddo
            ! do ir = 1,Ngrid
            !   write(*,*)ir,integc3(ir)
            ! enddo
            ! stop
          vx_u=vx_u - gc * u1 * dble(integc3) * dble(2*lpripri+1)
            ! !atomHF:
            ! call integC_BodesN_fun(Ngrid,r, tools, tools_info,1,  Bess_ik(:,k,lpripri+1,1) *u_all(:,ish)*u    ,integc1)
            ! integc1=integc1*rsfunC(k,2)*Bess_ik(:,k,lpripri+1,2)
            ! call integC_BodesN_fun(Ngrid,r, tools, tools_info,-1, Bess_ik(:,k,lpripri+1,2) *u_all(:,ish)*u    ,integc2)
            ! integc2=integc2*rsfunC(k,2)*Bess_ik(:,k,lpripri+1,1)
            ! rez=rez+rsfunC(k,1)*(integc1+integc2)
          

        else !Coulomb kernel
          if (mt) then      
            call integ_f(Ngrid,is,u2*u*r**lpripri,integ1,mt_integw)
          else
            call integ_f(Ngrid,is,u2*u*r**lpripri,integ1,atom_integw)
          endif
          integ1=integ1/r**(lpripri+1)
          if (mt) then
            call integ_f_rev(Ngrid,r,is,u2*u/r**(lpripri+1),integ2,mt_integw)
          else
            call integ_f_rev(Ngrid,r,is,u2*u/r**(lpripri+1),integ2,atom_integw)
          endif
          integ2=integ2*r**lpripri
          vx_u=vx_u + gc*u1*(-integ1-integ2)
        endif!erfc or Coulomb kernel


      endif !(gc.ne.0d0)
    enddo !lpripri
  

  end subroutine