Module mod_radial
   implicit none
! radial functions: apwfr and lofr in the same place
! indices: radial point, shell, atom
!          ir, ifun, ias
   real(8), allocatable :: radial(:,:,:)  
   integer :: maxradial                   ! maximum number radial functions among all species
   integer, allocatable :: nradial(:)     ! number of radial functions for each species
   integer, allocatable :: lrad(:,:)      ! angular momentum corresponding to each radial function

   real(8), allocatable :: radialproducts(:,:,:,:)
   real(8), allocatable :: radialmultipoles(:,:,:,:)
   real(8), allocatable :: radialpotential(:,:,:,:,:)

Contains

!
!
!
!BOP
! !ROUTINE: MTNullify
! !INTERFACE:
!
!
   subroutine init_radial
   use modinput
   use mod_APW_LO, only: apwfr, lofr, apword, nlorb, lorbl
   use mod_atoms, only: idxas, natmtot, nspecies, natoms
   use mod_muffin_tin, only: nrmtmax, nrmt
   implicit none

   integer :: is, ia, ias
   integer :: l, m, lm, io, ilo
   integer :: irad
   integer :: nlo,napw
   integer :: lmax

   lmax=input%groundstate%lmaxapw
   if (allocated(radial)) deallocate(radial)
   if (allocated(nradial)) deallocate(nradial)
   if (allocated(lrad)) deallocate(lrad)
! find the maximum number of radial functions per atom
   maxradial=0
   allocate(nradial(nspecies))

   do is=1, nspecies
     napw=0
     do l=0,lmax
       napw=napw+apword(l,is)
     enddo
     nradial(is) = nlorb(is)+napw
     maxradial = max(maxradial,nradial(is))
!     write(*,*) napw,nlorb(is),nradial(is)
   enddo 
!   write(*,*) maxradial

   allocate(radial(nrmtmax,maxradial,natmtot))
   write(*,*) 'radial', nrmtmax*maxradial*natmtot*8/1d6,' Mb allocated'
   allocate(lrad(maxradial,nspecies))
   lrad(:,:)=-1 ! for debugging purposes

   do is=1, nspecies
     do ia=1,natoms(is)   
       ias=idxas(ia,is)
       irad=0
       do l=0,input%groundstate%lmaxvr
         do io = 1, apword (l, is)
           irad=irad+1
           lrad(irad,is)=l
           radial(1:nrmt(is),irad,ias)=apwfr(1:nrmt(is),1,io,l,ias)
         enddo
       enddo

! local-orbital functions
       Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         irad=irad+1
         lrad(irad,is)=l
         radial(1:nrmt(is),irad,ias)=lofr(1:nrmt(is),1,ilo,ias)
       End Do
     enddo
   enddo

   end subroutine init_radial

   subroutine release_radial
   implicit none

     if (allocated(radial)) deallocate(radial)
     if (allocated(nradial)) deallocate(nradial)
     if (allocated(lrad)) deallocate(lrad)

   end subroutine release_radial

   subroutine init_radial_products
   use modinput
!   use mod_APW_LO, only: apwfr, lofr, apword, nlorb, lorbl
   use mod_atoms, only: idxas, natmtot, nspecies, natoms
   use mod_muffin_tin, only: nrmtmax, nrmt
   implicit none

   integer :: is, ia, ias
   integer :: l, m, lm, io, ilo
   integer :: irad, jrad
   integer :: lmax

!  if (allocated(radialproducts)) deallocate(radialproducts)
   call release_radial_products
   allocate(radialproducts(nrmtmax,maxradial,maxradial,natmtot))
   write(*,*) 'radialproducts', nrmtmax*maxradial*maxradial*natmtot*8/1d6,' Mb allocated'

   do is=1, nspecies
     do ia=1,natoms(is)   
       ias=idxas(ia,is)
       do jrad=1,nradial(is)
         do irad=1,nradial(is)
           radialproducts(1:nrmt(is),irad,jrad,ias)=radial(1:nrmt(is),irad,ias)*radial(1:nrmt(is),jrad,ias)
         enddo
       enddo
     enddo
   enddo
   end subroutine init_radial_products
   

   subroutine release_radial_products
   implicit none

   if (allocated(radialproducts)) deallocate(radialproducts)

   end subroutine release_radial_products

   subroutine init_vcoulradial
   use modinput
!   use mod_APW_LO, only: apwfr, lofr, apword, nlorb, lorbl
   use mod_atoms, only: idxas, natmtot, nspecies, natoms, spr
   use mod_muffin_tin, only: nrmtmax, nrmt
   use constants, only : fourpi
   use modinteg
   implicit none

   integer :: is, ia, ias
   integer :: l, m, lm, io, ilo
   integer :: irad, jrad
   integer :: lmax

   real(8), allocatable :: rl(:,:,:), ril1(:,:,:)
   real(8) :: ri(nrmtmax)
   real(8) :: f1(nrmtmax), g1(nrmtmax), f2(nrmtmax), g2(nrmtmax)

   real(8) :: t1
   real(8) :: ta,tb

   call timesec(ta)

   lmax=input%groundstate%lmaxapw
   call release_vcoulradial
   allocate(radialpotential(nrmtmax,0:lmax,maxradial,maxradial,natmtot))
   allocate(radialmultipoles(0:lmax,maxradial,maxradial,natmtot))
   write(*,*) 'vcoulradial', nrmtmax*maxradial*maxradial*(lmax+1)*natmtot*8/1d6,' Mb allocated'
   write(*,*) 'radialmultipoles', maxradial*maxradial*(lmax+1)*natmtot*8/1d6,' Mb allocated'


   allocate(rl(nrmtmax,0:lmax+2,nspecies))
   allocate(ril1(nrmtmax,0:lmax+2,nspecies))


   do is=1, nspecies
     rl(1:nrmt(is),0,is) = 1d0          !r^l
     ri(1:nrmt(is)) = 1d0/spr(1:nrmt(is),is)
     ril1(1:nrmt(is),0,is) = 1d0/spr(1:nrmt(is),is) ! r^(-l-1)
     do l = 0, lmax+1
       rl(1:nrmt(is),l+1,is)=rl(1:nrmt(is),l,is)*spr(1:nrmt(is),is)
       ril1(1:nrmt(is),l+1,is)=ril1(1:nrmt(is),l,is)*ri(1:nrmt(is))
     enddo
   enddo


   do is=1, nspecies


     do ia=1,natoms(is)   
       ias=idxas(ia,is)
       do jrad=1,nradial(is)
         do irad=1,nradial(is)
           do l=abs(lrad(irad,is)-lrad(jrad,is)),min(lrad(irad,is)+lrad(jrad,is),lmax),2
             t1 = fourpi/(2*l+1)
            
             f1(1:nrmt(is))=rl(1:nrmt(is),l+2,is)*radialproducts(1:nrmt(is),irad,jrad,ias)
             call integ_f (nrmt(is), is, f1, g1, mt_integw)
             f1(1:nrmt(is))=g1(1:nrmt(is))*ril1(1:nrmt(is),l,is)

             f2(1:nrmt(is))=ril1(1:nrmt(is),l,is)*rl(1:nrmt(is),2,is)*radialproducts(1:nrmt(is),irad,jrad,ias)
             call integ_f (nrmt(is), is, f2, g2, mt_integw)
             f2(1:nrmt(is))=(g2(nrmt(is))-g2(1:nrmt(is)))*rl(1:nrmt(is),l,is)

             radialmultipoles(l,irad,jrad,ias)= spr(nrmt(is),is)**(l+1) *(f1(nrmt(is))+f2(nrmt(is)))
             radialpotential(1:nrmt(is),l,irad,jrad,ias)=t1*((f1(1:nrmt(is))+f2(1:nrmt(is))) - (f1(nrmt(is))+f2(nrmt(is)))*rl(1:nrmt(is),l,is)/rl(nrmt(is),l,is))

!  zvclmt (lm, :nr)=t1*((f1+f2) - (f1(nr)+f2(nr))*rl(:)/rl(nr))
             
           enddo
         enddo
       enddo
     enddo
   enddo
   deallocate(rl,ril1)
   call timesec(tb)
   write(*,*) "init_vcoulradial:", tb-ta

   end subroutine init_vcoulradial

   subroutine release_vcoulradial
   implicit none

   if (allocated(radialpotential)) deallocate(radialpotential)
   if (allocated(radialmultipoles)) deallocate(radialmultipoles)

   end subroutine release_vcoulradial


End Module
