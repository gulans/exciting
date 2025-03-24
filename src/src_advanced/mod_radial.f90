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
!   real(8), allocatable :: vcoulradial(:,:,:,:,:)

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
   write(*,*) 'radialproducts', nrmtmax*nrmtmax*maxradial*natmtot*8/1d6,' Mb allocated'

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


End Module
