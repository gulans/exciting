Module mod_radial
   implicit none
! radial functions: apwfr and lofr in the same place
! indices: radial point, shell, atom
!          ir, ifun, ias
   real(8), allocatable :: radial(:,:,:)  
   integer :: maxradial                   ! maximum number radial functions among all species
   integer, allocatable :: nradial(:)     ! number of radial functions for each species
   integer, allocatable :: lrad(:,:)      ! angular momentum corresponding to each radial function

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
   use mod_muffin_tin, only: nrmtmax, rmt
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
   allocate(lrad(maxradial,nspecies))
   lrad(:,:)=-1 ! for debugging purposes
   write(*,*) 'radial', nrmtmax*maxradial*natmtot*8/1d6,' Mb allocated'
   do is=1, nspecies
     do ia=1,natoms(is)   
       ias=idxas(ia,is)
       irad=0
       do l=0,input%groundstate%lmaxvr
         do io = 1, apword (l, is)
           irad=irad+1
           lrad(irad,is)=l
           radial(:,irad,ias)=apwfr(:,1,io,l,ias)
         enddo
       enddo

! local-orbital functions
       Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         irad=irad+1
         lrad(irad,is)=l
         radial(:,irad,ias)=lofr(:,1,ilo,ias)
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

End Module
