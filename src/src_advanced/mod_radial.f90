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
! indices: ir, ipbf, l, ias
   real(8), allocatable :: productbasis(:,:,:,:) ! radialproducts are replaced with a linear combination of fiting functions
   integer :: maxpbf, maxpbfused                ! maximum number of product-basis functions among all atoms
   integer, allocatable :: npbf(:,:)      ! number of product-basis functions for every l and every atom

   real(8), allocatable :: radialmultipoles(:,:,:,:)
   real(8), allocatable :: radialpotential(:,:,:,:,:)

   real(8), allocatable :: uproducts(:,:,:,:,:) ! expansion coefficients for u-times-u in terms of product basis
   real(8), allocatable :: pbfpotential(:,:,:,:) ! radial Coulomb potential corresponding to productbasis
   real(8), allocatable :: pbfmultipoles(:,:,:) ! multipoles corresponding to productbasis
   

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
   radialmultipoles=0d0


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
             
           enddo
         enddo
       enddo
!test
       if (.false.) then
       l=3
       write(*,*) 'multipoles (l=',l, ') for ias=',ias
       do irad=1,nradial(is)
         write(*,*) radialmultipoles(l,irad,1,ias)
       enddo
       write(*,*) '-------------'
       endif
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

   subroutine init_productbasis
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
   integer :: ipbf,jpbf
   integer :: ir

   real(8), allocatable :: rl(:,:,:), ril1(:,:,:)
   real(8), allocatable :: pool(:,:,:), useful(:,:), norms(:,:), overlap(:,:)
   type :: arrtype
     real(8), allocatable :: chi(:,:,:)
   end type arrtype
   type(arrtype) :: pbf(natmtot)


   real(8) :: r2(nrmtmax)
   real(8) :: f1(nrmtmax), g1(nrmtmax), f2(nrmtmax), g2(nrmtmax)

   real(8) :: norm, dotp
   real(8) :: ta,tb
 
   real(8), parameter :: usefulness_thr=1d-8
   integer :: maxuseful,nuseful

   real(8), allocatable :: eval(:)
   real(8), allocatable :: work(:)
   integer :: lwork, liwork, info,liworkq
   integer, allocatable :: iwork(:)
   real(8) :: lworkq

   call timesec(ta)

   lmax=input%groundstate%lmaxapw
!   call release_vcoulradial
!  allocate(radialpotential(nrmtmax,0:lmax,maxradial,maxradial,natmtot))
!   allocate(radialmultipoles(0:lmax,maxradial,maxradial,natmtot))
!   write(*,*) 'vcoulradial', nrmtmax*maxradial*maxradial*(lmax+1)*natmtot*8/1d6,' Mb allocated'
!   write(*,*) 'radialmultipoles', maxradial*maxradial*(lmax+1)*natmtot*8/1d6,' Mb allocated'


!   allocate(rl(nrmtmax,0:lmax+2,nspecies))
!   allocate(ril1(nrmtmax,0:lmax+2,nspecies))

! safe estimate 
   maxpbf=maxradial*maxradial/4
   maxpbfused=0
   allocate(pool(nrmtmax,maxpbf,0:lmax))
   allocate(useful(nrmtmax,maxpbf))

   allocate(norms(maxpbf,0:lmax))
   allocate(overlap(maxpbf,maxpbf))
   allocate(eval(maxpbf))
! diagonalisation size query
   lwork=-1
   call dsyevd("V","U",maxpbf,overlap,maxpbf,eval,lworkq,lwork,liworkq,liwork,info)
   lwork=int(lworkq)
   liwork=liworkq
   write(*,*) info
   write(*,*) lwork, liwork
   allocate(work(lwork))
   allocate(iwork(liwork))
   
   call release_productbasis
!   allocate(productbasis(nrmtmax,maxpbf,0:lmax,natmtot))
   allocate(npbf(0:lmax,natmtot))
   npbf=0
   maxuseful=0
   
   do is=1, nspecies
     r2(1:nrmt(is))=spr(1:nrmt(is),is)*spr(1:nrmt(is),is)
     do ia=1,natoms(is)   
       ias=idxas(ia,is)
       do jrad=1,nradial(is)
         do irad=jrad,nradial(is)

           do l=abs(lrad(irad,is)-lrad(jrad,is)),min(lrad(irad,is)+lrad(jrad,is),lmax),2
             npbf(l,ias)=npbf(l,ias)+1
             pool(1:nrmt(is),npbf(l,ias),l)=radial(1:nrmt(is),irad,ias)*radial(1:nrmt(is),jrad,ias)
             f1(1:nrmt(is))=pool(1:nrmt(is),npbf(l,ias),l)*pool(1:nrmt(is),npbf(l,ias),l)*r2(1:nrmt(is))
             call integ_v(nrmt(is), is, f1, norms(npbf(l,ias),l), mt_integw)
             if (norms(npbf(l,ias),l).le.0d0) then
               write(*,*) "negative norm in init_productbasis"
               write(*,*) "terminating now"
               stop
             endif
             norms(npbf(l,ias),l)=sqrt(norms(npbf(l,ias),l))
             pool(1:nrmt(is),npbf(l,ias),l)=pool(1:nrmt(is),npbf(l,ias),l)/norms(npbf(l,ias),l)
           enddo
         enddo
       enddo
!       write(*,*) "ias=",ias
       do l=0,lmax
!         write(*,*) "l, npbf(l): ",l,npbf(l,ias)
! normalise PBFs
         do jpbf=1,npbf(l,ias)
!           overlap(ipbf,jpbf)=1d0
           do ipbf=jpbf,npbf(l,ias)
             f1(1:nrmt(is))=pool(1:nrmt(is),ipbf,l)*pool(1:nrmt(is),jpbf,l)*r2(1:nrmt(is))
             call integ_v(nrmt(is), is, f1, overlap(ipbf,jpbf), mt_integw)
             overlap(jpbf,ipbf)=overlap(ipbf,jpbf)
           enddo
         enddo
         call dsyevd("V","U",npbf(l,ias),overlap,maxpbf,eval,work,lwork,iwork,liwork,info)
         if (info.ne.0) then
           write(*,*) "Product-basis diagonalisation failed with info=",info
           write(*,*) "terminating now"
           stop
         endif
!         ipbf
         do ipbf=npbf(l,ias),1,-1
!           useful(npbf(l,ias)-ipbf+1)
!          write(*,*) ipbf,eval(ipbf)
           if (eval(ipbf).lt.usefulness_thr) exit
         enddo
!         write(*,*) "l,Nuseful",l,npbf(l,ias)-ipbf

         nuseful=npbf(l,ias)-ipbf
         call dgemm("N","N",nrmt(is),nuseful,npbf(l,ias),1d0,pool(1,1,l),nrmtmax,overlap(1,ipbf+1),maxpbf,0d0,useful,nrmtmax)
         npbf(l,ias)=nuseful

!         do ir=1, nrmt(is)
!           write(*,*) spr(ir,is),useful(ir,nuseful)
!         enddo
!         read(*,*)

!renormalise
!         norms(:,l)=0d0
         do ipbf=1,nuseful
           f1(1:nrmt(is))=useful(1:nrmt(is),ipbf)*useful(1:nrmt(is),ipbf)*r2(1:nrmt(is))
           norm=0d0
           call integ_v(nrmt(is), is, f1, norm, mt_integw)
           if (norm.le.0d0) then
             write(*,*) "negative norm in init_productbasis"
             write(*,*) "terminating now"
             stop
           endif
!           write(*,*) "ipbf",ipbf, "norm**2 = ", norm
           norm=sqrt(norm)
           pool(1:nrmt(is),ipbf,l)=useful(1:nrmt(is),ipbf)/norm
         enddo

!orthogonality test
         if (.false.) then
           do jpbf=1,nuseful
             do ipbf=1,nuseful
               dotp=0d0
               f1(1:nrmt(is))=pool(1:nrmt(is),ipbf,l)*pool(1:nrmt(is),jpbf,l)*r2(1:nrmt(is))
               call integ_v(nrmt(is), is, f1, dotp, mt_integw)
               write(*,*) "ipbf,jpbf",ipbf,jpbf, "norm**2 = ", dotp
             enddo
             read(*,*)
           enddo
         endif

         

!         
!         read(*,*)
         
!         productbasis(1:nrmt(is),1:npbf(l,ias),l,ias)=pool(1:nrmt(is),1:npbf(l,ias),l)
!         write(*,*) 'productbasis',l,ias,sum(productbasis(1:nrmt(is),1:npbf(l,ias),l,ias))
       enddo
       maxuseful=npbf(0,ias)
       do l=1,lmax
         maxuseful=max(npbf(l,ias),maxuseful)
       enddo
       maxpbfused=max(maxpbfused,maxuseful)
       write(*,*) maxuseful,ias
       allocate(pbf(ias)%chi(nrmtmax,maxuseful,0:lmax))
       
       do l=0,lmax
         pbf(ias)%chi(1:nrmt(is),1:npbf(l,ias),l)=pool(1:nrmt(is),1:npbf(l,ias),l)
!         write(*,*) 'pbf',l,ias,sum(pbf(ias)%chi(1:nrmt(is),1:npbf(l,ias),l))
       enddo
!       write(*,*)

     enddo
   enddo

   write(*,*) "maxpbf, maxpbfused", maxpbf,maxpbfused
   maxpbf=maxpbfused
   allocate(productbasis(nrmtmax,maxpbf,0:lmax,natmtot))
   write(*,*) 'productbasis', nrmtmax*maxpbf*(lmax+1)*natmtot*8/1d6,' Mb allocated'

   do is=1,nspecies
     do ia=1,natoms(is)
       ias=idxas(ia,is)
       do l=0,lmax
         productbasis(1:nrmt(is),1:npbf(l,ias),l,ias)=pbf(ias)%chi(1:nrmt(is),1:npbf(l,ias),l)
       enddo       
     enddo
   enddo
!   stop
 
   do ias=1,natmtot
!     do l=0,lmax
!       productbasis(1:nrmtmax,1:npbf(l,ias),l,ias)=pbf(ias)%chi(1:nrmtmax,1:npbf(l,ias),l)
!     enddo
     deallocate(pbf(ias)%chi)
   enddo  
   deallocate(pool)
   call timesec(tb)
   write(*,*) "init_productbasis:",tb-ta
!   stop
   end subroutine init_productbasis

   subroutine release_productbasis
   implicit none

   if (allocated(productbasis)) deallocate(productbasis)
   if (allocated(npbf)) deallocate(npbf)

   end subroutine release_productbasis

   subroutine init_uproducts
   use modinput
!   use mod_apw_lo, only: apwfr, lofr, apword, nlorb, lorbl
   use mod_atoms, only: idxas, natmtot, nspecies, natoms, spr
   use mod_muffin_tin, only: nrmtmax, nrmt
   use constants, only : fourpi
   use modinteg
   implicit none

   

   integer :: is, ia, ias
   integer :: l, m, lm, io, ilo
   integer :: irad, jrad
   integer :: lmax
   integer :: ipbf,jpbf
   integer :: ir

   real(8), allocatable :: pool(:,:,:), useful(:,:), norms(:,:), overlap(:,:)
   type :: arrtype
     integer, allocatable :: chi(:,:,:)
   end type arrtype
   type(arrtype) :: pbf(natmtot)


   real(8) :: r2(nrmtmax)
   real(8) :: f1(nrmtmax), g1(nrmtmax), f2(nrmtmax), g2(nrmtmax)

   real(8) :: norm, dotp, integ
   real(8) :: ta,tb

   lmax=input%groundstate%lmaxapw 
   call release_uproducts
   allocate(uproducts(maxpbfused,0:lmax,maxradial,maxradial,natmtot))
   write(*,*) 'uproducts', maxpbfused*maxradial*maxradial*(lmax+1)*natmtot*8/1d6,' Mb allocated'
   uproducts=0d0

   do is=1,nspecies
     r2(1:nrmt(is))=spr(1:nrmt(is),is)*spr(1:nrmt(is),is)
     do ia=1,natoms(is)
       ias=idxas(ia,is)
       do jrad=1,nradial(is)
         do irad=1,nradial(is)
           f2(1:nrmt(is))=radial(1:nrmt(is),irad,ias)*radial(1:nrmt(is),jrad,ias)*r2(1:nrmt(is))
           do l=abs(lrad(irad,is)-lrad(jrad,is)),min(lrad(irad,is)+lrad(jrad,is),lmax),2
             do ipbf=1,npbf(l,ias)
               f1(1:nrmt(is))=productbasis(1:nrmt(is),ipbf,l,ias)*f2(1:nrmt(is))
               call integ_v(nrmt(is), is, f1, integ, mt_integw) 
               uproducts(ipbf,l,irad,jrad,ias)=integ
             enddo
           enddo

         enddo
       enddo
     enddo
   enddo

!   f1=0d0
!   do ir=1,nrmt(is)
!     do ipbf=1,npbf(0,1)
!       f1(1:nrmt(1))=f1(1:nrmt(1))+productbasis(1:nrmt(1),ipbf,0,1)*uproducts(ipbf,0,1,1,1)
!     enddo
!   do ir=1,nrmt(1)
!     write(*,*) spr(ir,1),f1(ir),radial(ir,1,1)*radial(ir,1,1)
!   enddo
 !  enddo

!   write(*,*) 'done'
!   stop

   end subroutine init_uproducts

   subroutine release_uproducts
   implicit none
   
   if (allocated(uproducts)) deallocate(uproducts)

   end subroutine release_uproducts
   
   subroutine init_pbfpotential
   use modinput
!   use mod_apw_lo, only: apwfr, lofr, apword, nlorb, lorbl
   use mod_atoms, only: idxas, natmtot, nspecies, natoms, spr
   use mod_muffin_tin, only: nrmtmax, nrmt
   use constants, only : fourpi
   use modinteg
   implicit none

   

   integer :: is, ia, ias
   integer :: l, m, lm, io, ilo
   integer :: irad, jrad
   integer :: lmax
   integer :: ipbf,jpbf
   integer :: ir

   real(8), allocatable :: rl(:,:,:), ril1(:,:,:)
!   real(8), allocatable :: pool(:,:,:), useful(:,:), norms(:,:), overlap(:,:)
!   type :: arrtype
!     integer, allocatable :: chi(:,:,:)
!   end type arrtype
!   type(arrtype) :: pbf(natmtot)


   real(8) :: r2(nrmtmax), ri(nrmtmax)
   real(8) :: f1(nrmtmax), g1(nrmtmax), f2(nrmtmax), g2(nrmtmax)

   real(8) :: norm, dotp, integ, t1
   real(8) :: ta,tb

   lmax=input%groundstate%lmaxapw 
   call release_pbfpotential
   allocate(pbfpotential(nrmtmax,maxpbfused,0:lmax,natmtot))
   write(*,*) 'pbfpotential', nrmtmax*maxpbfused*(lmax+1)*natmtot*8/1d6,' Mb allocated'
   pbfpotential=0d0
   allocate(pbfmultipoles(maxpbfused,0:lmax,natmtot))
   write(*,*) 'pbfmultipoles', maxpbfused*(lmax+1)*natmtot*8/1d6,' Mb allocated'

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
       do l=0,lmax
         t1 = fourpi/(2*l+1)
         do ipbf=1,npbf(l,ias)
           f1(1:nrmt(is))=rl(1:nrmt(is),l+2,is)*productbasis(1:nrmt(is),ipbf,l,ias)   !radialproducts(1:nrmt(is),irad,jrad,ias)
           call integ_f (nrmt(is), is, f1, g1, mt_integw)
           f1(1:nrmt(is))=g1(1:nrmt(is))*ril1(1:nrmt(is),l,is)

           f2(1:nrmt(is))=ril1(1:nrmt(is),l,is)*rl(1:nrmt(is),2,is)*productbasis(1:nrmt(is),ipbf,l,ias)
           call integ_f (nrmt(is), is, f2, g2, mt_integw)
           f2(1:nrmt(is))=(g2(nrmt(is))-g2(1:nrmt(is)))*rl(1:nrmt(is),l,is)

           pbfmultipoles(ipbf,l,ias)= spr(nrmt(is),is)**(l+1) *(f1(nrmt(is))+f2(nrmt(is)))
           pbfpotential(1:nrmt(is),ipbf,l,ias)=t1*((f1(1:nrmt(is))+f2(1:nrmt(is))) - (f1(nrmt(is))+f2(nrmt(is)))*rl(1:nrmt(is),l,is)/rl(nrmt(is),l,is))
         enddo
       enddo

!test
!       if (.false.) then
!       l=3
!       write(*,*) '*multipoles (l=',l, ') for ias=',ias
!       do irad=1,nradial(is)
!         t1=0d0
!         do ipbf=1,npbf(l,ias)
!           t1=t1+pbfmultipoles(ipbf,l,ias)*uproducts(ipbf,l,irad,1,ias)
!         enddo  
!         write(*,*) t1 

!       enddo
!       write(*,*) '-------------'
!       endif

if (.false.) then       
       write(*,*) ias
       do jrad=1,nradial(is)
         do irad=1,nradial(is)
           do l=abs(lrad(irad,is)-lrad(jrad,is)),min(lrad(irad,is)+lrad(jrad,is),lmax),2
             t1=0d0
             do ipbf=1,npbf(l,ias)
               t1=t1+pbfmultipoles(ipbf,l,ias)*uproducts(ipbf,l,irad,jrad,ias)
             enddo
!             if (abs(t1-radialmultipoles(l,irad,jrad,ias)).gt.1d-10) then
!               write(*,*) l,irad,jrad,t1,radialmultipoles(l,irad,jrad,ias)
!             endif
             radialmultipoles(l,irad,jrad,ias)=t1
           enddo 
         enddo
       enddo
endif
!radialmultipoles(l,irad,irad,ias)
     enddo
   enddo
   deallocate(rl,ril1)




   write(*,*) 'pbfpotentials done'
!   stop

   end subroutine init_pbfpotential 

   subroutine release_pbfpotential 
   implicit none
   
   if (allocated(pbfpotential)) deallocate(pbfpotential)
   if (allocated(pbfmultipoles)) deallocate(pbfmultipoles)

   end subroutine release_pbfpotential

End Module
