subroutine WFprodcoul(ist1,wf1,ist2,wf2,prod,qlm)
      use modinput
      use mod_APW_LO
      use mod_atoms
      use mod_muffin_tin
      use mod_eigenvalue_occupancy
      use constants, only : zzero, zone, fourpi
      use mod_SHT
      use mod_Gvector, only : ngrtot
      use mod_eigensystem, only : WFType
      use weinert, only : poisson_and_multipoles_mt_yukawa2
      use mod_radial
      use mod_eigensystem, only : gntyyy
      use wigner3j_symbol, only : gaunt_yyy
! !USES:
! !DESCRIPTION:
! Evaluates a product of two WFs in the real space
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
      implicit none
      integer, intent(in) :: ist1,ist2
      type (WFType) :: wf1,wf2,prod
      complex(8), intent(out) :: qlm(lmmaxvr,natmtot)

      integer :: is,ia,ias
      integer :: l1,l2,m1,m2,lm1,lm2,io1,io2,if1,if2,ilo1,ilo2,lmmaxprod,lm,ir,l,m,ilo,io
      integer :: if1offset
      integer :: irad,jrad,irad1,irad2
      integer :: lmax
      integer :: blkstart,chunksize,iroffset
      integer, parameter :: blksize=64
      complex(8), allocatable :: factors(:), rho(:,:), fr(:)
      complex(8) :: mtmesh(ntpll,blksize), mtmesh1(ntpll,blksize), mtrlm1(blksize,lmmaxvr)
      complex(8) :: zfmt(lmmaxvr,nrmtmax),zfmttr(nrmtmax,lmmaxvr),zpot(lmmaxvr,nrmtmax)
      complex(8) :: zt
      real(8) :: ta,tb
      complex(8),allocatable :: H(:,:,:)
      real(8),external :: oldgaunt,oldwigner3j
 
      if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
      if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
 
! call timesec(ta)
      lmax= input%groundstate%lmaxvr

      allocate(H(lmmaxvr,maxradial,maxradial))

      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)

          H=0d0
if (.true.) then 
          if2=0
          do irad2=1,nradial(is)
            l2=lrad(irad2,is)
            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              if2=if2+1
              if1=0

              do irad1=1,nradial(is)
                l1=lrad(irad1,is)
!                write(*,*) 'L:',abs(l1-l2),min(lmax,l1+l2)
                do m1=-l1,l1
                  lm1=idxlm(l1,m1)
                  if1=if1+1
 
!                  do L=0,lmax           
                  do L=abs(l1-l2),min(lmax,l1+l2),2
!                    do M=-L,L
                    M=m1-m2
                    if ((M.le.L).and.(M.ge.-L)) then
                      LM=idxlm(L,M)
                      H(LM,irad1,irad2)=H(LM,irad1,irad2)+gntyyy(LM,lm2,lm1)*wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))!conjg(wf1%mt(if1,ist1,ias))*wf2%mt(if2,ist2,ias)
                 
!                     write(*,*) L,M,l1,m1,l2,m2,gntyyy(LM,lm1,lm2),oldgaunt(l2, l1, L, m2, m1, M),oldwigner3j (l1, l2, L, m1, -m2, M)  !gaunt_yyy (L, l2, l1, M, m2, m1) !gntyyy(LM,lm1,lm2)
                    endif 
!                    enddo
!                    do M=-L,L
!                      write(*,*) L,M,gntyyy(LM,lm1,lm2), oldwigner3j (l2, L, l1, -m2, M, m1) !oldwigner3j (l1, l2, L, m1, -m2, M)
!                      write(*,*) L,M,oldgaunt(l2, l1, L, m2, m1, M), oldwigner3j (l2, L, l1, -m2, M, m1) !oldwigner3j (l1, l2, L, m1, -m2, M)
!                    enddo
                  enddo
!                 read(*,*)


                enddo
              enddo

            enddo
          enddo
    
          zfmttr=0d0
          do irad2=1,nradial(is)
            do irad1=1,nradial(is)
!              do ir= 1,nrmt(is) !nrmt(is), nrmt(is) ! 1,nrmt(is)
!                zfmt(1:lmmaxvr,ir)=zfmt(1:lmmaxvr,ir)+H(1:lmmaxvr,irad1,irad2)*radialproducts(ir,irad1,irad2,ias)
!              enddo

              do lm=1,lmmaxvr !nrmt(is), nrmt(is) ! 1,nrmt(is)
                if (abs(H(lm,irad1,irad2)).gt.1d-10) then
!                  zfmttr(1:nrmt(is),lm)=zfmttr(1:nrmt(is),lm)+H(lm,irad1,irad2)*radial(1:nrmt(is),irad1,ias)*radial(1:nrmt(is),irad2,ias)
                   zfmttr(1:nrmt(is),lm)=zfmttr(1:nrmt(is),lm)+H(lm,irad1,irad2)*radialproducts(1:nrmt(is),irad1,irad2,ias)

!radialproducts(1:nrmt(is),irad1,irad2,ias)
                endif
              enddo

            enddo
          enddo
  
          do ir= 1,nrmt(is) 
            zfmt(1:lmmaxvr,ir)= zfmttr(ir,1:lmmaxvr)
          enddo

!H(1:lmmaxvr,irad1,irad2)*radialproducts(ir,irad1,irad2,ias)
!write(*,*) ist2,ist1,ias
!          do lm=1,16 !lmmaxvr
!            write(*,*) dble(zfmt(lm,nrmt(is)))
!          enddo
!write(*,*) '-----'
!stop
          
else
          chunksize=blksize
          do iroffset=1,nrmt(is),blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset

!-----------
! expand WFs in spherical harmonics
!-----------
            mtrlm1=0d0
            ias=idxas(ia,is)

            if1=0
            do irad = 1, nradial(is)
              l = lrad (irad,is)
              do m = - l, l
                if1=if1+1
                lm = idxlm (l, m)
                mtrlm1(1:chunksize,lm)=mtrlm1(1:chunksize,lm)+wf1%mt(if1,ist1,ias)*radial(iroffset:iroffset+chunksize-1,irad,ias)
              enddo
            enddo


            Call zgemm ('N', 'T', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
                       & mtrlm1, blksize, zzero, mtmesh, ntpll)
            do ir=1,chunksize
                 mtmesh(1:ntpll,ir)=mtmesh(1:ntpll,ir)*conjg(wf2%mtmesh(1:ntpll,iroffset+ir-1,ias,ist2))
            enddo
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, zfmt(1,iroffset),lmmaxvr) 
          enddo
endif
!          do lm=1,16 !lmmaxvr
!            write(*,*) dble(zfmt(lm,nrmt(is)))
!          enddo
!read(*,*)
          zpot=0d0
          call poisson_and_multipoles_mt_yukawa2( lmax, nrmt(is), spr(1:nrmt(is),is), zfmt, zpot, qlm(:,ias), is)
!          write(*,*) "poisson_mt"

!debug
!do lm=1,lmmaxvr
!  write(*,*) qlm(lm,ias)
!enddo
!stop


if(.false.)then
write(*,*) 'debug'
 open(11,file='mt_test.dat',status='replace')
 do ir=1, nrmt(1)
  write(11,*) spr(ir,is),dble(zfmt(1,ir)),dble(zpot(1,ir))  
 enddo
 close(11)
 stop
endif

          
          do iroffset=1,nrmt(is),blksize
            chunksize=blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset
            Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll, zpot(1,iroffset), lmmaxvr, zzero, mtmesh, ntpll)                        
            do ir=1,chunksize
              mtmesh(1:ntpll,ir)=mtmesh(1:ntpll,ir)*wf2%mtmesh(1:ntpll,iroffset+ir-1,ias,ist2)
            enddo
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, prod%mtrlm(1,iroffset,ias,1),lmmaxvr)
          enddo          
          
!prod%mtrlm(:,:,ias,1)
if(.false.)then
write(*,*) 'debug'
 open(11,file='mt_test.dat',status='replace')
 do ir=1, nrmt(1)
  write(11,*) spr(ir,1),dble(zpot(1,ir)), dble(prod%mtrlm(1,ir,1,1))
 enddo
 close(11)
 stop
endif
 
        enddo
      enddo
      



! call timesec(tb)
!write(*,*) tb-ta
 
end subroutine WFprodcoul
