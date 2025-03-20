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
      integer :: l1,l3,m1,m3,lm1,lm3,lm2,io1,io2,if1,if3,if1old,if3old,ilo1,ilo2,lmmaxprod,lm,ir,l,m,ilo,io
      integer :: lmax
      integer :: blkstart,chunksize,iroffset
      integer, parameter :: blksize=64
      complex(8), allocatable :: factors(:), rho(:,:), fr(:)
      complex(8) :: mtmesh(ntpll,blksize), mtmesh1(ntpll,blksize), mtrlm1(blksize,lmmaxvr)
      complex(8) :: zfmt(lmmaxvr,nrmtmax),zpot(lmmaxvr,nrmtmax)
      complex(8) :: zt
      real(8) :: ta,tb
 
      if (.not.allocated(prod%ir)) allocate(prod%ir(ngrtot,1))
      if (.not.allocated(prod%mtrlm)) allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
 
! call timesec(ta)
      lmax= input%groundstate%lmaxvr

      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)

 
          chunksize=blksize


          do iroffset=1,nrmt(is),blksize
            if (iroffset+blksize-1.gt.nrmt(is)) chunksize=nrmt(is)+1-iroffset

!-----------
! expand WFs in spherical harmonics
!-----------
            mtrlm1=0d0
            ias=idxas(ia,is)
            if1=0
! APW part
            do l=0,input%groundstate%lmaxvr
              do io = 1, apword (l, is)
                do m=-l,l
                  lm=idxlm(l,m)
                  if1=if1+1
                  mtrlm1(1:chunksize,lm)=mtrlm1(1:chunksize,lm)+wf1%mt(if1,ist1,ias)*apwfr(iroffset:iroffset+chunksize-1,1,io,l,ias)
!                  mtrlm2(lm,1:chunksize)=mtrlm2(lm,1:chunksize)+wf1%mt(if1,ist1,ias)*apwfr(iroffset:iroffset+chunksize-1,1,io,l,ias)
!                  mtrlm2(lm,1:chunksize)=mtrlm2(lm,1:chunksize)+wf1%mt(if1,ist1,ias)*apwfr(iroffset:iroffset+chunksize-1,1,io,l,ias)
                enddo
              enddo
            enddo

! local-orbital functions
            Do ilo = 1, nlorb (is)
              l = lorbl (ilo, is)
              Do m = - l, l
                if1=if1+1
                lm = idxlm (l, m)
                 mtrlm1(1:chunksize,lm)=mtrlm1(1:chunksize,lm)+wf1%mt(if1,ist1,ias)*lofr(iroffset:iroffset+chunksize-1,1,ilo,ias)
!                mtrlm2(lm,1:chunksize)=mtrlm2(lm,1:chunksize)+wf1%mt(if1,ist1,ias)*lofr(iroffset:iroffset+chunksize-1,1,ilo,ias)
!                 mtrlm2(lm,1:chunksize)=mtrlm2(lm,1:chunksize)+wf1%mt(if1,ist1,ias)*lofr(iroffset:iroffset+chunksize-1,1,ilo,ias)
              End Do
            End Do




!            Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
!                       & mtrlm2, lmmaxvr, zzero, mtmesh1, ntpll)
!            do lm=1,lmmaxvr
!              mtrlm1(1:chunksize,lm)=mtrlm2(lm,1:chunksize)
!            enddo
            Call zgemm ('N', 'T', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
                       & mtrlm1, blksize, zzero, mtmesh, ntpll)
!           Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
!                      & mtrlm2, lmmaxvr, mtmesh1, ntpll)

!           Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
!                      & wf1%mtrlm(1,iroffset,ias,ist1), lmmaxvr, zzero, mtmesh1, ntpll)
!           Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
!                      & mtrlm2, lmmaxvr, zzero, mtmesh2, ntpll)

!            Call zgemm ('N', 'N', ntpll, chunksize, lmmaxvr, zone, zbshthf, ntpll,  &
!                       & wf2%mtrlm(1,iroffset,ias,ist2), lmmaxvr, zzero, mtmesh2, ntpll)
            do ir=1,chunksize
!              do lm=1,ntpll
                 mtmesh(1:ntpll,ir)=mtmesh(1:ntpll,ir)*conjg(wf2%mtmesh(1:ntpll,iroffset+ir-1,ias,ist2)) !conjg(mtmesh2(lm,ir))   !conjg(wf1%mtmesh(lm,ir,ias,ist1))*wf2%mtmesh(lm,ir,ias,ist2)
!                 mtmesh(lm,ir)=wf1%mtmesh(lm,iroffset+ir-1,ias,ist1)*conjg(wf2%mtmesh(lm,iroffset+ir-1,ias,ist2)) !conjg(mtmesh2(lm,ir))   !conjg(wf1%mtmesh(lm,ir,ias,ist1))*wf2%mtmesh(lm,ir,ias,ist2)
!                 mtmesh(lm,ir)=conjg(wf1%mtmesh(lm,iroffset+ir-1,ias,ist1))*wf2%mtmesh(lm,iroffset+ir-1,ias,ist2) !conjg(mtmesh2(lm,ir))   !conjg(wf1%mtmesh(lm,ir,ias,ist1))*wf2%mtmesh(lm,ir,ias,ist2)
!              enddo
            enddo
            Call zgemm ('N', 'N', lmmaxvr, chunksize, ntpll, zone, zfshthf, lmmaxvr, mtmesh, ntpll, zzero, zfmt(1,iroffset),lmmaxvr) 
!prod%mtrlm(1,iroffset,ias,1) , lmmaxvr) ! Genshtmat3
          enddo
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
