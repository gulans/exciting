subroutine WFprodcoul3(ia,is,ist1,wf1,ist2,wf2,FF,qlm)
!subroutine WFprodcoul3(ist1,wf1,ist2,wf2,FF,qlm)
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
      use mod_eigensystem, only : gntyyy,gntyyyT,gntlyy,gntyyl
!      use wigner3j_symbol, only : gaunt_yyy
! !USES:
! !DESCRIPTION:
! Evaluates a product of two WFs in the real space
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
      implicit none
      integer, intent(in) :: ist1,ist2 !,ia,is
      type (WFType) :: wf1,wf2
      complex(8), intent(out) :: qlm(lmmaxvr) !,natmtot)
      complex(8), intent(out) :: FF(maxpbf2,lmmaxvr) !,natmtot)
    

      integer :: is,ia
      integer :: ias
      integer :: l1,l2,m1,m2,lm1,lm2,io1,io2,if1,if2,ilo1,ilo2,lmmaxprod,lm,ir,l,m,ilo,io,if2offset
      integer :: if1offset
      integer :: irad,jrad,irad1,irad2,ipbf,radoffset
      integer :: lmax,llow,lhi
      integer :: blkstart,chunksize,iroffset,LMoffset
      integer, parameter :: blksize=64
      complex(8), allocatable :: factors(:), rho(:,:), fr(:)
      complex(8) :: mtmesh(ntpll,blksize), mtmesh1(ntpll,blksize), mtrlm1(blksize,lmmaxvr)
      complex(8) :: zfmt(lmmaxvr,nrmtmax),zfmttr(nrmtmax,lmmaxvr),zpot(lmmaxvr,nrmtmax),qlm2(lmmaxvr)!,H3(lmmaxvr)
      real(8) :: refmttr(nrmtmax,lmmaxvr),imfmttr(nrmtmax,lmmaxvr)
      complex(8) :: zt,zt2
      real(8) :: ta,tb
      complex(8),allocatable :: H(:,:,:), HH(:,:,:)
      complex(8),allocatable :: T(:,:,:)
      complex(8) :: F(maxpbf,lmmaxvr),U(maxpbf2,lmmaxvr)
 
 
! call timesec(ta)
!qlm=0d0 !(:,ias)
      lmax= input%groundstate%lmaxvr

      allocate(H(lmmaxvr,maxradial,maxradial))
      allocate(T(maxpbf2,lmmaxvr,lmmaxvr))

!      do is=1,nspecies
!        do ia=1,natoms(is)
          ias=idxas(ia,is)

          H=0d0
          F=0d0
          
          if2=0
          do irad2=1,nradial(is)
            l2=lrad(irad2,is)
            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              if2=if2+1
              if1=0

              do irad1=1,nradial(is)
                l1=lrad(irad1,is)
                do m1=-l1,l1
                  lm1=idxlm(l1,m1)
                  if1=if1+1
                  zt=wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))

                  do L=abs(l1-l2),min(lmax,l1+l2),2
                    M=m1-m2
                    if ((M.le.L).and.(M.ge.-L)) then
                      LM=idxlm(L,M)
                      H(LM,irad1,irad2)=H(LM,irad1,irad2)+gntyyl(L,lm2,lm1)*zt !wf1%mt(if1,ist1,ias)*conjg(wf2%mt(if2,ist2,ias))                
                    endif 
                  enddo
                enddo
              enddo
            enddo
          enddo

          do irad2=1,nradial(is)
            l2=lrad(irad2,is)
            do irad1=1,nradial(is)
              l1=lrad(irad1,is)
              do L=abs(l1-l2),min(lmax,l1+l2),2
                do M=-L,L
                  LM=idxlm(L,M)
                  F(1:npbf(L,ias),LM)=F(1:npbf(L,ias),LM)+H(LM,irad1,irad2)*uproducts(1:npbf(L,ias),L,irad1,irad2,ias) 
                enddo
              enddo
            enddo
          enddo

!-----------------

!          qlm(:,ias)=0d0
          qlm=0d0
          do l=0,lmax
            do m=-l,l
              lm=idxlm(l,m)
              zt=0d0
              do ipbf=1,npbf(l,ias)
                qlm(lm)=qlm(lm)+pbfmultipoles(ipbf,l,ias)*F(ipbf,lm)
              enddo
            enddo
          enddo

!-----------------


          T=0d0
          if2=0
          do irad2=1,nradial(is)
            U=0d0
            l2=lrad(irad2,is)
            radoffset=0
            do l1=0,lmax
              do m1=-l1,l1
                lm1=idxlm(l1,m1)
                do irad1=1,npbf(l1,ias)
                  U(1:npbf2(0,ias),lm1)=U(1:npbf2(0,ias),lm1)+F(irad1,lm1)*uproducts3(1:npbf2(0,ias),radoffset+irad1,irad2,ias)
                enddo
              enddo
              radoffset=radoffset+npbf(l1,ias)
            enddo
            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              if2=if2+1
              call zaxpy(npbf2(0,ias)*lmmaxvr, wf2%mt(if2,ist2,ias),U(1,1),1,T(1,1,lm2),1)
!              do lm1=1,lmmaxvr
!                T(1:npbf2(0,ias),lm1,lm2)=T(1:npbf2(0,ias),lm1,lm2)+U(1:npbf2(0,ias),lm1)*wf2%mt(if2,ist2,ias)
!              enddo
            enddo

          enddo

          do l2=0,lmax
            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              do l1=0,lmax
                do m1=-l1,l1
                  lm1=idxlm(l1,m1)
                  do L=abs(l1-l2),min(lmax,l1+l2),2
                    M=m1+m2 
                    if ((M.le.L).and.(M.ge.-L)) then
                      LM=idxlm(L,M)
                      FF(1:npbf2(0,ias),LM)=FF(1:npbf2(0,ias),LM)+T(1:npbf2(0,ias),lm1,lm2)*gntlyy(L,lm1,lm2)
                    endif
                  enddo
                enddo
              enddo 
            enddo
          enddo


          
 
!        enddo
!      enddo
      
deallocate(H,T)


! call timesec(tb)
!write(*,*) tb-ta
 
end subroutine WFprodcoul3
