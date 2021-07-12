!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepvnlk2 (ikp, vnlcv, vnlvv)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (Out) :: vnlcv (ncrmax, natmtot, nstsv)
      Complex (8), Intent (Out) :: vnlvv (nstsv, nstsv)
! local variables
      Integer :: ngknr, ik, ist1, ist2, ist3
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm, ir
      Integer :: nrc, iq, ig, iv (3), igq0
      Real (8) :: v (3), cfq, ta,tb
      Complex (8) zrho01, zrho02, zt1, zt2
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Integer, Allocatable :: igkignr (:)
      Real (8), Allocatable :: vgklnr (:, :)
      Real (8), Allocatable :: vgkcnr (:, :)
      Real (8), Allocatable :: gkcnr (:)
      Real (8), Allocatable :: tpgkcnr (:, :)
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalsvp (:)
      Real (8), Allocatable :: evalsvnr (:)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: sfacgknr (:, :)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :, :, :)
      Complex (8), Allocatable :: wfir1 (:, :, :)
      Complex (8), Allocatable :: wfir2 (:, :, :)
      Complex (8), Allocatable :: wfcr1 (:, :, :)
      Complex (8), Allocatable :: wfcr2 (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
      Complex (8), Allocatable :: zvcltp (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
      type (WFType) :: wf1,wf2,prod,pot
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp
! allocate local arrays
      if (.not.allocated(gntyyy)) call gengntyyy
      Allocate (igkignr(ngkmax))
      Allocate (vgklnr(3, ngkmax))
      Allocate (vgkcnr(3, ngkmax))
      Allocate (gkcnr(ngkmax))
      Allocate (tpgkcnr(2, ngkmax))
      Allocate (vgqc(3, ngvec))
      Allocate (tpgqc(2, ngvec))
      Allocate (gqc(ngvec))
      Allocate (jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, ngvec, nspecies))
      Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
!      Allocate (evalsvp(nstsv))
      Allocate (evalsvnr(nstsv))
      Allocate (sfacgknr(ngkmax, natmtot))
!      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
!      Allocate (evecfv(nmatmax, nstfv))
!      Allocate (evecsv(nstsv, nstsv))
      Allocate (ylmgq(lmmaxvr, ngvec))
      Allocate (sfacgq(ngvec, natmtot))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir1(ngrtot, nspinor, nstsv))
      Allocate (wfir2(ngrtot, nspinor, nstsv))
      Allocate (wfcr1(lmmaxvr, nrcmtmax, 2))
      Allocate (wfcr2(lmmaxvr, nrcmtmax, 2))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zvclir(ngrtot))
      Allocate (zvcltp(lmmaxvr, nrcmtmax))
      Allocate (zfmt(lmmaxvr, nrcmtmax))

      call WFInit(wf1)
      call WFInit(wf2)
      call WFInit(prod)
      call WFInit(pot)

      call genWF(ikp,wf1)
      call genWFinMT(wf1)
      call genWFonMesh(wf1)


! factor for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the nuclear charges to zero
      zn (:) = 0.d0
      vnlcv (:, :, :) = 0.d0
      vnlvv (:, :) = 0.d0
! get the eigenvalues/vectors from file for input k-point
!      Call getevalsv (vkl(:, ikp), evalsvp)
!      Call getevecfv (vkl(:, ikp), vgkl(:, :, :, ikp), evecfv)
!      Call getevecsv (vkl(:, ikp), evecsv)
! find the matching coefficients
!      Call match (ngk(1, ikp), gkc(:, 1, ikp), tpgkc(:, :, 1, ikp), sfacgk(:, :, 1, ikp), apwalm)
! calculate the wavefunctions for all states for the input k-point
!      Call genwfsv (.False., ngk(1, ikp), igkig(:, 1, ikp), evalsvp, apwalm, evecfv, evecsv, wfmt1, wfir1)

      allocate(pot%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
      allocate(pot%ir(ngrtot,1))


do ist2= 1,nstfv 
! start loop over non-reduced k-point set
      zvclir=0d0
      zvclmt=0d0

      Do ik = 1, nkptnr
         
! generate G+k-vectors
         Call gengpvec (vklnr(:, ik), vkcnr(:, ik), ngknr, igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-points
         Call getevalsv (vklnr(:, ik), evalsvnr)
!         Call getevecfv (vklnr(:, ik), vgklnr, evecfv)
!         Call getevecsv (vklnr(:, ik), evecsv)
! generate the structure factors
         Call gensfacgp (ngknr, vgkcnr, ngkmax, sfacgknr)
! find the matching coefficients
!         Call match (ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm)
! determine q-vector
         iv (:) = ivk (:, ikp) - ivknr (:, ik)
         iv (:) = modulo (iv(:), input%groundstate%ngridk(:))
         iq = iqmap (iv(1), iv(2), iv(3))
         v (:) = vkc (:, ikp) - vkcnr (:, ik)
         Do ig = 1, ngvec
! determine G+q vectors
            vgqc (:, ig) = vgc (:, ig) + v (:)
! G+q-vector length and (theta, phi) coordinates
            Call sphcrd (vgqc(:, ig), gqc(ig), tpgqc(:, ig))
! spherical harmonics for G+q-vector
            Call genylm (input%groundstate%lmaxvr, tpgqc(:, ig), ylmgq(:, ig))
         End Do
! structure factors for G+q
         Call gensfacgp (ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
         Call findigp0 (ngvec, gqc, igq0)
         sfacgq0 (:) = sfacgq (igq0, :)
! compute the required spherical Bessel functions
         lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
!         Call genjlgpr (lmax, gqc, jlgqr)
         Call genjlgpr (lmax, gc, jlgqr)
         Call genjlgq0r (gqc(igq0), jlgq0r)
! calculate the wavefunctions for occupied states

         call genWF(ik,wf2)
         call genWFinMT(wf2)
         call genWFonMesh(wf2)

         Do ist3 = 1, nstsv
            If (evalsvnr(ist3) .Lt. efermi) Then

! calculate the complex overlap density

call timesec(ta)
                     call WFprod(ist3,wf2,ist2,wf1,prod)
call timesec(tb)
!write(*,*) 'WFprod',tb-ta
call timesec(ta)
                     Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                    & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prod%ir(:,1), zrho01)
call timesec(tb)
!write(*,*) 'zrhogp',tb-ta
call timesec(ta)
                     prod%ir(:,1)=prod%ir(:,1)-zrho01
                     prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00
call timesec(tb)
!write(*,*) 'remove average',tb-ta
call timesec(ta)
! calculate the Coulomb potential
                     Call zpotcoul (nrcmt, nrcmtmax, nrcmtmax, rcmt, &
                    & igq0, gqc, jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:,:,:,1), &
                    & prod%ir(:,1), pot%mtrlm(:,:,:,1), pot%ir(:,1), zrho02)


!do ir=1,nrmt(1)
 ! write(*,*) spr(ir,1),dble(pot%mtrlm(1,ir,1,1))
!enddo
!read(*,*)

                     Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                    & igq0), sfacgq0, pot%mtrlm(:,:,:,1), pot%ir(:,1), zrho01)

                     pot%ir(:,1)=pot%ir(:,1)-zrho01
                     pot%mtrlm(1,:,:,1)=pot%mtrlm(1,:,:,1)-zrho01/y00

                     call genWFonMeshOne(pot)
                     pot%ir=conjg(pot%ir)
                     pot%mtmesh=conjg(pot%mtmesh)
                     call WFprodrs(1,pot,ist3,wf2,prod)

                     zvclir=zvclir+wkptnr(ik)*prod%ir(:,1)
                     zvclmt=zvclmt+wkptnr(ik)*prod%mtrlm(:,:,:,1)

!write(*,*) sum(zvclir)                     

!write(*,*) 'loop begins'
!                     do ist1=1,nstfv
!                       zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),prod%mtrlm(:,:,:,1), wf1%ir(:,ist1), prod%ir(:,1))

!                       zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),zvclmt, wf1%ir(:,ist1), zvclir)


!                       call WFprodrs(ist3,wf2,ist1,wf1,prod)

!                       call WFprodrs(ist3,wf2,ist1,wf1,prod)
!                       Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
!                      & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prod%ir(:,1), zrho01)
!                       prod%ir(:,1)=prod%ir(:,1)-zrho01
!                       prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00

!                       zt1 = zfinp (.True., prod%mtrlm(:,:,:,1), pot%mtrlm(:,:,:,1), prod%ir(:,1), pot%ir(:,1))
!                       vnlvv (ist1, ist2) = vnlvv (ist1, ist2) - wkptnr(ik)*zt1
!                     enddo



!                     Do ist1 = 1, nstsv
! calculate the complex overlap density
!                            call WFprod(ist3,wf2,ist1,wf1,prod)

 !                           Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
 !                          & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prod%ir(:,1), zrho01)

!                            prod%ir(:,1)=prod%ir(:,1)-zrho01
!                            prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00

!                           zt1 = zfinp (.True., prod%mtrlm(:,:,:,1), pot%mtrlm(:,:,:,1), prod%ir(:,1), pot%ir(:,1))
!-------------------------------------------------------------------
! compute the density coefficient of the smallest G+q-vector
!                           vnlvv (ist1, ist2) = vnlvv (ist1, ist2) - &
!                          & (wkptnr(ik)*zt1)
!                     End Do

! end loop over ist3
            End If
         End Do
! end loop over non-reduced k-point set
      call WFRelease(wf2)

      End Do
!write(*,*) sum(zvclir)
                     do ist1=1,nstfv
!                       zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),prod%mtrlm(:,:,:,1), wf1%ir(:,ist1), prod%ir(:,1))

                       zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),zvclmt, wf1%ir(:,ist1), zvclir)


!                       call WFprodrs(ist3,wf2,ist1,wf1,prod)

!                       call WFprodrs(ist3,wf2,ist1,wf1,prod)
!                       Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
!                      & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prod%ir(:,1), zrho01)
!                       prod%ir(:,1)=prod%ir(:,1)-zrho01
!                       prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00

!                       zt1 = zfinp (.True., prod%mtrlm(:,:,:,1), pot%mtrlm(:,:,:,1), prod%ir(:,1), pot%ir(:,1))
                       vnlvv (ist1, ist2) = vnlvv (ist1, ist2) - zt1
                     enddo

!      do ist1=1,nstfv
!        zt1 = zfinp (.True., wf1%mtrlm(:,:,:,ist1),zvclmt, wf1%ir(:,ist1), zvclir)
!        vnlvv (ist1, ist2) = vnlvv (ist1, ist2) - wkptnr(ik)*zt1       
!      enddo

enddo ! end loop over ist2

      Deallocate (igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
!      Deallocate (evalsvp)
      Deallocate (evalsvnr) 
!     Deallocate (evecfv, evecsv)
!      Deallocate (apwalm)
      Deallocate (sfacgknr, ylmgq, sfacgq)
      Deallocate (wfmt1, wfmt2, wfir1, wfir2, wfcr1, wfcr2)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir, zvcltp, zfmt)
      call WFRelease(wf1)
      call WFRelease(prod)
      call WFRelease(pot)

write(*,*) 'WFRelease done'
      deallocate(gntyyy)
      Return
End Subroutine
!EOC
