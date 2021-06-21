!
!
!
! Copyright (C) 2021 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: zpotcoul
! !INTERFACE:
!
!
Subroutine pscharge (nr, nrmax, ld, r, igp0, gpc, jlgpr, ylmgp, sfacgp, &
& zn, zrhomt, zrhoir)
      Use modinput
! !USES:
      Use modmain
#ifdef USEOMP
      use omp_lib
#endif
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial points for each species (in,integer(nspecies))
!   nrmax  : maximum nr over all species (in,integer)
!   ld     : leading dimension of r (in,integer)
!   r      : radial mesh for each species (in,real(ld,nspecies))
!   igp0   : index of the shortest G+p-vector (in,integer)
!   gpc    : G+p-vector lengths (in,real(ngvec))
!   jlgpr  : spherical Bessel functions for evergy G+p-vector and muffin-tin
!            radius (in,real(0:lmaxvr+npsden+1,ngvec,nspecies))
!   ylmgp  : spherical harmonics of the G+p-vectors (in,complex(lmmaxvr,ngvec))
!   sfacgp : structure factors of the G+p-vectors (in,complex(ngvec,natmtot))
!   zn     : nuclear charges at the atomic centers (in,real(nspecies))
!   zrhomt : muffin-tin charge density (in,complex(lmmaxvr,nrmax,natmtot))
!   zrhoir : interstitial charge density (in,complex(ngrtot))
!   zvclmt : muffin-tin Coulomb potential (out,complex(lmmaxvr,nrmax,natmtot))
!   zvclir : interstitial Coulomb potential (out,complex(ngrtot))
!   zrho0  : G+p=0 term of the pseudocharge density (out,complex)
! !DESCRIPTION:
! Calculates soft pseudocharge basically repeating the first steps in zpotcoul. The result is stored in PSEDODENSITY.OUT. 
! The code reuses parts of zpotcoul.f90. 
! !REVISION HISTORY:
!   Created June 2021 (Andris)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: nr (nspecies)
      Integer, Intent (In) :: nrmax
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: r (ld, nspecies)
      Integer, Intent (In) :: igp0
      Real (8), Intent (In) :: gpc (ngvec)
      Real (8), Intent (In) :: jlgpr(0:input%groundstate%lmaxvr+ &
      & input%groundstate%npsden+1, ngvec, nspecies)
      Complex (8), Intent (In) :: ylmgp (lmmaxvr, ngvec)
      Complex (8), Intent (In) :: sfacgp (ngvec, natmtot)
      Real (8), Intent (In) :: zn (nspecies)
      Complex (8), Intent (In) :: zrhomt (lmmaxvr, nrmax, natmtot)
      Complex (8), Intent (In) :: zrhoir (ngrtot)
      Complex (8) :: zvclmt (lmmaxvr, nrmax, natmtot)
      Complex (8) :: zvclir (ngrtot)
! local variables
      Integer :: is, ia, ias, l, m, lm
      Integer :: ir, ig, ifg
      Real (8) :: fpo, t1, t2, t3
      Complex (8) zsum, zt1, zt2
! automatic arrays
      Real (8) :: rmtl (0:input%groundstate%lmaxvr+3, nspecies)
      Real (8) :: rl (nrmax, 0:input%groundstate%lmaxvr)
      Complex (8) vilm (lmmaxvr)
      Complex (8) qmt (lmmaxvr, natmtot)
      Complex (8) qi (lmmaxvr, natmtot)
      Complex (8) qilocal (lmmaxvr)
      Complex (8) zrp (lmmaxvr)
      real(8) :: vn(nrmax)
      real(8) :: third 
      parameter (third=0.3333333333333333333333d0)
      real(8), allocatable :: rho(:)
      complex(8), allocatable :: psrho(:)
 
#ifdef USEOMP
      integer ithr,nthreads,whichthread
#endif
      real(8), allocatable :: vdplmt(:,:,:), vdplir(:) 

! external functions
      Real (8) :: factnm
      External factnm
      fpo = fourpi / omega
! solve Poisson's equation for the isolated charge in the muffin-tin
      Do is = 1, nspecies
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias)
!$OMP DO
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call zpotclmt (input%groundstate%ptnucl, &
            & input%groundstate%lmaxvr, nr(is), r(:, is), zn(is), &
            & lmmaxvr, zrhomt(:, :, ias), zvclmt(:, :, ias))
            ias = idxas (ia, is)
         End Do
!$OMP END DO
!$OMP END PARALLEL
      End Do

! compute (R_mt)^l
      Do is = 1, nspecies
         rmtl (0, is) = 1.d0
         Do l = 1, input%groundstate%lmaxvr + 3
            rmtl (l, is) = rmtl (l-1, is) * rmt (is)
         End Do
      End Do
! compute the multipole moments from the muffin-tin potentials
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lm = 0
            Do l = 0, input%groundstate%lmaxvr
               t1 = dble (2*l+1) * rmtl (l+1, is) / fourpi
               Do m = - l, l
                  lm = lm + 1
                  qmt (lm, ias) = t1 * zvclmt (lm, nr(is), ias)
               End Do
            End Do
!            write(*,*) qmt (:, ias)
!            write(*,*)
         End Do
      End Do

! Fourier transform density to G-space and store in zvclir
      zvclir (:) = zrhoir (:)
      Call zfftifc (3, ngrid,-1, zvclir)
! find the multipole moments of the interstitial charge density
      qi (:, :) = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(qilocal,ig,ifg,zt1,t1,lm,t2,zt2,m,l,nthreads,whichthread,ithr) SHARED(input,gpc,sfacgp,zvclir,rmt,qi,ias,ngvec,is,ylmgp,jlgpr,zil,rmtl,igfft)
            qilocal=0d0
!$OMP DO    
#else
            qilocal=0d0
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
               If (gpc(ig) .Gt. input%structure%epslat) Then
                  zt1 = fourpi * zvclir (ifg) * sfacgp (ig, ias)
                  t1 = 1.d0 / (gpc(ig)*rmt(is))
                  lm = 0
                  Do l = 0, input%groundstate%lmaxvr
                     t2 = t1 * rmtl (l+3, is) * jlgpr (l+1, ig, is)
                     zt2 = t2 * zt1 * zil (l)
                     Do m = - l, l
                        lm = lm + 1
                        qilocal (lm) = qilocal (lm) + zt2 * conjg (ylmgp(lm, ig))
                     End Do
                  End Do
               Else
                  t1 = fourpi * y00 * rmtl (3, is) / 3.d0
                  qilocal (1) = qilocal (1) + t1 * zvclir (ifg)
               End If
            End Do
#ifdef USEOMP
!$OMP END DO
            nthreads=omp_get_num_threads()
            whichthread=omp_get_thread_num()
            do ithr=0,nthreads-1
              if (ithr.eq.whichthread) then
                qi(:,ias)=qi(:,ias)+qilocal(:) 
              endif
!$OMP BARRIER
            enddo
!$OMP END PARALLEL
#else
            qi(:,ias)=qi(:,ias)+qilocal(:)
#endif
         End Do
      End Do
! find the smooth pseudocharge within the muffin-tin whose multipoles are the
! difference between the real muffin-tin and interstitial multipoles
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lm = 0
            Do l = 0, input%groundstate%lmaxvr
! note the factor 2^N*N! is omitted because of reciprocal term in the
! form-factor
               t1 = factnm (2*(l+input%groundstate%npsden)+3, 2) / &
              & factnm (2*l+1, 2)
               Do m = - l, l
                  lm = lm + 1
                  zrp (lm) = (qmt(lm, ias)-qi(lm, ias)) * t1
               End Do
            End Do
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,ifg,zt1,t1,t2,zsum,m,l,lm,t3) SHARED(input,gpc,zvclir,rmt,qi,ias,ngvec,is,ylmgp,jlgpr,zil,rmtl,igfft,zrp,fpo,sfacgp)
!$OMP DO    
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
               If (gpc(ig) .Gt. input%structure%epslat) Then
                  zt1 = fpo * conjg (sfacgp(ig, ias))
                  t1 = gpc (ig) * rmt (is)
                  t2 = t1 ** (input%groundstate%npsden+1)
                  lm = 0
                  Do l = 0, input%groundstate%lmaxvr
                     lm = lm + 1
                     zsum = zrp (lm) * ylmgp (lm, ig)
                     Do m = - l + 1, l
                        lm = lm + 1
                        zsum = zsum + zrp (lm) * ylmgp (lm, ig)
                     End Do
                     t3 = jlgpr (input%groundstate%npsden+l+1, ig, is) / (t2*rmtl(l, is))
                     zvclir (ifg) = zvclir (ifg) + t3 * zt1 * zsum * conjg (zil(l))
                  End Do
               Else
                  t1 = fpo * y00 / factnm (2*input%groundstate%npsden+3, 2)
                  zvclir (ifg) = zvclir (ifg) + t1 * zrp (1)
               End If
            End Do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL 
#endif

         End Do
      End Do


       allocate(psrho(ngrtot))
       psrho=zvclir
       Call zfftifc (3, ngrid, 1, psrho)
       allocate(rho(ngrtot))
       rho=dble(psrho)*omega
       open(777,file='PSEUDODENSITY.OUT',action='write')
       write(777,*) "pseudodensity"
       write(777,*) "0.529177249"
       write(777,"(3F12.6)") input%structure%crystal%basevect(1:3, 1)
       write(777,"(3F12.6)") input%structure%crystal%basevect(1:3, 2)
       write(777,"(3F12.6)") input%structure%crystal%basevect(1:3, 3)
       write(777,*) natoms(1:nspecies)
       write(777,*) "Direct"  !"Cartesian"
       do is=1,nspecies
         do ia=1,natoms(is)
           write(777,"(3F13.6)") input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)  !    atposc(1,ia,is),atposc(2,ia,is),atposc(3,ia,is)
         enddo
       enddo
       write(777,*)
       write(777,"(3I8)") ngrid(1),ngrid(2), ngrid(3)
       write(777,"(E16.9)") (rho(ir),ir=1,ngrtot)
       close(777)
       deallocate(rho,psrho)
     
      Return
End Subroutine
!EOC
