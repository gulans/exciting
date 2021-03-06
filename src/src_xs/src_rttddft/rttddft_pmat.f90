! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! !REVISION HISTORY:
!
!  Created April 2019 (Ronaldo Rodrigues Pela)
! Reference: https://arxiv.org/abs/2102.02630

!> This module deals with the momentum matrix considering as basis (L)APW+lo
module rttddft_pmat
  use modmpi
  use modxs, only: ripaa, ripalo, riploa, riplolo
  use mod_gkvector, only: ngk, ngkmax, gkc, vgkc, igkig
  use mod_kpoint, only: nkpt
  use mod_eigensystem, only: nmat, nmatmax, idxlo
  use mod_atoms, only: nspecies, natoms, idxas, natmtot
  use mod_APW_LO, only: apword, apwordmax, nlorb, lorbl, nlotot, nlomax, lolmax
  use mod_muffin_tin, only: idxlm, lmmaxapw
  use modinput, only: input
  use mod_gvector, only: ivg, ivgig, cfunig
  use rttddft_GlobalVariables, only: pmat, apwalm
  use constants, only: zzero, zone, zi
  use precision, only: dp

  implicit none

  private :: genpmatbasisik
  public :: Obtain_Pmat_LAPWLOBasis
contains

!> Here, we calculat the momentum matrix elements considering as basis (L)APW+lo
!> We copied most of the code from genpmatxs.f90, but there the basis are
!> the KS-wavefunctions
subroutine Obtain_Pmat_LAPWLOBasis
  implicit none

  integer                   :: ik

  pmat(:,:,:,:) = zzero

  if(allocated(ripaa)) deallocate(ripaa)
  allocate(ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw,natmtot, 3))
  if(nlotot .gt. 0) then
    if(allocated(ripalo)) deallocate(ripalo)
    allocate(ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, natmtot, 3))
    if(allocated(riploa)) deallocate(riploa)
    allocate(riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, natmtot, 3))
    if(allocated(riplolo)) deallocate(riplolo)
    allocate(riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, natmtot, 3))
  end if

  ! Calculate gradient of radial functions times spherical harmonics
  call pmatrad

#ifdef MPI
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik) SHARED(rank,apwalm,pmat,ripaa,ripalo,riploa,riplolo)
!$OMP DO
#endif
  do ik = firstk(rank), lastk(rank)
#else
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik) SHARED(nkpt,apwalm,pmat,ripaa,ripalo,riploa,riplolo)
!$OMP DO
#endif
  do ik = 1, nkpt
#endif
    call genpmatbasisik(ik,apwalm(:,:,:,:,ik),pmat(:,:,:,ik))
  end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

  deallocate(ripaa)
  if(nlotot .gt. 0) then
    deallocate(ripalo)
    deallocate(riploa)
    deallocate(riplolo)
  end if

end subroutine Obtain_Pmat_LAPWLOBasis

  !> Here, we obtain the momentum matrix elements for a given a k-point
  !> @param[in]   ik        Index of k-point
  !> @param[in]   apwalmk   The matching coefficients for the (L)APW's
  !>                        (coefficient that smoothly connect LAPW's and
  !>                        their MT counterparts)
  !>                        Dimensions: ngkmax,apwordmax,lmmaxapw,natmtot
  !> @param[out]  pmatk     Momentum matrix elements
  !>                        Dimensions: nmatmax, nmatmax, 3
subroutine genpmatbasisik( ik, apwalmk, pmatk )
  implicit none

  !> ik        Index of k-point
  integer, intent(in)       :: ik
  !> apwalmk:  The matching coefficients for the (L)APW's (coefficient that
  !>           smoothly connect LAPW's and their MT counterparts)
  !>           Dimensions: ngkmax,apwordmax,lmmaxapw,natmtot
  complex(dp), intent(in)   :: apwalmk(:, :, :, :)
  !> pmatk     Momentum matrix elements. Dimensions: nmatmax, nmatmax, 3
  complex(dp), intent(out)  :: pmatk(:, :, :)

  integer                   :: is, ia, ias, io, io1, io2, ilo, ilo1, ilo2
  integer                   :: ig1, igp1, ig2, igp2, ig, iv1 (3), iv (3)
  integer                   :: i, j, l1, m1, lm1, l2, m2, lm2
  integer                   :: nmatp, ng
  complex(dp), allocatable  :: zv (:), aux(:,:)


  nmatp = nmat(1,ik)
  ng = ngk(1,ik)

  allocate (zv(ng))
  allocate (aux(ng,ng))

  ! loop over species and atoms
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      !---------------------------!
      !     APW-APW contribution  !
      !---------------------------!
      Do j = 1, 3 ! (loop over the x,y,z components)
        aux(:,:) = zzero
        Do l1 = 0, input%groundstate%lmaxapw
          Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
              zv (:) = zzero
              Do l2 = 0, input%groundstate%lmaxapw
                Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do io2 = 1, apword (l2, is)
                    zv(1:ng) = zv(1:ng) + ripaa(io1,lm1,io2,lm2,ias,j)*apwalmk(1:ng,io2,lm2,ias)
                  End Do
                End Do
              End Do
              Call zoutpr(ng,ng,zone,apwalmk(1:ng,io1,lm1,ias),zv(1:ng),aux(1:ng,1:ng))
            End Do
          End Do
        End Do
        pmatk(1:ng,1:ng,j) = pmatk(1:ng,1:ng,j) + aux(1:ng,1:ng)
      End Do !Do j = 1, 3
      if (nlotot .gt. 0) then
      !--------------------------------------!
      !     APW-local-orbital contribution   !
      !--------------------------------------!
        Do j = 1, 3
          Do l1 = 0,input%groundstate%lmaxapw ! (loop over lapw)
            Do m1 = - l1, l1
              lm1 = idxlm(l1, m1)
              Do io = 1, apword(l1, is)
                Do ilo = 1, nlorb(is) ! (loop over local orbitals)
                  l2 = lorbl(ilo, is)
                  Do m2 = - l2, l2
                    lm2 = idxlm(l2, m2)
                    zv(1:ng) = ripalo(io,lm1,ilo,m2,ias,j)*conjg(apwalmk(1:ng,io,lm1,ias))
                    pmatk(1:ng,ng+idxlo(lm2,ilo,ias),j) = pmatk(1:ng,ng+idxlo(lm2,ilo,ias),j) + &
                      & zv(1:ng)
                  End Do ! Do m2 = - l2, l2
                End Do ! Do ilo = 1, nlorb(is)
              End Do ! Do io = 1, apword(l1, is)
            End Do ! Do m1 = - l1, l1
          End Do ! Do l1 = 0,input%groundstate%lmaxapw
        End Do ! Do j = 1, 3
        !--------------------------------------!
        !     local-orbital-APW contribution   !
        !--------------------------------------!
        Do j = 1, 3
          Do ilo = 1, nlorb (is)
            l1 = lorbl (ilo, is)
            Do m1 = - l1, l1
              lm1 = idxlm (l1, m1)
              Do l2 = 0, input%groundstate%lmaxapw
                Do m2 = - l2, l2
                  lm2 = idxlm(l2, m2)
                  Do io = 1, apword(l2, is)
                    zv(1:ng) = riploa(ilo,m1,io,lm2,ias,j)*apwalmk(1:ng,io,lm2,ias)
                    pmatk(ng+idxlo(lm1,ilo,ias),1:ng,j) = &
                      & pmatk(ng+idxlo(lm1,ilo,ias),1:ng,j) + zv(1:ng)
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
        !------------------------------------------------!
        !     local-orbital-local-orbital contribution   !
        !------------------------------------------------!
        Do j = 1, 3
          Do ilo1 = 1, nlorb (is)
            l1 = lorbl (ilo1, is)
            Do m1 = - l1, l1
                lm1 = idxlm (l1, m1)
                Do ilo2 = 1, nlorb (is)
                   l2 = lorbl (ilo2, is)
                   Do m2 = - l2, l2
                      lm2 = idxlm (l2, m2)
                      pmatk(ng+idxlo(lm1,ilo1,ias),ng+idxlo(lm2,ilo2,ias),j) = &
                        & riplolo(ilo1,m1,ilo2,m2,ias,j)
                   End Do
                End Do
             End Do
          End Do
        End Do
        ! end case of local orbitals
      End If
    ! end loop over atoms and species
    End Do
  End Do

  pmatk(:,:,1) = -zi*pmatk(:,:,1)
  pmatk(:,:,3) = -zi*pmatk(:,:,3)

  !  calculate momentum matrix elements in the interstitial region
  do j = 1, 3
    do igp1 = 1, ng
       ig1 = igkig(igp1,1,ik)
       iv1 (:) = ivg (:, ig1)
       do igp2 = 1, ng
          ig2 = igkig(igp2,1,ik)
          iv (:) = iv1 (:) - ivg (:, ig2)
          ig = ivgig (iv(1), iv(2), iv(3))
          pmatk(igp1,igp2,j) = pmatk(igp1,igp2,j) + vgkc(j,igp2,1,ik)*cfunig(ig)
       end do
    end do
  end do

end subroutine genpmatbasisik

end module rttddft_pmat
