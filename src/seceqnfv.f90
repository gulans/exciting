!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
!
!
Subroutine seceqnfv(ispn, ik, nmatp, ngp, igpig, vgpc, apwalm, evalfv, evecfv)
  ! !USES:
      Use modinput
      Use modmain
      Use modfvsystem
      Use mod_hartreefock, only: ihyb
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
  !   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
  !   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
  !   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  !   Solves the secular equation,
  !   $$ (H-\epsilon O)b=0, $$
  !   for the all the first-variational states of the input $k$-point.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !EOP
  !BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: ispn
      Integer, Intent (In) :: ik
      Integer, Intent (In) :: nmatp
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
  ! local variables
      Type (evsystem) :: system
      Logical :: packed

  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!

      packed = input%groundstate%solver%packedmatrixstorage

      Call newsystem(system,packed,nmatp)
      Call hamiltonandoverlapsetup(system,ngp,apwalm,igpig,vgpc)

  !------------------------------------------------------------------------!
  !     If Hybrid potential is used apply the non-local exchange potential !
  !------------------------------------------------------------------------!
      If (associated(input%groundstate%Hybrid)) Then
         if (input%groundstate%Hybrid%exchangetypenumber == 1) then
            if (ihyb > 1) call add_vxnl(system,ik,nmatp)
         end if
      End If

  !------------------------------------!
  !     solve the secular equation     !
  !------------------------------------!
     Call solvewithlapack(system,nstfv,evecfv,evalfv)

End Subroutine seceqnfv
!EOC
