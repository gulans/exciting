!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genlofr
! !INTERFACE:
!
!
Subroutine genlofr
! !USES:
      use modmpi, only: mpiglobal
      Use modinput
      Use modmain
! !DESCRIPTION:
!> Generates the local-orbital radial functions. This is done by integrating
!> the scalar relativistic Schr\"{o}dinger equation or Dirac equation
!> (or its energy deriatives) at the current linearisation energies using the
!> spherical part of the effective potential. Dirac-type local orbitals are
!> useful in the context of spin-orbit coupling, but should only be used along
!> with second variation with local orbitals. For more details see:
!> arXiv:2306.02965 [cond-mat.mtrl-sci].
!> For each local-orbital, a linear combination of {\tt lorbord} radial functions is
!> constructed such that its radial derivatives up to order ${\tt lorbord}-1$ are zero
!> at the muffin-tin radius.This function is normalised and the radial Hamiltonian applied
!> to it. The results are stored in the global array {\tt lofr}.
!
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, nr, ir
      Integer :: ilo, io1, io2
      Integer :: j, l, nn, info, np
      Real (8) :: t1
! automatic arrays
      Real (8) :: vr (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
      Real (8) :: p0 (nrmtmax, maxlorbord), p1 (nrmtmax,maxlorbord)
      Real (8) :: q0 (nrmtmax, maxlorbord), q1 (nrmtmax, maxlorbord)
      Real (8) :: p0s (nrmtmax), p1s (nrmtmax) ,q0s (nrmtmax), q1s (nrmtmax)
      Real (8) :: hp0 (nrmtmax)
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: xa (:), ya (:)
      Real (8), Allocatable :: a (:, :), b (:), c (:)
! external functions
      Real (8) :: polynom
      External polynom
! variables for the lo recommendation
      Real (8) energy,energyp,tmp,tmp2,ens(0:20),elo,ehi,flo,fhi,emi,fmi
      integer nodes
      character(len=1024) :: filename
      Logical :: done (natmmax,nspecies)
      integer :: ja,jas
      call stopwatch("exciting:genlofr", 1)
      
      np = Max (maxlorbord+1, 4)
      done (:,:) = .False.
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(input,ex_coef,np,nspecies,nrmt,done,natoms,idxas,veffmt,nlorb,lorbl,lorbord,wfkappa,lorbdm,lorbe,spr,lofr_old,lofr_new,lofr,eqatoms,rmt) !y00 un polynom 
      Allocate (ipiv(np))
      Allocate (xa(np), ya(np), c(np))
      Allocate (a(np, np), b(np))
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ia = 1, natoms (is)
            If ( .Not. done(ia,is)) Then  
            ias = idxas (ia, is)
            vr (1:nr) = veffmt (1, 1:nr, ias) * y00
!$OMP DO
            Do ilo = 1, nlorb (is)
               l = lorbl (ilo, is)
               Do io2 = 1, lorbord (ilo, is)
                  if (wfkappa(io2, ilo, is) /= 0) then
                     ! integrate the radial Dirac equation                                                                                                           
                     Call rdiracdme (lorbdm(io2, ilo, is), wfkappa(io2, ilo, is), &
                          & lorbe(io2, ilo, ias), nr, &
                          & spr(:, is), vr, nn, p0(:, io2), p1(:, io2), q0(:, io2), &
                          & q1(:, io2), .false.)
                  else
                     ! integrate the radial Schrodinger equation
                     if(associated(input%groundstate%Hybrid).and.input%groundstate%Hybrid%updateRadial.and.(ex_coef.ne.0d0)) then
                        Call rschroddme2 (is,ia,lorbdm(io2, ilo, is), l, 0, &
                                 & lorbe(io2, ilo, ias), nr, &
                                 & spr(:nr, is), vr(:nr), nn, p0(:nr, io2), p1(:nr, io2), q0(:nr, io2), &
                                 & q1(:nr, io2))
                     else
                        Call rschroddme (lorbdm(io2, ilo, is), l, 0, &
                           & lorbe(io2, ilo, ias), nr, &
                           & spr(:, is), vr, nn, p0(:, io2), p1(:, io2), q0(:, io2), &
                           & q1(:, io2))
                     endif
                  endif
! normalise radial functions
                  Do ir = 1, nr
                     fr (ir) = p0 (ir, io2) ** 2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  t1 = 1.d0 / Sqrt (Abs(gr(nr)))
                  p0 (1:nr, io2) = t1 * p0 (1:nr, io2)
                  p1 (1:nr, io2) = t1 * p1 (1:nr, io2)
                  q0 (1:nr, io2) = t1 * q0 (1:nr, io2)
                  q1 (1:nr, io2) = t1 * q1 (1:nr, io2)

               if (lorbdm(io2, ilo, is).eq.0) then
!$OMP CRITICAL
               if(.false.) then

                  if(ex_coef.gt.0d0)then
                     WRITE(filename, '(a3,A2,a2,i1,F6.2,a6)')'rf-',input%structure%speciesarray(is)%species%chemicalSymbol,"-l",l,lorbe(io2, ilo, ias),'HF.dat'
                  else
                     WRITE(filename, '(a3,A2,a2,i1,F6.2,a7)')'rf-',input%structure%speciesarray(is)%species%chemicalSymbol,"-l",l,lorbe(io2, ilo, ias),'PBE.dat'
                  endif
                  open (11, file = filename, status = 'replace')
                  Do ir = 1, nr
                  write(11,*)spr(ir, is),",",p0(ir, io2)
                  enddo
                  close(11)
               endif
!$OMP END CRITICAL
          endif

! set up the matrix of radial derivatives
                  Do j = 1, np
                     ir = nr - np + j
                     xa (j) = spr (ir, is)
                     ya (j) = p0 (ir, io2) / spr (ir, is)
                  End Do
                  Do io1 = 1, lorbord (ilo, is)
                     a (io1, io2) = polynom (io1-1, np, xa, ya, c, &
                    & rmt(is))
                  End Do
               End Do !io2
! set up the target vector
               b (:) = 0.d0
               b (lorbord(ilo, is)) = 1.d0
               Call dgesv (lorbord(ilo, is), 1, a, np, ipiv, b, np, &
              & info)
               If (info .Ne. 0) Then
                  Write (*,*)
                  Write (*, '("Error(genlofr): degenerate local-orbital&
                 & radial functions")')
                  Write (*, '(" for species ", I4)') is
                  Write (*, '(" atom ", I4)') ia
                  Write (*, '(" and local-orbital ", I4)') ilo
                  Write (*, '(" ZGESV returned INFO = ", I8)') info
                  Write (*,*)
                  Stop
               End If
! generate linear superposition of radial functions
               p0s (:) = 0.d0
               p1s (:) = 0.d0
               q0s (:) = 0.d0
               q1s (:) = 0.d0
               Do io1 = 1, lorbord (ilo, is)
                 t1 = b (io1)
                 p0s (1:nr) = p0s (1:nr) + t1 * p0 (1:nr, io1)
                 p1s (1:nr) = p1s (1:nr) + t1 * p1 (1:nr, io1)
                 q0s (1:nr) = q0s (1:nr) + t1 * q0 (1:nr, io1)
                 q1s (1:nr) = q1s (1:nr) + t1 * q1 (1:nr, io1)
               End Do
! normalise radial functions
               Do ir = 1, nr
                 fr (ir) = p0s (ir) ** 2
               End Do
               Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
               t1 = 1.d0 / Sqrt (Abs(gr(nr)))
               p0s (1:nr) = t1 * p0s (1:nr)
               p1s (1:nr) = t1 * p1s (1:nr)
               q0s (1:nr) = t1 * q0s (1:nr)
               q1s (1:nr) = t1 * q1s (1:nr)
               Do ir = 1, nr
                 t1 = 1.d0 / spr (ir, is)
                 lofr_old (ir, 1, ilo, ias) = lofr (ir, 1, ilo, ias)
                 lofr_old (ir, 2, ilo, ias) = lofr (ir, 2, ilo, ias) 
                 lofr_new (ir, 1, ilo, ias) = t1 * p0s (ir)
                 lofr_new (ir, 2, ilo, ias) = (p1s(ir)-p0s(ir)*t1) * t1
               End Do
            End Do !ilo
!$OMP END DO
            done (ia,is) = .True.
! copy to equivalent atoms
            Do ja = 1, natoms (is)
               jas = idxas (ja, is)
               If (( .Not. done(ja,is)) .And. (eqatoms(ia, ja, is))) Then
                  !write(*,*)"genlofr ekvivalenti atomi"
                  lofr_old (:, :, :, jas) = lofr_old (:, :, :, ias)
                  lofr_new (:, :, :, jas) = lofr_new (:, :, :, ias)
                  done (ja,is) = .True.
               End If
            End Do !ja
         endif! if done

         End Do! ia
      End Do! is
      Deallocate (ipiv, xa, ya, a, b, c)
     
!$OMP END PARALLEL
      if(.not.(associated(input%groundstate%Hybrid).and.input%groundstate%Hybrid%updateRadial.and.(ex_coef.ne.0d0))) then
         lofr =lofr_new
         !write(*,*)"*** atjaunojam lo"
      else
         !write(*,*)"*** neatjaunojam lo" 
      endif



      call stopwatch("exciting:genlofr", 0)
      Return
End Subroutine
!EOC
