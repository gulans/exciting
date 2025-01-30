!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genapwfr
! !INTERFACE:
!
!
Subroutine genapwfr
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates the APW radial functions. This is done by integrating the scalar
!   relativistic Schr\"{o}dinger equation (or its energy deriatives) at the
!   current linearisation energies using the spherical part of the effective
!   potential. The number of radial functions at each $l$-value is given by the
!   variable {\tt apword} (at the muffin-tin boundary, the APW functions have
!   continuous derivatives up to order ${\tt apword}-1$). Within each $l$, these
!   functions are orthonormalised with the Gram-Schmidt method. The radial
!   Hamiltonian is applied to the orthonormalised functions and the results are
!   stored in the global array {\tt apwfr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, nr, ir
      Integer :: nn, l, io1, io2
      Real (8) :: t1
! automatic arrays
      Real (8) :: vr (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
      Real (8) :: p0 (nrmtmax, apwordmax), p1 (nrmtmax, apwordmax), p1s &
     & (apwordmax)
      Real (8) :: q0 (nrmtmax, apwordmax), q1 (nrmtmax, apwordmax)
      Real (8) :: hp0 (nrmtmax)
      character(len=1024) :: filename
      Logical :: done (natmmax,nspecies)
      integer :: ja,jas
      call stopwatch("exciting:genapwfr", 1)
      done (:,:) = .False.
!$OMP PARALLEL DEFAULT(none) PRIVATE(is,nr,ia,ias,vr,l,io1,nn,p0,p1,q0,q1,ir,fr,gr,cf,t1,p1s,io2,ja,jas) SHARED(nspecies,done,nrmt,natoms,idxas,veffmt,input,apword,ex_coef,apwdm,apwe,spr,apwfr_old,apwfr_new,apwfr,eqatoms)
      Do is = 1, nspecies   
         nr = nrmt (is)
         Do ia = 1, natoms (is)
            If ( .Not. done(ia,is)) Then  
            ias = idxas (ia, is)
            vr (1:nr) = veffmt (1, 1:nr, ias) * y00
!$OMP DO
            Do l = 0, input%groundstate%lmaxapw
               Do io1 = 1, apword (l, is)
! integrate the radial Schrodinger equation

               if(associated(input%groundstate%Hybrid).and.input%groundstate%Hybrid%updateRadial.and.(ex_coef.ne.0d0)) then
                  Call rschroddme2 (is,ia,apwdm(io1, l, is), l, 0, apwe(io1, &
                    & l, ias), nr, spr(:nr, is), &
                    & vr(:nr), nn, p0(:nr, io1), p1(:nr, io1), q0(:nr, io1), q1(:nr, io1))
               else

                  Call rschroddme (apwdm(io1, l, is), l, 0, apwe(io1, &
                 & l, ias), nr, spr(:, is), &
                 & vr, nn, p0(:, io1), p1(:, io1), q0(:, io1), q1(:, io1))
               endif
! normalise radial functions
                  Do ir = 1, nr
                     fr (ir) = p0 (ir, io1) ** 2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  t1 = 1.d0 / Sqrt (Abs(gr(nr)))
                  p1s (io1) = t1 * p1 (nr, io1)
                  p0 (1:nr, io1) = t1 * p0 (1:nr, io1)
                  p1 (1:nr, io1) = t1 * p1 (1:nr, io1)
                  q0 (1:nr, io1) = t1 * q0 (1:nr, io1)
                  q1 (1:nr, io1) = t1 * q1 (1:nr, io1)


                  ! if(ex_coef.gt.0d0)then
                  !    WRITE(filename, '(a2,F5.2,a2,i1,a2,i1,a6)')'rf', apwe(io1,l,ias),"-o",apwdm(io1, l, is),"-l",l,'HF.dat'   
                  ! else
                  !    WRITE(filename, '(a2,F5.2,a2,i1,a2,i1,a4)')'rf', apwe(io1,l,ias),"-o",apwdm(io1, l, is),"-l",l,'.dat'
                  ! endif
                                       
                  ! open (11, file = filename, status = 'replace')
                  ! Do ir = 1, nr
                  !    write(11,*)spr(ir, is),",",p0(ir, io1)
                  ! enddo
                  ! close(11)


! subtract linear combination of previous vectors
                  Do io2 = 1, io1 - 1
                     Do ir = 1, nr
                        fr (ir) = p0 (ir, io1) * p0 (ir, io2)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     t1 = gr (nr)
                     p1s (io1) = p1s (io1) - t1 * p1s (io2)
                     p0 (1:nr, io1) = p0 (1:nr, io1) - t1 * p0 (1:nr,io2)
                     p1 (1:nr, io1) = p1 (1:nr, io1) - t1 * p1 (1:nr,io2)
                     q0 (1:nr, io1) = q0 (1:nr, io1) - t1 * q0 (1:nr,io2)
                     q1 (1:nr, io1) = q1 (1:nr, io1) - t1 * q1 (1:nr,io2)
                  End Do
! normalise radial functions
                  Do ir = 1, nr
                     fr (ir) = p0 (ir, io1) ** 2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  t1 = Abs (gr(nr))
                  If (t1 .Lt. 1.d-20) Then
                     Write (*,*)
                     Write (*, '("Error(genapwfr): degenerate APW radia&
                    &l functions")')
                     Write (*, '(" for species ", I4)') is
                     Write (*, '(" atom ", I4)') ia
                     Write (*, '(" angular momentum ", I4)') l
                     Write (*, '(" and order ", I4)') io1
                     Write (*,*)
                     Stop
                  End If
                  t1 = 1.d0 / Sqrt (t1)
                  p1s (io1) = t1 * p1s (io1)
                  p0 (1:nr, io1) = t1 * p0 (1:nr, io1)
                  p1 (1:nr, io1) = t1 * p1 (1:nr, io1)
                  q0 (1:nr, io1) = t1 * q0 (1:nr, io1)
                  q1 (1:nr, io1) = t1 * q1 (1:nr, io1)
                  Do ir = 1, nr
                     t1 = 1.d0 / spr (ir, is)
! !$OMP CRITICAL                     
                     apwfr_old (ir, 1, io1, l, ias)=apwfr (ir, 1, io1, l, ias)
                     apwfr_old (ir, 2, io1, l, ias)=apwfr (ir, 2, io1, l, ias)

                     apwfr_new (ir, 1, io1, l, ias) = t1 * p0 (ir, io1)
                     apwfr_new (ir, 2, io1, l, ias) = (p1(ir,io1)-p0(ir, io1)*t1) * t1
! !$OMP END CRITICAL                      
                  End Do
               End Do! io1
            End Do! l
!$OMP END DO
            done (ia,is) = .True.
! copy to equivalent atoms
            Do ja = 1, natoms (is)
               jas = idxas (ja, is)
               If (( .Not. done(ja,is)) .And. (eqatoms(ia, ja, is))) Then
                  apwfr_old (:, :, :, :, jas) = apwfr_old (:, :, :, :, ias)
                  apwfr_new (:, :, :, :, jas) = apwfr_new (:, :, :, :, ias)
                  done (ja,is) = .True.
               End If
            End Do !ja
         endif! if done
         End Do! ia
      End Do! is
!$OMP END PARALLEL
      if(.not.(associated(input%groundstate%Hybrid).and.input%groundstate%Hybrid%updateRadial.and.(ex_coef.ne.0d0))) then
         apwfr =apwfr_new
         !write(*,*)"*** atjaunojam apw"
      else
         !write(*,*)"*** neatjaunojam apw" 
      endif

      
      
      call stopwatch("exciting:genapwfr", 0)
      Return
End Subroutine
!EOC
