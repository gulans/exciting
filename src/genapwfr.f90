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
Use mod_hybrids, only: ex_coef

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


!if((ex_coef.gt.0d0).and.input%groundstate%Hybrid%updateRadial)then
!open (2, file = 'base_apw_hf.dat', status = 'old')
! Do is = 1, nspecies
!         Do ia = 1, natoms (is)
!read(2,*)
!read(2,*)!is ia
!read(2,*)
!            ias = idxas (ia, is)
!            Do l = 0, input%groundstate%lmaxapw
!               Do io1 = 1, apword (l, is)
!read(2,*)i1,i2,e1
!write(*,*)i1,i2,e1
!        apwe(io1, l, ias)=e1
!
!               enddo!io1
!            enddo!l
!         enddo!ia
!       enddo!is
!close(2)
!endif


write(*,*)"updating apw base, ex=",ex_coef, apwe(1, 0,1) 
open (12, file = 'base_apw.dat', status = 'replace')


      Do is = 1, nspecies
         nr = nrmt (is)
         Do ia = 1, natoms (is)
write(12,*)"is ia"
write(12,*)is, ia
write(12,*)"l order energy"
            ias = idxas (ia, is)
            vr (1:nr) = veffmt (1, 1:nr, ias) * y00
            Do l = 0, input%groundstate%lmaxapw
               Do io1 = 1, apword (l, is)
               if (apwe(io1, l, ias).gt.1e6) cycle
! integrate the radial Schrodinger equation
write(12,*)l,apwdm(io1, l, is),apwe(io1, l, ias)

!if Hybrid is not associated then updateRadial is not defined, but it does not give error:
if(associated(input%groundstate%Hybrid).and.input%groundstate%Hybrid%updateRadial.and.(ex_coef.ne.0d0)) then
                  Call rschroddme2 (is,ia,apwdm(io1, l, is), l, 0, apwe(io1, &
                 & l, ias), nr, spr(:, is), &
                 & vr, nn, p0(:, io1), p1(:, io1), q0(:, io1), q1(:, io1))
 else
                  Call rschroddme (apwdm(io1, l, is), l, 0, apwe(io1, &
                 & l, ias), nr, spr(:, is), &
                 & vr, nn, p0(:, io1), p1(:, io1), q0(:, io1), q1(:, io1))
 endif

 write(*,*)"rfun l=",l,"e=",apwe(io1,l, ias),"dot=",apwdm(io1, l, is),"nodes=",nn

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

                  
if(ex_coef.gt.0d0)then
WRITE(filename, '(a2,F5.2,a2,i1,a2,i1,a6)')'rf', apwe(io1,l,ias),"-o",apwdm(io1, l, is),"-l",l,'HF.dat'

else
WRITE(filename, '(a2,F5.2,a2,i1,a2,i1,a4)')'rf', apwe(io1,l,ias),"-o",apwdm(io1, l, is),"-l",l,'.dat'
endif
                  

open (11, file = filename, status = 'replace')
Do ir = 1, nr

write(11,*)spr(ir, is),",",p0(ir, io1)
enddo
close(11)


                  
                  
                  
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
                     apwfr_old (ir, 1, io1, l, ias)=apwfr (ir, 1, io1, l, ias)
                     apwfr_old (ir, 2, io1, l, ias)=apwfr (ir, 2, io1, l, ias)

                     apwfr_new (ir, 1, io1, l, ias) = t1 * p0 (ir, io1)
                     apwfr_new (ir, 2, io1, l, ias) = (p1(ir,io1)-p0(ir, io1)*t1) * t1
 
                  End Do
               End Do
            End Do
         End Do
      End Do
!      write(*,*)
close(12)



 if (.not.(associated(input%groundstate%Hybrid).and.input%groundstate%Hybrid%updateRadial.and.(ex_coef.ne.0d0))) then
    apwfr =apwfr_new
    write(*,*)"*** atjaunojam apw"
else
    write(*,*)"*** neatjaunojam apw" 
 endif


Return
End Subroutine
!EOC
