!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gencore
! !INTERFACE:
!
!
Subroutine gencore
! !USES:
      use modinput, only: input
      use mod_atoms, only: natoms, idxas, spvr, spnr, nspecies, &
        & spnst, spcore, spn, spr, spk, spl, spocc, natmmax, &
        & spnrmax
      use mod_symmetry, only: eqatoms
      use mod_muffin_tin, only: nrmt
      use mod_potential_and_density, only: veffmt
      use constants, only: y00, fourpi
      use mod_corestate, only: rhocr, rwfcr, evalcr
      !Use modmain
! !DESCRIPTION:
!   Computes the core radial wavefunctions, eigenvalues and densities. The
!   radial Dirac equation is solved in the spherical part of the effective
!   potential to which the atomic potential has been appended for
!   $r>R^{\rm MT}$. In the case of spin-polarised calculations, the effective
!   magnetic field is ignored.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ja, ias, jas, ist, ir
      Real (8) :: t1
      logical dirac_eq
      logical ls_iter
! automatic arrays
      Logical :: done (natmmax)
      Real (8) :: vr (spnrmax)
      
!new varieables
	
      integer :: d_order,i_order,tools_info(3)
      integer ::lmax, il, spin, mspin,k,i, number_of_states,countl0
      real(8), allocatable ::tools(:,:), wf(:,:,:),eig(:)
      integer, allocatable :: count_l(:)



      dirac_eq=(input%groundstate%CoreRelativity.eq."dirac")
      ls_iter=(input%groundstate%CoreLS.eq."ls")
write(*,*) "dirac", dirac_eq, "ls", ls_iter
      Do is = 1, nspecies
	write(*,*)"nspecies", nspecies
         done (:) = .False.
         Do ia = 1, natoms (is)
		!write(*,*)"natoms", natoms(is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
               If (input%groundstate%frozencore) Then
! use atomic potential for the frozen core approximation
                  vr (1:nrmt(is)) = spvr (1:nrmt(is), is)
               Else
! else use the spherical part of the crystal effective potential
                  vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
               End If
! append the effective potential from the atomic calculation
               t1 = vr (nrmt(is)) - spvr (nrmt(is), is)
               Do ir = nrmt (is) + 1, spnr (is)
                  vr (ir) = spvr (ir, is) + t1
               End Do
               rhocr (:, ias) = 0.d0
!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ir,t1)

	if (dirac_eq) then
	write(*,*) "dirac"
!!$OMP DO
               Do ist = 1, spnst (is)
                  If (spcore(ist, is)) Then
! solve the Dirac equation



                     Call rdirac (0,spn(ist, is), spl(ist, is), spk(ist, &
                    & is), spnr(is), spr(:, is), vr, evalcr(ist, ias), rwfcr(:, 1, ist, ias), &
                    & rwfcr(:, 2, ist, ias),dirac_eq,.false.)
                     t1 = spocc (ist, is)
!!$OMP CRITICAL
                     if (dirac_eq) then
                       Do ir = 1, spnr (is)
! add to the core density
                           rhocr (ir, ias) = rhocr (ir, ias) + t1 * &
                         & (rwfcr(ir, 1, ist, ias)**2+rwfcr(ir, 2, ist, &
                         & ias)**2)
                       End Do
                     else
                       Do ir = 1, spnr (is)
                           rhocr (ir, ias) = rhocr (ir, ias) + t1 * &
                         & rwfcr(ir, 1, ist, ias)**2
                       Enddo
                     endif
!!$OMP END CRITICAL




                  End If
               End Do
!!$OMP END DO
!!$OMP END PARALLEL
	else if (ls_iter) then
	write(*,*)"LS iteration"

	d_order=9
	i_order=9
        tools_info=(/40,d_order,i_order/)
	 !info about tools_info array:
	 !tools_info(1) - (2-nd dimenstion size of tools array 1-10 - coefficinets for Lagrange interp,
	 !                11-20 coefficints for (ri+1/4) interpolation, and 21-30 for (ri+2/4), 31-40 for (ri+3/4) 
	 !tools_info(2) - order of intrpolation polynom for Lagrange interpolation when calculating derivative,
	 !tools_info(3) - order for inperolation for integration with Bodes folmula)

	Allocate(tools(spnr(is),tools_info(1)))
	call generate_tools(spnr(is),spr(:,is),tools,tools_info)

	lmax=0
	do i = 1,spnst(is)
		if ((spl(i,is).ge.lmax).and.(spcore(i,is))) then
		lmax = spl(i,is)
		end if
	end do
	allocate(count_l(lmax+1))



		do ist = 1, lmax+1!counts how many shells of equal l
			countl0 = 0
			do il = 1, size(spl(:,is))
				if (((ist-1).eq.spl(il,is)).and.(spcore(il,is))) then	
					countl0=countl0+1
					if ((spn(il,is).eq.spn(il+1,is)).and.(spl(il,is).eq.(spl(il+1,is)))) then
					!takes care of the equal n,l but different kappa case
						countl0 = countl0-1
					end if
				end if
			end do
			count_l(ist) = countl0
		end do!end shell counting

		allocate(wf(spnr(is),maxval(count_l(:)),1), eig(maxval(count_l(:))))



		number_of_states = sum(count_l)
!-------start of iteration process, iterating through each value of l
		Do il = 1, lmax+1
			k=1
!!!------get wavefunctions of previous iterations
			do i = 1,spnst(is)!putting rwfcr in correct positions for L
				if ((spl(i,is).eq.(il-1)).and.(spcore(i,is))) then
					wf(:,k,1) = rwfcr(:,1,i,ias)/spr(:,is)
					eig(k) = evalcr(i,ias)
					k=k+1
					if (spn(i,is).eq.spn(i+1,is)) then
						cycle
					end if
					if (k.eq.number_of_states) then
						exit
					end if
				end if
			end do
!!!-------
			

			mspin= 1!right now spin not taken into account
			Do spin = 1,mspin
				call LS_iteration(spnr(is), spr(:,is),tools, tools_info,il-1,spin,&
					count_l(il),mspin,vr,&
					wf(:,:,spin),eig(:))

			end Do

			
!!!---------put the corresponding wf in correct places in rwfcr

			k = 1
			do ist = 1,spnst(is)

				If ((spcore(ist, is)).and.(spl(ist,is).eq.(il-1))) Then
					write(*,*)"l", spl(ist,is), "n",spn(ist,is),"k", ist
					rwfcr(:,1,ist,ias)=wf(:,k,mspin)*spr(:,is)
					evalcr(ist,ias)=eig(k)
				

					if ((spn(ist,is).eq.spn(ist+1,is)).and.(spl(ist,is).eq.spl(ist+1,is))) then
						rwfcr(:,1,ist,ias)=wf(:,k,mspin)*spr(:,is)
						write(*,*)"n", spn(ist,is)
						cycle
					end if
					k=k+1
				end IF
			end do
!!!-----------
		
		end Do!il
!-----------end iteration process

!-------- add to core density
		do ist = 1, spnst(is)
		t1 = spocc (ist, is)
		    If (spcore(ist, is)) Then
                       Do ir = 1, spnr (is)
                           rhocr (ir, ias) = rhocr (ir, ias) + t1 * &
                         & rwfcr(ir, 1, ist, ias)**2
                       Enddo
		     end IF
		end do

	deallocate(count_l,wf,tools,eig)
!___________________________
!end of LS iteration process
	end if

               Do ir = 1, spnr (is)
                  rhocr (ir, ias) = rhocr (ir, ias) / (fourpi*spr(ir, &
                 & is)**2)
               End Do
               done (ia) = .True.
! copy to equivalent atoms
               Do ja = 1, natoms (is)
                  If (( .Not. done(ja)) .And. (eqatoms(ia, ja, is))) &
                 & Then
                     jas = idxas (ja, is)
                     Do ist = 1, spnst (is)
                        If (spcore(ist, is)) Then
                           evalcr (ist, jas) = evalcr (ist, ias)
                           rwfcr (1:spnr(is), :, ist, jas) = rwfcr &
                          & (1:spnr(is), :, ist, ias)
                        End If
                     End Do
                     rhocr (1:spnr(is), jas) = rhocr (1:spnr(is), ias)
                     done (ja) = .True.
                  End If
               End Do
			
           End If
         End Do
		
      End Do
      Return
End Subroutine
!EOC
