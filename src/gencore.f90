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
      use mod_corestate, only: rhocr, rwfcr, evalcr, engy_exnl_core
      !Use modmain
use modmain, only: mt_dm
Use mod_hybrids, only: ex_coef
      
      use modinteg
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
! automatic arrays
      Logical :: done (natmmax)
      Real (8) :: vr (spnrmax)


!new varieables
        integer ::countl0,ish,lmax, l, il, k,i, number_of_states
        real(8), allocatable :: wf(:,:,:),vx_wf(:,:,:),wf0(:,:,:),eig(:),occ(:,:)
        integer, allocatable :: count_l(:),number_n(:),number_l(:),nn,l_n
logical ::file_exists
!new dummy or extra variables
real(8) :: ftemp1(spnrmax),e_kin,e1,e2,hybx_coef

hybx_coef = ex_coef !mod_hybrids variable
e_kin=0d0
engy_exnl_core=0d0




      dirac_eq=(input%groundstate%CoreRelativity.eq."dirac")
     
      Do is = 1, nspecies
         done (:) = .False.
         Do ia = 1, natoms (is)
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

               
if (.false.) then               


 Do ist = 1, spnst (is)

     !call integ_v(spnr(is),is,rwfcr(:,1,ist,ias)**2,t1,atom_integw)
     !write(*,*)"norm before",t1
 If (spcore(ist, is)) Then
!    write(*,*)"eval: ",evalcr(ist, ias)
 endif
 enddo



!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ir,t1)
!$OMP DO
               Do ist = 1, spnst (is)
                  If (spcore(ist, is)) Then
! solve the Dirac equation
                     Call rdirac (0,spn(ist, is), spl(ist, is), spk(ist, &
                    & is), spnr(is), spr(:, is), vr, evalcr(ist, ias), rwfcr(:, 1, ist, ias), &
                    & rwfcr(:, 2, ist, ias),dirac_eq,.false.)
                     t1 = spocc (ist, is)
!$OMP CRITICAL
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
!$OMP END CRITICAL
                  End If
               End Do
!$OMP END DO
!$OMP END PARALLEL




else !new solver

        if (dirac_eq) then
               write(*,*) "gencore.f90 can't do CoreRelativity='dirac'"
               stop
        endif

!all orbitals with all quantum numbers and eigenvalues         
!write(*,*)"spnst",spnst(is)
!do k=1, spnst(is)
!  write(*,*)"spl",k,". l=",spl(k,is)," n=",spn(k,is)," kappa=",spk(k,is)," core=",spcore(k,is),"e=",evalcr(k,ias)
!enddo
!


!find the largest l between core orbitals
   lmax=0
   do i = 1,spnst(is)
      if ((spl(i,is).ge.lmax).and.(spcore(i,is))) then
         lmax = spl(i,is)
      end if
   end do
!write(*,*)"lmax",lmax

   allocate(count_l(lmax+1))


!counts how many diferent n for equal l
   do il = 1, lmax+1
      countl0 = 0
      nn=il
      do k = 1, spnst(is)
         if ( ((il-1).eq.spl(k,is)) .and. (spcore(k,is)) .and. (spn(k,is).eq.nn ) )  then 
 !            write(*,*)"OK: l=",il-1," n=",nn, " k=",k,"eig=",evalcr(k,ias)
             countl0=countl0+1
             nn=nn+1            
         endif
      end do
      count_l(il) = countl0
   end do
   number_of_states = sum(count_l)

!   write(*,*)"count_l", count_l
!   write(*,*)"number_of_states", number_of_states



   allocate( wf0(spnr(is),number_of_states,1), eig(number_of_states) )
   allocate( wf(spnr(is),number_of_states,1))
   allocate( vx_wf(spnr(is),number_of_states,1))
   allocate( number_n(number_of_states),number_l(number_of_states))
    allocate( occ(number_of_states,1) )

!   write(*,*)"count_l", count_l
!   write(*,*)"number_of_states", number_of_states



!stores the wave fuctions in array wf0
   l_n=1
   do il = 1, lmax+1 
      nn=il
      do k = 1, spnst(is)
         if ( ((il-1).eq.spl(k,is)) .and. (spcore(k,is)) .and. (spn(k,is).eq.nn ) )  then
!             write(*,*)"OK: l=",il-1," n=",nn, " k=",k,"eig=",evalcr(k,ias)
             wf0(:,l_n,1)=rwfcr(:,1,k,ias)/spr(:,is)
             eig(l_n)=evalcr(k,ias)
             number_n(l_n)=spn(k,is)
             number_l(l_n)=il-1
             nn=nn+1
             l_n=l_n+1
         endif
      end do
   end do

!sum occ for each l,n   
write(*,*)"Final:"
  l_n=1
  do il = 1, lmax+1 !count occ
     do nn=1,count_l(il)
        occ(l_n,1)=0d0
        do k = 1, spnst(is)
           if ( ((il-1).eq.spl(k,is)) .and. (spcore(k,is)) .and. (number_n(l_n).eq. spn(k,is)) )  then
             occ(l_n,1)=occ(l_n,1)+spocc (k, is)   
           endif
        end do
      write(*,*)"l_n=",l_n,"l=",number_l(l_n)," n=",number_n(l_n), "eig=",eig(l_n),"occ=",occ(l_n,1)
      l_n=l_n+1
      end do
   end do




!open(11,file='wf.out',status='replace')
!  do ir=1,spnr(is)
!  write(11,*) spr(ir,is),wf0(ir,6,1)
!  enddo
!  close(11)


!
!  do l_n=1, number_of_states
!     call integ_v(spnr(is),is,spr(:,is)**2*wf0(:,l_n,1)**2,t1,atom_integw)
!     write(*,*)"norm",t1
!  enddo





!-------start of iteration process, iterating through each value of l
wf=wf0
l_n=0

!write(*,*)"gencore.f90"
!read(*,*)

Do il = 1, lmax+1
       call LS_iteration(spnr(is),is,ia,hybx_coef,spr(:,is),vr,il-1,number_l,occ(:,1),1,count_l(il),l_n,& !in
            number_of_states,1,.false.,lmax, wf0,& !in
            wf(:,l_n+1:l_n+count_l(il),1),vx_wf(:,l_n+1:l_n+count_l(il),1),&
            eig(l_n+1:l_n+count_l(il)))  

!    write(*,*)"eig",eig(l_n+1:l_n+count_l(il))

    l_n=l_n+count_l(il)
enddo
       ! call LS_iteration(spnr(is),is, spr(:,is),tools, tools_info,il-1,spin,&
       !                                 count_l(il),mspin,vr,&
       !                                 wf(:,:,spin),eig(:))





!!----------engy_exnl_core
! integral <psi|vx_nl|psi>

if(hybx_coef.gt.1e-10)then

l_n=0
Do il = 1, lmax+1
  do nn=1,count_l(il)
  l_n=l_n+1
  call integ_v(spnr(is),is,wf(:,l_n,1)*vx_wf(:,l_n,1)*spr(:,is)**2,t1,atom_integw)
write(*,*)"l_n=",l_n,"Ex_nl_core=",occ(l_n, 1)*t1/2d0
engy_exnl_core = engy_exnl_core + occ(l_n, 1)*t1

      
   enddo
enddo
!read(*,*)


!----------engy_exnl_core

endif




!!!---------put the corresponding wf in correct places in rwfcr

l_n=0

Do il = 1, lmax+1
  do nn=1,count_l(il)
    l_n=l_n+1
    do ist=1,spnst(is)
      if ( (spcore(ist, is)).and.(spl(ist,is).eq.(il-1)).and.(spn(ist,is).eq.number_n(l_n)) ) then
              write(*,*)evalcr(ist,ias),"<-",eig(l_n)
              evalcr(ist,ias)=eig(l_n) 
write(*,*)shape(spr(:,is)),shape(wf(:,l_n,1))
              rwfcr(:,1,ist,ias)=wf(:,l_n,1)*spr(:,is)
      endif
    enddo 
  enddo
enddo

!Kinetic energy

do ist = 1, spnst (is)
if (spcore(ist, is)) Then

  call deriv_f(spnr(is),is,rwfcr(:,1,ist,ias)/spr(:,is),ftemp1,atom_integw)
  call integ_v(spnr(is),is,0.5d0*ftemp1**2*spr(:,is)**2,e1,atom_integw)
  call integ_v(spnr(is),is,0.5d0*dble(spl(ist,is))*dble(spl(ist,is)+1)*&
          (rwfcr(:,1,ist,ias)/spr(:,is))**2,e2,atom_integw)
  e_kin=e_kin+spocc(ist, is)*(e1+e2)
write(*,*)"kinetic new:",e_kin
endif
enddo

!read(*,*)















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

  deallocate(count_l,wf,vx_wf,eig)

   deallocate( wf0 )
   deallocate( number_n)
   deallocate( number_l)

   deallocate( occ )
        
        
!___________________________
!end of LS iteration process










endif !new solver solved





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
                     End Do !ist
                     rhocr (1:spnr(is), jas) = rhocr (1:spnr(is), ias)
                     done (ja) = .True.
                  End If
               End Do !ja
            End If
         End Do!ia
      End Do!is

      Return
End Subroutine
!EOC
