!
!
!
! Copyright (C) 2021 D. Zavickis, A. Gulans
!
!
Subroutine FockExchange (ikp, q0corr, vnlvv, vxpsiirgk, vxpsimt)
      Use modmain 
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Use potentials, only: coulomb_potential
      USE OMP_LIB
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Real (8), Intent (In) :: q0corr
      Complex (8), Intent (Out) :: vnlvv (nstsv, nstsv)
      Complex (8), Intent (Out) :: vxpsiirgk (ngkmax, nstsv)
      Complex (8), Intent (Out) :: vxpsimt (lmmaxvr, nrcmtmax, natmtot, nstsv)

! local variables
      Integer :: ngknr, ik, jk, ist1, ist2, ist3
      Integer :: is, ia, ias, ic, m1, m2, lmax, lm
      Integer :: nrc, iq, ig, iv (3), igq0, igk
      Integer :: ilo, loindex
      Integer :: info
      Integer :: ifg
      Logical :: solver

      Real (8) :: v (3), cfq, ta,tb, t1, norm, uir
      Complex (8) zrho01, zrho02, ztmt, zt2, ztir
      Integer :: nr, l, m, io1, lm2, ir, if3
! automatic arrays
      Real (8) :: zn (nspecies)
      Complex (8) sfacgq0 (natmtot)
! allocatable arrays
      Real (8), Allocatable :: vgqc (:, :)
      Real (8), Allocatable :: vtest (:)
      Real (8), Allocatable :: tpgqc (:, :)
      Real (8), Allocatable :: gqc (:)
      Real (8), Allocatable :: jlgqr (:, :, :)
      Real (8), Allocatable :: jlgq0r (:, :, :)
      Real (8), Allocatable :: evalfv (:,:)
      Complex (8), Allocatable :: prodir (:)
      Complex (8), Allocatable :: potir (:)
      Complex (8), Allocatable :: ylmgq (:, :)
      Complex (8), Allocatable :: sfacgq (:, :)
      Complex (8), Allocatable :: vxpsiirtmp (:)
      Complex (8), Allocatable :: wfcr1 (:, :)
      Complex (8), Allocatable :: wf1ir (:)
      Complex (8), Allocatable :: wf2ir (:)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :, :)
      Complex (8), Allocatable :: zvcltp (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: zwfir (:)
      Complex (8), Allocatable :: ztir1(:,:)
      type (WFType) :: wf1,wf2,prod,pot
! external functions
      Complex (8) zfinp, zfmtinp, zfinpir, zfinpmt
      External zfinp, zfmtinp, zfinpir, zfinpmt
      integer :: thread
! allocate local arrays
      Allocate (vtest(3))
      Allocate (vgqc(3, ngvec))
      Allocate (tpgqc(2, ngvec))
      Allocate (gqc(ngvec))
      Allocate (jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, ngvec, nspecies))
      Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
      Allocate (ylmgq(lmmaxvr, ngvec))
      Allocate (sfacgq(ngvec, natmtot))
      Allocate (vxpsiirtmp(ngrtot))
      Allocate (wfcr1(ntpll, nrcmtmax))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvcltp(ntpll, nrcmtmax))
      Allocate (zfmt(lmmaxvr, nrcmtmax))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot, nstsv))
      Allocate (zwfir(ngkmax))

write(*,*)nrcmtmax, ngkmax, "parametri"
      if (allocated(evalfv)) deallocate(evalfv)
      allocate(evalfv(nstfv,kset%nkpt))
      
      ! test comment for push
write(*,*)nrcmtmax, ngkmax, nrmtmax
      evalfv(:,:) = 0.d0
      Do ik = 1, nkpt
         Call getevalfv(kset%vkl(:,ik), evalfv(:,ik))
      End Do
      call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)

      call WFInit(wf1)

      ik  = kset%ikp2ik(ikp) ! 1d reduced index -> 1d non-reduced k-point index
      call genWF(ik,wf1)
      call genWFinMT(wf1)
      call genWFonMesh(wf1)


      call WFInit(wf2)

      t1 = 1/sqrt(omega)
! factor for long-range term
      cfq = 0.5d0 * (omega/pi) ** 2
! set the nuclear charges to zero
      zn (:) = 0.d0
      vxpsiirgk (:, :) = 0.d0
      vxpsimt (:, :, :, :) = 0.d0
      vnlvv (:, :) = 0.d0
! calculate the wavefunctions for all states for the input k-point

! if (.true.) then
! start loop over non-reduced k-point set
      Do iq = 1, kqset%nkpt

call timesec(ta)
         jk  = kqset%kqid(ik,iq) ! ID(k') = ID(k-q)-> ID(k) set

! determine q-vector
         v (:) = kqset%vkc (:, ik) - kqset%vkc (:, jk)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig) SHARED(vgqc,vgc,v,ngvec,gqc,tpgqc,ylmgq,input) 
!$OMP DO
         Do ig = 1, ngvec
! determine G+q vectors
            vgqc (:, ig) = vgc (:, ig) + v (:) ! Checked: vgc == Gset%vgc

! G+q-vector length and (theta, phi) coordinates
            Call sphcrd (vgqc(:, ig), gqc(ig), tpgqc(:, ig))
! spherical harmonics for G+q-vector
            Call genylm (input%groundstate%lmaxvr, tpgqc(:, ig), ylmgq(:, ig))
         End Do
!$OMP END DO
!$OMP END PARALLEL

! structure factors for G+q
         Call gensfacgp (ngvec, vgqc, ngvec, sfacgq)
! find the shortest G+q-vector
         Call findigp0 (ngvec, gqc, igq0)
         sfacgq0 (:) = sfacgq (igq0, :)
! compute the required spherical Bessel functions
         lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
         Call genjlgpr (lmax, gqc, jlgqr)
         Call genjlgq0r (gqc(igq0), jlgq0r)
call timesec(tb)
write(*,*) 'qpt init', tb-ta

! calculate the wavefunctions for occupied states

call timesec(ta)

         call genWF(jk,wf2)
         call genWFinMT(wf2)
         call genWFonMesh(wf2)
call timesec(tb)
write(*,*) 'genWFs',tb-ta


         solver = (input%groundstate%hybrid%singularity.ne."exc")
         
         zvclmt (:, :, :, :) = 0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist3,wf1ir,wf2ir,igk,ifg,prod,prodir,zrho01,pot,potir,vxpsiirtmp) 
         call WFInit(prod)
         call WFInit(pot)
!write(*,*)"after init wf"
         Allocate(pot%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
         Allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
         !write(*,*)"after po and prod"
         Allocate (wf1ir(ngrtot))
         Allocate (wf2ir(ngrtot))
         Allocate (prodir(ngrtot))
         Allocate (potir(ngrtot))
         pot%mtrlm(:,:,:,:) = 0.d0
         zrho01 = 0.d0
!$OMP DO
         Do ist3 = 1, nstfv

            wf1ir(:) = 0.d0
            Do igk = 1, Gkqset%ngk (1, ik)
               ifg = igfft (Gkqset%igkig(igk, 1, ik))
               wf1ir(ifg) = t1*wf1%gk(igk, ist3)
            End Do
            Call zfftifc (3, ngrid, 1, wf1ir(:))

 ! !$OMP DO
            Do ist2 = 1, nomax

               vxpsiirtmp(:) = 0.d0
               wf2ir(:) = 0.d0 
               Do igk = 1, Gkqset%ngk (1, jk)
                  ifg = igfft (Gkqset%igkig(igk, 1, jk))
                  wf2ir(ifg) = t1*wf2%gk(igk, ist2)
               End Do
               Call zfftifc (3, ngrid, 1, wf2ir(:))

   ! calculate the complex overlap density
   !-----------------------------------------------------------------------------------

               call WFprodrs(ist2,wf2,ist3,wf1,prod)
               prodir(:)=conjg(wf2ir(:))*wf1ir(:)

   if ((ik.eq.jk).and.(.not.solver)) then
                  Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                  & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prodir(:), zrho01) 
                  prodir(:)=prodir(:)-zrho01
                  prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-zrho01/y00
   endif
                  Call coulomb_potential (nrcmt, rcmt, ngvec, gqc, igq0, &
                  & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:,:,:,1), &
                  & prodir(:), pot%mtrlm(:,:,:,1), potir(:), zrho02, &
                  & cutoff=input%groundstate%hybrid%singularity.eq."exc0d")

!               call WFprodrs(ist2,wf2,ist3,wf1,prod)
!               prodir(:)=conjg(wf2ir(:))*wf1ir(:)

   if ((ik.eq.jk).and.(.not.solver)) then
                  Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                  & igq0), sfacgq0, pot%mtrlm(:,:,:,1), potir(:), zrho01)

                  potir(:)=potir(:)-zrho01
                  pot%mtrlm(1,:,:,1)=pot%mtrlm(1,:,:,1)-zrho01/y00
   endif

   !-----------------------------------------------------------------------------------
               call genWFonMeshOne(pot)
               pot%mtmesh=conjg(pot%mtmesh)
               call WFprodrs(1,pot,ist2,wf2,prod)
               vxpsiirtmp(:) = potir(:)*wf2ir(:)*wkptnr(jk)*cfunir(:)

   ! ----------------------------------------------------------------------------------
               ! Calculate Fourier transform of vxpsiirtmp(r)
               Call zfftifc (3, ngrid,-1,vxpsiirtmp)
               Do igk=1, Gkqset%ngk (1, ik)
                  vxpsiirgk(igk, ist3)=vxpsiirgk(igk, ist3)+vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega) ! pace IR
               End Do
 
!write(*,*)omp_get_thread_num(), "ist2=", ist2, "ist3=", ist3
               zvclmt(:,:,:,ist3)=zvclmt(:,:,:,ist3)+prod%mtrlm(:,:,:,1)*wkptnr(jk)
               
            End Do ! ist2

             !Do igk=1, Gkqset%ngk (1, ik)
             !     vxpsiirgk(igk, ist3)=vxpsiirgk(igk, ist3)+vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega) ! pace IR
             !  End Do
!!$OMP END DO NOWAIT
!!$OMP END DO
!zvclmt(:,:,:,ist3)=zvclmt(:,:,:,ist3)+prod%mtrlm(:,:,:,1)*wkptnr(jk)

         End Do ! ist3
        !$OMP END DO         
         call WFRelease(prod)
         call WFRelease(pot)
         Deallocate (wf1ir)
         Deallocate (wf2ir)
         Deallocate (prodir)
         Deallocate (potir)

!$OMP END PARALLEL

call timesec(ta)
write(*,*) 'omp loop',ta-tb


         vxpsimt=vxpsimt+zvclmt
      End Do ! non-reduced k-point set

!----------------------------------------------!
!     valence-core-valence contribution        !
!----------------------------------------------!
call timesec(ta)
      zvclmt (:, :, :, :) = 0.d0
            
If (.true.) Then
      Do is = 1, nspecies
         nrc = nrcmt(is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist2 = 1, spnst (is) !This is essentially ncore(is)
               If (spcore(ist2, is)) Then
                  l = spl(ist2, is)
                  norm = sqrt(0.5d0*spocc(ist2,is)/dble(2*l+1))
                  Do m = -l, l
                     lm = idxlm (l, m)
                     Do ir = 1, nrmt(is)
                        uir = norm * rwfcr (ir, 1, ist2, ias) / spr (ir, is)
                        wfcr1 (1:ntpll, ir) = uir * zbshthf (1:ntpll, lm)
                     End Do ! ir

                     ! Begin loop over occupied and empty states
                     Do ist3 = 1, nstsv

! calculate the complex overlap density
                        Call vnlrhomt2 (.true., is, wfcr1(:, :), &
                        & wf1%mtmesh(:, :, ias, ist3), zrhomt(:, :, &
                        & ias)) ! Psi*_{a}.Psi_{nk} = rho_{a;nk}; Returns in SH)
      
! calculate the Coulomb potential
                        Call zpotclmt (input%groundstate%ptnucl, &
                        & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                        & 0.d0, lmmaxvr, zrhomt(:, :, ias), zfmt) ! Returns SH
      
                        Call zgemm ('N', 'N', ntpll, nrc, lmmaxvr, &
                        & zone, zbshthf, ntpll, zfmt, lmmaxvr, &
                        & zzero, zvcltp, ntpll) ! Returns zvcltp in SC
      
                        zvcltp=conjg(zvcltp)
      
! calculate the complex overlap density
                        Call vnlrhomt2 (.true., is, zvcltp, wfcr1(:, :), zrhomt(:, :, ias)) ! Returns in SH
      
                        zvclmt(:,:,ias,ist3)=zvclmt(:,:,ias,ist3)+zrhomt(:, :, ias)
                        
                     End Do ! ist3
                  End Do ! m
               End If ! spcore(ist2, is)
            End do ! ist2
         End Do ! ia
      End Do ! is
End If
call timesec(tb)
write(*,*) 'vcv',tb-ta

      vxpsimt=vxpsimt+zvclmt
      if (.true.) then
      Allocate (wf1ir(ngrtot))
call timesec(ta)
      Do ist1 = 1, nstsv

         If ((ist1.le.nomax).and.(q0corr.ne.0.d0)) Then
            ! Evaluate wavefunction in real space
            wf1ir(:) = 0.d0
            Do igk = 1, Gkqset%ngk (1, ik)
               ifg = igfft (Gkqset%igkig(igk, 1, ik))
               wf1ir(ifg) = t1*wf1%gk(igk, ist1)
            End Do
            Call zfftifc (3, ngrid, 1, wf1ir(:))

            ! Apply q=0 correction to MT part
            vxpsimt(:,:,:,ist1) = vxpsimt(:,:,:,ist1) + q0corr*wf1%mtrlm(:,:,:,ist1)

            ! Apply correction to IR part and roll back correction to momentum space
            vxpsiirtmp(:) = q0corr*wf1ir(:)*cfunir(:)
            Call zfftifc (3, ngrid,-1,vxpsiirtmp)
            Do igk=1, Gkqset%ngk (1, ik)
               vxpsiirgk(igk, ist1)=vxpsiirgk(igk, ist1)+vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega) ! pace IR
            End Do
         End If 

         ! Write(*,*) 'ist1=',ist1
         Do ist3 = 1, nstsv

            ztir = 0.d0
            Do igk = 1, Gkqset%ngk (1, ik)
               ztir = ztir + conjg(wf1%gk(igk, ist1))*vxpsiirgk(igk,ist3)
            End Do
            ztmt = zfinpmt (.True., wf1%mtrlm(:,:,:,ist1),vxpsimt(:,:,:,ist3))
            vnlvv (ist1, ist3) = vnlvv (ist1, ist3) - ztmt - ztir

         End Do ! ist3
      End Do ! ist1
call timesec(tb)

write(*,*) 'Matrix',tb-ta
end if
      if(.false.) then
call timesec(ta)
      If (q0corr.ne.0.d0) Then
        Allocate (wf1ir(ngrtot))
        Do ist1 = 1, nomax
            ! Evaluate wavefunction in real space
            wf1ir(:) = 0.d0
            Do igk = 1, Gkqset%ngk (1, ik)
               ifg = igfft (Gkqset%igkig(igk, 1, ik))
               wf1ir(ifg) = t1*wf1%gk(igk, ist1)
            End Do
            Call zfftifc (3, ngrid, 1, wf1ir(:))

            ! Apply q=0 correction to MT part
            vxpsimt(:,:,:,ist1) = vxpsimt(:,:,:,ist1) + q0corr*wf1%mtrlm(:,:,:,ist1)

            ! Apply correction to IR part and roll back correction to momentum space
            vxpsiirtmp(:) = q0corr*wf1ir(:)*cfunir(:)
            Call zfftifc (3, ngrid,-1,vxpsiirtmp)
            Do igk=1, Gkqset%ngk (1, ik)
               vxpsiirgk(igk, ist1)=vxpsiirgk(igk, ist1)+vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega) ! pace IR
            End Do
        End Do ! ist1
        Deallocate (wf1ir)
      End If 
call timesec(tb)
write(*,*)'q0corr', tb-ta







!----------------------------------------
call timesec(ta)
! IR part
! replace with ZGEMM
!if (.t.) then
      allocate(ztir1(nstsv, nstsv))
      call zgemm('C', 'N', nstsv, nstsv, Gkqset%ngk (1, ik), dcmplx(1.0D0,0.0), wf1%gk, &
              & nmatmax, vxpsiirgk, ngkmax, dcmplx(0.0D0,0.0d0), ztir1, nstsv)
      vnlvv = vnlvv - ztir1
     ! do ist1 = 1, 6
     !    write(*,'(14F13.9)') dble(vnlvv(ist1,1:6))
     ! end do

      deallocate(ztir1)
call timesec(tb)
      if (.false.) then
write(*,*)"zgemm", tb-ta
      Do ist1 = 1, nstsv
         Do ist3 = 1, nstsv

            ztir = 0.d0
            Do igk = 1, Gkqset%ngk (1, ik)
               ztir = ztir + conjg(wf1%gk(igk, ist1))*vxpsiirgk(igk,ist3)
            End Do
            vnlvv (ist1, ist3) = vnlvv (ist1, ist3) - ztir
                write(*,*)ist1, ist3, ztir
         End Do ! ist3
      End Do ! ist1

 !     do ist1 = 1, 6
 !        write(*,'(14F13.9)') dble(vnlvv(ist1,1:6))
 !     end do
   end if
!stop
      
! MT part
call timesec(ta)
!write(*,*)"before omp"
!!$omp parallel default(shared), private(ist3, ztmt)
      Do ist1 = 1, nstsv
         ! Write(*,*) 'ist1=',ist1
!!$omp do
         Do ist3 = 1, nstsv

            ztmt = zfinpmt (.True., wf1%mtrlm(:,:,:,ist1),vxpsimt(:,:,:,ist3))
            vnlvv (ist1, ist3) = vnlvv (ist1, ist3) - ztmt

         End Do ! ist3
!!$omp end do
      End Do ! ist1
!!$omp end parallel
call timesec(tb)
end if
do ist1 = 1, 6
         write(*,'(14F13.9)') dble(vnlvv(ist1,1:6))
 end do
!stop
!write(*,*) 'omp mt',tb-ta

if (.false.) then
      Write(*,*) "ikp, ik, memopt:", ikp, ik
      write(*,*) 'vnlvv real (1:14,1:14)'
      do ist1 = 1, 14
         write(*,'(14F13.9)') dble(vnlvv(ist1,1:14))
      end do
      
      write(*,*) 'vnlvv imag (1:14,1:14)--'
      do ist1 = 1, 14
         write(*,'(14F13.9)') dimag(vnlvv(ist1,1:14))
      end do
end if
      
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
      Deallocate (ylmgq, sfacgq)
      Deallocate (wfcr1)
      Deallocate (zrhomt, zrhoir, zvcltp, zfmt)
      call WFRelease(wf1)
      call WFRelease(wf2)
      call WFRelease(prod)
      Return
End Subroutine
!EOC
