!
!
!
! Copyright (C) 2021 D. Zavickis, A. Gulans
!
!
Subroutine FockExchange (ikp, q0corr, vnlvv, vxpsiirgk, vxpsimt)
      use modbess, only: nfit, zbessi,zbessk,erfc_fit,zilmt,lambda
      use modinteg
      Use modmain 
      Use modinput
      Use modgw, only : kqset,Gkqset, kset, nomax, numin, ikvbm, ikcbm, ikvcm, Gset
      Use potentials, only: coulomb_potential
      use weinert, only: poisson_mt_yukawa
      USE OMP_LIB

      use poterf

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
      Integer :: ifg, ngvec1, ifit
      Logical :: solver, cutoff, handleG0

      Real (8) :: v (3), cfq, ta,tb, t1, norm, uir, x, gmax_pw_method,kvec(3)
      Complex (8) zrho01, zrho02, ztmt,zt1,zt2,zt3,zt4, ztir
      Integer :: nr, l, m, io1, lm2, ir, if3, j, lmaxvr

      Complex (8) ::  potmt0(lmmaxvr, nrcmtmax, natmtot), potir0(ngrtot),rhoG0,potG0

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
      Complex (8), Allocatable :: zfmt (:, :),zfmt0 (:, :)
      Complex (8), Allocatable :: zwfir (:)

      real(8), allocatable :: jlgqsmallr(:,:,:,:),jlgrtmp(:)
      

      type (WFType) :: wf1,wf2,prod,pot
! external functions
      Complex (8) zfinp, zfmtinp, zfinpir, zfinpmt
      External zfinp, zfmtinp, zfinpir, zfinpmt

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
      Allocate (zfmt(lmmaxvr, nrcmtmax),zfmt0(lmmaxvr, nrcmtmax))
      Allocate (zvclmt(lmmaxvr, nrcmtmax, natmtot, nstsv))
      Allocate (zwfir(ngkmax))


      write(*,*) "erfcapprox=",input%groundstate%hybrid%erfcapprox

      lmaxvr=input%groundstate%lmaxvr
      write(*,*)"lambda",lambda




     

      if (allocated(evalfv)) deallocate(evalfv)
      allocate(evalfv(nstfv,kset%nkpt))
      
      ! test comment for push

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

         kvec=v
         write(*,*)"kvec",kvec



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

!!!variables for the erfcapprox="PW"
      if (input%groundstate%hybrid%erfcapprox.eq."PW")then 
         write(*,*)"shortest g+q vec",igq0
         gmax_pw_method=6d0*dsqrt(2d0)*lambda ! 6*sigmas of a Gaussian
         write(*,*)"min(g+q): G(",igq0,")=", gqc(igq0),"Largest G vector for PW method:",gmax_pw_method
         do ig=1, ngvec
               if (gqc(ig).lt.gmax_pw_method) then ! 6*sigmas of a Gaussian
                  ngvec1=ig+1
                  write(*,*)"G(",ig,")=",gqc(ig)
                  !exit
               endif
         enddo
   write(*,*)"ngvec1=",ngvec1


         if (allocated(jlgqsmallr)) then
            deallocate(jlgqsmallr)
         endif
         if (allocated(jlgrtmp)) then
            deallocate(jlgrtmp)
         endif
         allocate(jlgqsmallr(nrcmtmax,0:lmaxvr,ngvec1,nspecies))
         allocate(jlgrtmp(0:lmaxvr))
         
         do is=1,nspecies
            do ig=1,ngvec1
              do ir=1,nrcmt(is)
                x=gqc(ig)*rcmt(ir,is)
                call sbessel(lmaxvr,x,jlgrtmp)
                jlgqsmallr(ir,:,ig,is)=jlgrtmp(:)
              enddo
            enddo   
          enddo  

         endif


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


         if ((input%groundstate%hybrid%erfcapprox.eq."truncatedYukawa").or.&
         & ((input%groundstate%hybrid%erfcapprox.eq."none").and.(input%groundstate%hybrid%singularity.eq."exc0d")))then 
            cutoff=.True.
         else 
            cutoff=.False.
         endif
         if ((input%groundstate%hybrid%erfcapprox.eq."truncatedYukawa").or.&
            & ((input%groundstate%hybrid%erfcapprox.eq."none").and.(input%groundstate%hybrid%singularity.eq."exc0d")).or.&
            & (input%groundstate%hybrid%erfcapprox.eq."Yukawa") ) then
            handleG0=.false.
         else
            if (ik.eq.jk) then
               handleG0=.true.
            else
               handleG0=.false.
            endif
         endif

         write(*,*)"cutoff", cutoff,"ik, jk",ik,jk,"handleG0",handleG0

!write(*,*)"pirms", OMP_GET_THREAD_NUM()
write(*,*)"nomax",nomax,"nstfv",nstfv
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist3,wf1ir,wf2ir,igk,ifg,prod,prodir,zrho01,pot,potir,vxpsiirtmp,potmt0,potir0,j,rhoG0) REDUCTION(+:zvclmt,vxpsiirgk)
!write(*,*)"pēc", OMP_GET_THREAD_NUM()
         call WFInit(prod)
         call WFInit(pot)

         Allocate(pot%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
         Allocate(prod%mtrlm(lmmaxvr,nrmtmax,natmtot,1))
         Allocate (wf1ir(ngrtot))
         Allocate (wf2ir(ngrtot))
         Allocate (prodir(ngrtot))
         Allocate (potir(ngrtot))


         Do ist2 = 1, nomax

            wf2ir(:) = 0.d0
            Do igk = 1, Gkqset%ngk (1, jk)
               ifg = igfft (Gkqset%igkig(igk, 1, jk))
               wf2ir(ifg) = t1*wf2%gk(igk, ist2)
            End Do
            Call zfftifc (3, ngrid, 1, wf2ir(:))

!$OMP DO SCHEDULE(DYNAMIC)
            Do ist3 = 1, nstfv

               vxpsiirtmp(:) = 0.d0
               wf1ir(:) = 0.d0 
               Do igk = 1, Gkqset%ngk (1, ik)
                  ifg = igfft (Gkqset%igkig(igk, 1, ik))
                  wf1ir(ifg) = t1*wf1%gk(igk, ist3)
               End Do
               Call zfftifc (3, ngrid, 1, wf1ir(:))

   ! calculate the complex overlap density
   !-----------------------------------------------------------------------------------

               call WFprodrs(ist2,wf2,ist3,wf1,prod)
               prodir(:)=conjg(wf2ir(:))*wf1ir(:)


               


               !if ((ik.eq.jk).and.(.not.solver)) then
               if (handleG0) then
                  Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                  & igq0), sfacgq0, prod%mtrlm(:,:,:,1), prodir(:), rhoG0) 
                  prodir(:)=prodir(:)-rhoG0
                  prod%mtrlm(1,:,:,1)=prod%mtrlm(1,:,:,1)-rhoG0/y00
               endif


  
               pot%mtrlm(:,:,:,1)=zzero
               potir=zzero
               potmt0=zzero
               potir0=zzero
               if ((input%groundstate%hybrid%erfcapprox.eq."truncatedYukawa").or.(input%groundstate%hybrid%erfcapprox.eq."Yukawa")) then !Yukawa case
                  do j=1, nfit
                     Call coulomb_potential (nrcmt, rcmt, ngvec, gqc, igq0, &
                     & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:,:,:,1), &
                     & prodir(:), potmt0, potir0, zrho02, &
                     & cutoff=cutoff,hybrid_in=.true.,yukawa_in=.true., &
                     & zlambda_in=erfc_fit(j,2),zbessi=zbessi(:,j,:,:),zbessk=zbessk(:,j,:,:),zilmt=zilmt(j,:,:),kvec_in=kvec)
                     pot%mtrlm(:,:,:,1)=pot%mtrlm(:,:,:,1)+potmt0 * erfc_fit(j,1)
                     potir=potir+potir0 * erfc_fit(j,1)
                  enddo
                  do j=2, nfit
                     Call coulomb_potential (nrcmt, rcmt, ngvec, gqc, igq0, &
                     & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:,:,:,1), &
                     & prodir(:), potmt0, potir0, zrho02, &
                     & cutoff=cutoff,hybrid_in=.true.,yukawa_in=.true., &
                     & zlambda_in=conjg(erfc_fit(j,2)),zbessi=conjg(zbessi(:,j,:,:)),zbessk=conjg(zbessk(:,j,:,:)),zilmt=conjg(zilmt(j,:,:)),kvec_in=kvec)
                     pot%mtrlm(:,:,:,1)=pot%mtrlm(:,:,:,1)+potmt0 * conjg(erfc_fit(j,1))
                     potir=potir+potir0 * conjg(erfc_fit(j,1))
                  enddo
               endif
   

               if ((input%groundstate%hybrid%erfcapprox.eq."PW").or.(input%groundstate%hybrid%erfcapprox.eq."none")) then
                  Call coulomb_potential (nrcmt, rcmt, ngvec, gqc, igq0, &
                  & jlgqr, ylmgq, sfacgq, zn, prod%mtrlm(:,:,:,1), &
                  & prodir(:), potmt0, potir0, zrho02, &
                  & cutoff=cutoff, hybrid_in=.true.,kvec_in=kvec)

                  if (input%groundstate%hybrid%erfcapprox.ne."PW")then 
                     pot%mtrlm(:,:,:,1)=potmt0
                     potir=potir0
                  endif
               endif

if (input%groundstate%hybrid%erfcapprox.eq."PW")then  
   call poterfpw(ngvec1, prodir,prod%mtrlm(:,:,:,1),igfft,sfacgq,ylmgq,gqc,jlgqsmallr,potir, pot%mtrlm(:,:,:,1))
   potir=potir0 - potir !Coulomb - erf
   pot%mtrlm(:,:,:,1)=potmt0 - pot%mtrlm(:,:,:,1)
endif




              ! call WFprodrs(ist2,wf2,ist3,wf1,prod)
              ! prodir(:)=conjg(wf2ir(:))*wf1ir(:)

               !if ((ik.eq.jk).and.(.not.solver)) then
               if (handleG0) then
                  Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:, &
                  & igq0), sfacgq0, pot%mtrlm(:,:,:,1), potir(:), zrho01)
                  potir(:)=potir(:)-zrho01
                  pot%mtrlm(1,:,:,1)=pot%mtrlm(1,:,:,1)-zrho01/y00
               endif

               if ((input%groundstate%hybrid%erfcapprox.eq."PW").and.(handleG0))then
                  potG0 = rhoG0*pi/lambda**2
                  potir(:)=potir(:)+potG0
                  pot%mtrlm(1,:,:,1)=pot%mtrlm(1,:,:,1)+potG0/y00
               endif

   !-----------------------------------------------------------------------------------1
               call genWFonMeshOne(pot)
               pot%mtmesh=conjg(pot%mtmesh)
               call WFprodrs(1,pot,ist2,wf2,prod)
               vxpsiirtmp(:) = potir(:)*wf2ir(:)*wkptnr(jk)*cfunir(:)

   ! ----------------------------------------------------------------------------------
               ! Calculate Fourier transform of vxpsiirtmp(r)
               Call zfftifc (3, ngrid,-1,vxpsiirtmp)


do j=1,10
!write(*,*) ist3,",",OMP_GET_THREAD_NUM(),",",j,"," ,dble(vxpsiirtmp(j)),",",imag(vxpsiirtmp(j))
enddo


               Do igk=1, Gkqset%ngk (1, ik)
                  vxpsiirgk(igk, ist3)=vxpsiirgk(igk, ist3)+vxpsiirtmp(igfft(Gkqset%igkig(igk, 1, ik)))*sqrt(Omega) ! pace IR
               End Do

               zvclmt(:,:,:,ist3)=zvclmt(:,:,:,ist3)+prod%mtrlm(:,:,:,1)*wkptnr(jk)

            End Do ! ist3
!$OMP END DO NOWAIT
         End Do ! ist2

         call WFRelease(prod)
         call WFRelease(pot)
         Deallocate (wf2ir)
         Deallocate (prodir)
         Deallocate (potir)
         Deallocate (wf1ir)
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
                     if (input%groundstate%hybrid%erfcapprox.eq."none") then
                        Call zpotclmt (input%groundstate%ptnucl, &
                        & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                        & 0.d0, lmmaxvr, zrhomt(:, :, ias), zfmt) ! Returns SH     
                     else
                        zfmt=zzero
                        do ifit=1, nfit
                           call poisson_mt_yukawa( lmaxvr, nrc, rcmt(:, is), zrhomt(:, :, ias), zfmt0, &
                              & erfc_fit(ifit,2), zbessi(:,ifit,:,is), zbessk(:,ifit,:,is), is)
                           zfmt=zfmt+zfmt0 * erfc_fit(ifit,1)
                        enddo
                        do ifit=2, nfit
                           call poisson_mt_yukawa( lmaxvr, nrc, rcmt(:, is), zrhomt(:, :, ias), zfmt0, &
                              & conjg(erfc_fit(ifit,2)), conjg(zbessi(:,ifit,:,is)), conjg(zbessk(:,ifit,:,is)), is)
                           zfmt=zfmt+zfmt0 * conjg(erfc_fit(ifit,1))
                        enddo
                     endif


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
      Allocate (wf1ir(ngrtot))
call timesec(ta)
      Do ist1 = 1, nstsv
         if(.false.)then
         write(*,*)"q0corr",q0corr
         !If ((ist1.le.nomax).and.(q0corr.ne.0.d0)) Then
            ! Evaluate wavefunction in real space
            wf1ir(:) = 0.d0
            Do igk = 1, Gkqset%ngk (1, ik)
               ifg = igfft (Gkqset%igkig(igk, 1, ik))
               wf1ir(ifg) = t1*wf1%gk(igk, ist1)
            End Do
            Call zfftifc (3, ngrid, 1, wf1ir(:))

!!!!!!!!!!!!!!!!!!!!!!!!!! KOREKCIJA ko vajadzētu atslēgt Aux funct method
write(*,*)"FockExchange korekcija notiek"

            ! Apply q=0 correction to MT part
            vxpsimt(:,:,:,ist1) = vxpsimt(:,:,:,ist1) + q0corr*wf1%mtrlm(:,:,:,ist1)

            ! Apply correction to IR part and roll back correction to momentum space
            vxpsiirtmp(:) = q0corr*wf1ir(:)*cfunir(:)

!!!!!!!!!!!!!!!!!!!!!!!!!! KOREKCIJA ko vajadzētu atslēgt Aux funct method//

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

if (.false.) then
      write(*,*)"element(1,1)",dble(vnlvv(1,1)),",",imag(vnlvv(1,1))
      Write(*,*) "ikp, ik, memopt:", ikp, ik
      write(*,*) 'vnlvv real (1:nstfv,1:nstfv)'
      do ist1 = 1, min(nstfv,14)
         write(*,'(14F13.9)') dble(vnlvv(ist1,1:min(nstfv,14)))
      end do
      
      write(*,*) 'vnlvv imag (1:nstfv,1:nstfv)--'
      do ist1 = 1, min(nstfv,14)
         write(*,'(14F13.9)') dimag(vnlvv(ist1,1:min(nstfv,14)))
      end do
end if
      
      Deallocate (vgqc, tpgqc, gqc, jlgqr, jlgq0r)
      Deallocate (ylmgq, sfacgq)
      Deallocate (wfcr1)
      Deallocate (wf1ir)
      Deallocate (zrhomt, zrhoir, zvcltp, zfmt, zfmt0)

      if (allocated(jlgqsmallr)) then
         deallocate(jlgqsmallr)
      endif
      if (allocated(jlgrtmp)) then
         deallocate(jlgrtmp)
      endif


      call WFRelease(wf1)
      call WFRelease(wf2)
      call WFRelease(prod)

 ! write(*,*)"FockExchange.f90 stop"
 ! stop
      Return
End Subroutine
!EOC
