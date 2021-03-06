!
!
!
! Copyright (C) 2002-2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine overlapsetup (system, ngp, apwalm, igpig, vgpc)
      Use modfvsystem
      Use modinput
      Use mod_eigensystem
      Use mod_atoms
      Use mod_timing
      Use mod_muffin_tin
      Use mod_APW_LO
      Use mod_gkvector
!
      Implicit None
      Type (evsystem) :: system
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Integer :: n
      Character (256) :: prefix
!local variables
      Integer, Save :: ikc
      Real (8), Save :: cputot
      Real (8) :: cpuaa, cpualo, cpulolo, cpui, cpu00, cpu01,ts0,ts1
      Integer :: i, is, ia
      Complex (8) v(1),viens
      Real (8) :: cpu0, cpu1
      Real (8) :: threshold
      Complex (8), allocatable :: apwi(:,:),zm(:,:),apwi2(:,:)
      integer if1,if3,l3,m3,lm3,io1,io2,ias,maxnlo,ilo,j1,j2,j3,lm1,lm2,j,io,l,ilo1,ilo2,l1


!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
!
!
      Call timesec (cpu0)
! set the matrices to zero
!
! muffin-tin contributions
      allocate(apwi(haaijSize,ngp))
      allocate(apwi2(ngp,haaijSize))
      allocate(zm(haaijSize,ngp))
      Do is = 1, nspecies
        Do ia = 1, natoms (is)
!--Hamiltonian--
! APW-APW part 
          Call timesec (ts0)
          ias = idxas (ia, is)
          apwi=dcmplx(0d0,0d0)
          apwi2=dcmplx(0d0,0d0)
          if3=0
          Do l3 = 0, input%groundstate%lmaxmat
            Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
              Do io2 = 1, apword (l3, is)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngp, io2, lm3, ias)
                apwi2(:,if3)=conjg(apwalm(1:ngp, io2, lm3, ias))
              End Do
            End Do
          End Do
          zm=zzero
          viens=dcmplx(1d0,0)


!--Overlap--
! APW-APW part
          Call timesec (ts0)

            deallocate(zm)
            allocate(zm(ngp,haaijSize))
            zm=zzero
            if3=0
            Do l3 = 0, input%groundstate%lmaxmat
              Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)
               
                Do io2 = 1, apword (l3, is)
                  Do io1 = 1, apword (l3, is)
                    zm(:,if3+io2)=zm(:,if3+io2)+h1aa(io1,io2,l3,ias)*apwi2(:,if3+io1)
                  enddo
                  zm(:,if3+io2)=zm(:,if3+io2)+apwi2(:,if3+io2)
                End Do
                if3=if3+apword (l3, is)
              End Do
            End Do
            call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'C', &           ! TRANSB = 'N'  op( B ) = B.
                        ngp, &          ! M ... rows of op( A ) = rows of C
                        ngp, &           ! N ... cols of op( B ) = cols of C
                        haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                        viens, &          ! alpha
                        apwi, &           ! A
                        haaijSize,&           ! LDA ... leading dimension of A
                        zm, &           ! B
                        ngp, &          ! LDB ... leading dimension of B
                        viens, &          ! beta
                        system%overlap%za, &  ! C
                        system%overlap%rank &      ! LDC ... leading dimension of C
                       )
             deallocate(zm)
             allocate(zm(haaijSize,ngp))

          Call timesec (ts1)
          time_olpaan=ts1-ts0+time_olpaan

!What if it is, say, LAPW calculation without any local orbitals?
        if (nlorb(is).ne.0) then
! APW-LO part
          Call timesec (ts0)
          Do ilo = 1, nlorb (is)
            l = lorbl (ilo, is)
            lm1 = idxlm (l,-l)
            lm2 = idxlm (l, l)
            j1 = ngp + idxlo (lm1, ilo, ias)
            j2 = ngp + idxlo (lm2, ilo, ias)
            Do io = 1, apword (l, is)
              system%overlap%za(1:ngp,j1:j2)=system%overlap%za(1:ngp,j1:j2)+conjg(apwalm(1:ngp, io, lm1:lm2, ias) * (oalo (io, ilo, ias)+h1loa(io, ilo, ias)))
            End Do
            do j=j1,j2
              system%overlap%za(j,1:ngp)=conjg(system%overlap%za(1:ngp,j))
            End Do
          End Do
          Call timesec (ts1)
          time_olpalon=ts1-ts0+time_olpalon
! LO-LO part
          Call timesec (ts0)
          Do ilo1 = 1, nlorb (is)
            l = lorbl (ilo1, is)
            Do ilo2 = 1, nlorb (is)
              If (lorbl(ilo2, is) .Eq. l) Then
                lm1=idxlm (l,-l)
                j1= ngp + idxlo (lm1, ilo1, ias)
                j2= ngp + idxlo (lm1, ilo2, ias)
                do lm2=idxlm (l,-l),idxlm (l, l)
                  system%overlap%za(j1+lm2-lm1,j2+lm2-lm1)=system%overlap%za(j1+lm2-lm1,j2+lm2-lm1)+dcmplx(ololo(ilo1, ilo2, ias)+h1lolo(ilo1, ilo2, ias),0d0)
                enddo
              End If
            End Do
          End Do
          Call timesec (ts1)
          time_olplolon=ts1-ts0+time_olplolon
        endif

        End Do
      End Do
      deallocate(apwi,apwi2,zm)


! interstitial contributions
       Call timesec (ts0)
      Call olpistln (system%overlap, ngp, igpig, vgpc)
       Call timesec (ts1)
       time_olpistln=ts1-ts0+time_olpistln

       If ( .Not. ispacked(system%hamilton)) Then
          Call hamiltonoverlapocopy_UL (system)
       End If


#ifdef DEBUGHO
      Write (*,*) "apwalm", apwalm
      prefix = "H"
      Call HermitianMatrixToFiles (system%hamilton, prefix)
      prefix = "O"
      Call HermitianMatrixToFiles (system%overlap, prefix)
      Write (*,*) "wrote"
      Stop
#endif
!
      Call timesec (cpu1)
      timemat = timemat + cpu1 - cpu0

End Subroutine
