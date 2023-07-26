!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: findband
! !INTERFACE:
!
Subroutine findband1 (findlinentype, l, k, nr, r, vr, de0, etol, e0, tfnd,is,ia)
! !INPUT/OUTPUT PARAMETERS:
!   l   : angular momentum quantum number (in,integer)
!   k   : quantum number k, zero if Dirac eqn. is not to be used (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   de0 : default energy step size (in,real)
!   e   : input energy and returned band energy (inout,real)
! !DESCRIPTION:
!   Finds the band energies for a given radial potential and angular momentum.
!   This is done by first searching upwards in energy until the radial
!   wavefunction at the muffin-tin radius is zero. This is the energy at the top
!   of the band, denoted $E_{\rm t}$. A downward search is now performed from
!   $E_{\rm t}$ until the slope of the radial wavefunction at the muffin-tin
!   radius is zero. This energy, $E_{\rm b}$, is at the bottom of the band. The
!   band energy is taken as $(E_{\rm t}+E_{\rm b})/2$. If either $E_{\rm t}$ or
!   $E_{\rm b}$ cannot be found then the band energy is set to the input value.
!
! !REVISION HISTORY:
!   Created September 2004 (JKD)
!EOP
!BOC
      use modmain,      only : iscl, efermi, tlast
      Implicit None
  ! arguments
      Character (*), Intent (In) :: findlinentype
      Integer, Intent(In)    :: l
      Integer, Intent(In)    :: k
      Integer, Intent(In)    :: nr

      Integer, Intent(In)    :: is
      Integer, Intent(In)    :: ia

      Real(8), Intent(In)    :: r (nr)
      Real(8), Intent(In)    :: vr (nr)
      Real(8), Intent(In)    :: de0
      real(8), intent(in)    :: etol
      Real(8), Intent(Inout) :: e0
      logical, intent(out)   :: tfnd
  ! local variables
  ! maximum number of steps
      Integer, Parameter :: maxstp = 1000
      Character * 4, Parameter :: emain = 'continue'
      Real (8), Parameter :: etoolow = - 1000.d0
      Real (8), Parameter :: ecutlow = - 100.d0
      Real (8), Parameter :: efermibands = 0.5d0
      Real (8), Parameter :: erangebands = 20.d0 !default 1d0 es nomainīju
      Real (8), Parameter :: ediffusebands = efermibands + erangebands
      Real (8), Parameter :: erange1 = 3.d0
      Real (8), Parameter :: edefault1 = 1.d0
      Real (8), Parameter :: eupcut = 30.d0
      Real (8), Parameter :: edncut = -10.d0
      Integer :: ir,ie, nn
      Real (8) :: de, e, t0, t1, t00, t10, dl, dl0
      Real (8) :: u, uup, udn, upold, ddnold, dnold, du, dudn, dupold, &
     & duup, utst, dutst
      Real (8) :: e1, e2, eidn, eiup
  ! automatic arrays
      Real (8) :: p0 (nr), p1 (nr), q0 (nr), q1 (nr)
      character(10) :: fname
      character(1024) :: message
     character(len=1024) :: filename
real(8) :: estart
estart=e0
      tfnd=.false.
      de = Abs(de0)


      write(*,*)"findband1 begins l=",l,"e0=",e0


      Select Case (trim(findlinentype))

!-------------------------------------------
      Case ('logderiv')
!-------------------------------------------
         
!        find the Linearization Energy from the equation D_{l}(E)=-(l+1)/R_{MT}
         e = e0
         call rschroddme(0, l, k, e, nr, r, vr, nn, p0, p1, q0, q1)

         t00 = p0(nr)
         t10 = p1(nr)
         if ( dabs(t00)>1.0d-8 ) then
           dl0 = r(nr)*t10/t00+(l+1)
         else
           dl0 = 1000.d0*dsign(1.d0,t00)*dsign(1.d0,t10)
         end if
         
         ! search for the solutions
         if ( dabs(dl0)<=1.0d-6 ) then
           continue
         else if ( dl0>0.d0 ) then
           e = e + de
         else if ( dl0<0.d0 ) then
           e = e - de
         end if 
           
         do while ( dabs(e) <= 100.d0 )
           call rschroddme(0, l, k, e, nr, r, vr, nn, p0, p1, q0, q1)
           t0 = p0(nr)
           t1 = p1(nr)
           if ( dabs(t0)>1.d-8 ) then
             ! shifted value of the logarithmic derivative
             dl = r(nr)*t1/t0+(l+1)
             if ( dl*dl0 <= 0.d0 ) then
               ! the node is found
               tfnd=.true.
               e0 = e
               exit
             else if ( dl>0.d0 ) then
               e = e + de
             else if ( dl<0.d0 ) then 
               e = e - de
             end if
           end if
         end do
         
!-------------------------------------------
      Case ('Wigner_Seitz')
!-------------------------------------------

! ! do ie=0,10
! ! write(*,*)e0+ie*0.1d0
! !    Call rschroddme2(is,ia,0, l, k, e0+ie*0.1d0, nr, r, vr, nn, p0, p1, q0, q1)
! !    WRITE(filename, '(a1,F5.2,a4)')'u',e0+ie*0.1d0,'.dat'
! !    open (11, file = filename, status = 'replace')
! !    Do ir = 1, nr
! !       write(11,*)r(ir),",",p0(ir)
! !    enddo
! !    close(11)

! !    WRITE(filename, '(a3,F5.2,a4)')'ugr',e0+ie*0.1d0,'.dat'
! !    open (11, file = filename, status = 'replace')
! !    Do ir = 1, nr
! !       write(11,*)r(ir),",",p1(ir)
! !    enddo
! !    close(11)
! ! enddo



! ! stop



!write(*,*)"e0",e0
         de = de0
         de= 0.02d0
         Call rschroddme2(is,ia,0, l, k, e0, nr, r, vr, nn, p0, p1, q0, q1)
         write(*,*)"findband1 e=",e0 ,"nn=",nn, "p0(nr)=",p0 (nr),"p1(nr)=", p1(nr)
         u = p0 (nr)
         du = p1 (nr)
         eiup = e0
         eidn = e0
         e1 = etoolow
         e2 = etoolow
         upold = u
         dupold = du
         dnold = u
         ddnold = du
11       Continue
!write(*,*)"eiup",eiup
         eiup = eiup + de
         Call rschroddme2(is,ia,0, l, k, eiup, nr, r, vr, nn, p0, p1, q0, q1)
         write(*,*)"findband1 up e=",eiup ,"nn=",nn, "p0(nr)=",p0 (nr),"p1(nr)=", p1(nr)
         uup = p0 (nr)
         duup = p1 (nr)
         utst = upold * uup
         dutst = dupold * duup
         upold = uup
         dupold = duup
         If (utst .Lt. 0.d0) Then
            e2 = eiup
            If (e1 .Gt. ecutlow) Go To 31
         End If
         If (dutst .Lt. 0.d0) Then
            e1 = eiup
            If (e2 .Gt. ecutlow) Go To 31
         End If

         !write(*,*)"findband1 e1,eiup,ecutlow,ediffusebands", e1,eiup,ecutlow,ediffusebands

         If ((e1 .Lt. ecutlow) .And. (eiup .Lt. ediffusebands)) Then
21          Continue
!write(*,*)"eidn",eidn

            eidn = eidn - de
            Call rschroddme2(is,ia,0, l, k, eidn, nr, r, vr, nn, p0, p1, q0, q1)
            write(*,*)"findband1 down e=",eidn ,"nn=",nn, "p0(nr)=",p0 (nr),"p1(nr)=", p1(nr)
            udn = p0 (nr)
            dudn = p1 (nr)
            utst = dnold * udn
            dutst = ddnold * dudn
            dnold = udn
            ddnold = dudn
            If (utst .Lt. 0.d0) Then
               e2 = eidn
               If (e1 .Gt. ecutlow) Go To 31
            End If
            If (dutst .Lt. 0.d0) Then
               e1 = eidn
               If (e2 .Gt. ecutlow) Go To 31
            End If
            If (e2 .Lt. ecutlow) Then
               Go To 11
            Else If (eidn .Gt. (e0-erange1)) Then
               Go To 21
            End If
         Else If (eiup .Lt. efermibands) Then
            Go To 11
         End If
31       Continue
         If ((e1 .Lt. ecutlow) .And. (e2 .Lt. ecutlow)) Then
            If (emain .Eq. 'stop') Then
               Go To 40
            Else
               e0 = edefault1
            End If
         Else If (e2 .Lt. ecutlow) Then
            If (emain .Eq. 'stop') Then
               Go To 40
            Else
               e0 = Max (e1, e0)
            End If
         Else If (e1 .Lt. ecutlow) Then
            Go To 40
         Else
            e0 = (e1+e2) * 0.5d0
            write(*,*)"findband1 result e1=",e1,"e2=",e2,"e=",e0
            Call rschroddme2(is,ia,0, l, k, e0, nr, r, vr, nn, p0, p1, q0, q1) !just to get nn to the write in the next line
            write(*,*)"findband1 final l=",l,"e0=",estart,"e=",e0,"nn=",nn
         End If
         tfnd=.true.
         Return
40       Continue
!        Print the warning
         call warning('Warning(findband):')
         Write(message, '(" No energy limits found for l=",i2)') l
         call warning(message)
         Write(message, '(" E-bottom ",g18.10,4x,"E-top ",g18.10)') e1, e2
         call warning(message)
         tfnd=.false.
         return
      Case Default
         Write (*,*)
         Write (*, '("Error(findband): No such method for search of lin&
        &earization energies: ",a)') trim (findlinentype)
         Write (*,*)
         Stop
      End Select

End Subroutine findband1
!EOC
