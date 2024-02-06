module poterf
        
      implicit none
      private
    
      public :: poterfpw
    
      contains

subroutine poterfpw (ngvec1, rhoir,rhomt,igfft,sfacgq,ylmgq,gqc, jlgqsmallr,potir, potmt)
      use modbess, only: lambda
      use mod_atoms, only: nspecies,natoms, natmtot, idxas
      use mod_lattice, only: omega
      use modinput
      use mod_muffin_tin, only: rcmt,nrmtmax,nrmt,idxlm,lmmaxvr
      use modinteg
      use constants, only: fourpi,zil,pi,zzero
      use mod_Gvector, only: cfunir,ngrid,ngrtot,ngvec
      Implicit None
! arguments
      integer, Intent (In)      :: ngvec1
      complex (8), Intent (In)  :: rhoir(:)     !ngrtot
      complex (8), Intent (In)  :: rhomt(:,:,:) !lmmaxvr,nrmtmax,natmtot
      integer, intent(in)       :: igfft(:)     !ngrtot
      complex (8), Intent (In)  :: sfacgq(:,:)  !ngvec, natmtot
      real (8), Intent (In)     :: gqc(:)
      complex (8), Intent (In)  :: ylmgq(:, :)  !lmmaxvr, ngvec
      real (8), Intent (In)     :: jlgqsmallr(:,0:,:,:)  !nrmtmax,0:input%groundstate%lmaxvr,ngvec1,nspecies
      complex (8), Intent (Out) :: potir(:)      ! ngrtot
      complex (8), Intent (Out) :: potmt(:,:,:)  !lmmaxvr,nrmtmax,natmtot
      

      complex (8) :: zfmt1(nrmtmax),zfmt2(nrmtmax),zt1,zt2,zt3,zt4
      integer :: lmaxvr,is,ia,ias,ig,ifg,l,m,lm,ir

      lmaxvr=input%groundstate%lmaxvr

potir(:) = cfunir(:)*rhoir(:)

Call zfftifc (3, ngrid, -1, potir(:))  !to G space

!!! Obtain Fourier coeficients of the density in the MT and add it in porir
do is=1,nspecies
   do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ig=1,ngvec1
         ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq)) 
         do l=0,lmaxvr
            do m=-l,l 
               lm=idxlm(l,m)
               zfmt1=jlgqsmallr(:,l,ig,is)*rcmt(:,is)**2* rhomt(lm,:,ias)
               call integ_cf (nrmt(is), is, zfmt1, zfmt2, mt_integw)
               zt3=zfmt2(nrmt(is))
               zt4=zt3*4d0*pi*ylmgq(lm,ig)*conjg(sfacgq(ig, ias))/(omega*zil(l)) !!!Fāzes reizinātājs sfacgq(ig, ias) ??
               potir(ifg)=potir(ifg)+zt4
            enddo ! m
         enddo ! l
      enddo ! ig
   enddo ! ia 
enddo ! is



do ig=1, ngvec1
      ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq))   
      If (gqc(ig) .Gt. input%structure%epslat) Then  
         potir(ifg)=potir(ifg)*4d0*pi*exp(-gqc(ig)**2/(4d0*lambda**2))/gqc(ig)**2
      else !!!! erfc kernels G=0 *(-1)  
         potir(ifg)=-potir(ifg)*pi/lambda**2
      endif
   enddo

  

   do ig=ngvec1,ngrtot!ngvec
      ifg=igfft(ig)
      potir(ifg)=zzero
   enddo
   

!!!obtain radial MT functions from potir and store in pot%mtrlm(:,:,:,1) (lm,ir,ias)
   potmt(:,:,:)=zzero
   do is=1,nspecies
      do ia=1,natoms(is)
         ias=idxas(ia,is)
         do ir=1,nrmt(is)
            do ig=1, ngvec1
               ifg =igfft(ig)! igfft(Gkqset%igkig(ig, 1, iq))
               do l=0,lmaxvr
                  zt1=4d0*pi*potir(ifg)*zil(l)*jlgqsmallr(ir,l,ig,is) * sfacgq(ig, ias)
                  do m=-l,l                      
                     lm=idxlm(l,m)                     
                     zt2=zt1*conjg(ylmgq(lm,ig))
                     potmt(lm,ir,ias)=potmt(lm,ir,ias)+zt2
                 
                  enddo 
               enddo
            enddo
         enddo
      enddo
   enddo

   Call zfftifc (3, ngrid, 1, potir(:)) !to realspace 



end subroutine

end module poterf