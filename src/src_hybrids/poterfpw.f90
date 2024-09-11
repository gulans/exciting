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
      integer :: lmaxvr,is,ia,ias,ig,ifg,l,m,lm,ir,nr
      complex (8), allocatable :: rhorcmt2(:,:)

      lmaxvr=input%groundstate%lmaxvr

potir(:) = cfunir(:)*rhoir(:)

Call zfftifc (3, ngrid, -1, potir(:))  !to G space

!!! Obtain Fourier coeficients of the density in the MT and add it in porir
Allocate(rhorcmt2(lmmaxvr,nrmtmax))
do is=1,nspecies
   nr=nrmt(is)
   do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nr
        do lm=1,lmaxvr
          rhorcmt2(lm,ir)=rcmt(ir,is)**2 * rhomt(lm,ir,ias)
        enddo
!       zt1=4*pi*rcmt(ir,is)**2/omega
!       do l=0,lmaxvr
!         zt2=zt1/zil(l)
!         do m=-l,l
!           lm=idxlm(l,m)
!           rhorcmt2(lm,ir)=zt1 * rhomt(lm,ir,ias)
!         enddo
!        enddo
      enddo

      do ig=1,ngvec1
         ifg=igfft(ig)!(Gkqset%igkig(ig, 1, iq)) 
         do l=0,lmaxvr
            do m=-l,l 
               lm=idxlm(l,m)
               zfmt1(:nr)=jlgqsmallr(:nr,l,ig,is)*rhorcmt2(lm,ir) !rcmt(:nr,is)**2* rhomt(lm,:nr,ias)
               call integ_cv (nr, is, zfmt1(:nr), zt3, mt_integw)
               zt4=zt3*4d0*pi*ylmgq(lm,ig)*conjg(sfacgq(ig, ias))/(omega*zil(l)) !!!Fāzes reizinātājs sfacgq(ig, ias) ??
!              zt4=zt3*ylmgq(lm,ig)*conjg(sfacgq(ig, ias)) !!!Fāzes reizinātājs sfacgq(ig, ias) ??
               potir(ifg)=potir(ifg)+zt4
            enddo ! m
         enddo ! l
      enddo ! ig
   enddo ! ia 
enddo ! is
Deallocate(rhorcmt2)

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
   

!!!obtain radial MT functions from potir and store in pot%mtrlm(:,:,:) (lm,ir,ias)
   potmt(:,:,:)=zzero
   do is=1,nspecies
      nr=nrmt(is)
      do ia=1,natoms(is)
         ias=idxas(ia,is)
         do ir=1,nr
            do ig=1, ngvec1
               ifg =igfft(ig)! igfft(Gkqset%igkig(ig, 1, iq))
               zt1=potir(ifg)*sfacgq(ig, ias)
               do l=0,lmaxvr
                  zt2=zt1*jlgqsmallr(ir,l,ig,is) 
                  do m=-l,l                      
                     lm=idxlm(l,m)                     
                     potmt(lm,ir,ias)=potmt(lm,ir,ias)+zt2*conjg(ylmgq(lm,ig))
                  enddo 
               enddo
            enddo
            do l=0,lmaxvr
              do m=-l,l
                lm=idxlm(l,m)
                potmt(lm,ir,ias)=4d0*pi*zil(l)*potmt(lm,ir,ias)
              enddo
            enddo
         enddo
      enddo
   enddo

   Call zfftifc (3, ngrid, 1, potir(:)) !to realspace 



end subroutine

end module poterf
