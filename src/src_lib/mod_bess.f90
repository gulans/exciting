module modbess 
use modinput
implicit none

integer :: nfit   
complex (8),allocatable :: erfc_fit(:,:)
complex (8),allocatable :: zbessi(:,:,:,:),zbessk(:,:,:,:)
complex (8),allocatable :: zilmt(:,:,:) ! (nfit, l, nspecies)
!> Stores all the needed arrays in IntegWeigthType data type variables for each species for MT grid and atom grid.

contains

subroutine init_bess(nrmtmax,nspecies,nrmt,r)
  use constants, only: zzero
  implicit none
  integer, intent(in) :: nrmtmax
  integer, intent(in) :: nspecies
  integer, intent(in) :: nrmt(nspecies)
  real(8), intent(in) :: r(nrmtmax,nspecies)
  real(8) :: lambda
  integer :: ii,ir, is
  nfit=9
  lambda=0.1d0
  allocate (erfc_fit(nfit,2))
  call errfun(nfit,lambda,erfc_fit)

  allocate (zbessi(nrmtmax,nfit,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies))
  allocate (zbessk(nrmtmax,nfit,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies))
  allocate (zilmt(nfit,0:input%groundstate%lmaxvr+input%groundstate%npsden+1,nspecies))
  zbessi=zzero
  zbessk=zzero

  do is=1, nspecies
    call get_Bess_fun( nrmt(is), r(1:nrmt(is),is),input%groundstate%lmaxvr+input%groundstate%npsden+1,&
      &nfit,erfc_fit,zbessi(1:nrmt(is),:,:,is),zbessk(1:nrmt(is),:,:,is))
  
  enddo

  do ii=1, nfit
    do is=1, nspecies
      zilmt(ii,:,is)=zbessi(nrmt(is),ii,:,is)
    enddo
  enddo



end subroutine

subroutine get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,bessi,bessk)
  implicit none
  real(8), PARAMETER :: Pi = 3.1415926535897932384d0
  integer, intent(in) :: Ngrid,Nrsfun,lmax
  real(8), intent(in) :: r(Ngrid)
  complex(8),intent(in) :: rsfunC(Nrsfun,2) 
  complex(8), intent(out) ::bessi(Ngrid,Nrsfun,0:lmax),bessk(Ngrid,Nrsfun,0:lmax)
  complex(8) ::z,zmin,zmax, bi,bk
  integer :: ir,l,fun
  integer :: ii
  real(8) :: remin,immin,remax,immax,zminabs,zmaxabs,zabs
    
    do l=0,lmax
    do fun=1,Nrsfun
    do ir = 1,Ngrid
      z=rsfunC(fun,2)*r(ir)
      call msbesselic (l, z, bi)
      call msbesselkc (l, z, bk)
      bessi(ir,fun,l) = bi
      bessk(ir,fun,l) = bk
  
    end do
    end do
    end do
  

  ! Debug
  ! write all argument values and bessel functions to files and print out range of Bessel function argument
  if(.false.) then
         
  ii=0
  zminabs=1d100
  zmaxabs=0d0
  open(11,file='besi_all.dat',status='replace')
  open(12,file='besk_all.dat',status='replace')
  
  write(11,*)"npk l r Re(z) Im(z) Re(Bessi) Im(Bessi)"
  write(12,*)"npk l r Re(z) Im(z) Re(Bessk) Im(Bessk)"
  
    do l=0,lmax
    do fun=1,Nrsfun
    do ir = 1,Ngrid
      ii=ii+1
      z=rsfunC(fun,2)*r(ir)
      zabs=dsqrt(realpart(z)**2+imagpart(z)**2)
      if (zabs.gt.zmaxabs)then
              zmaxabs=zabs
              zmax=z
      endif
      if (zabs.lt.zminabs)then
              zminabs=zabs
              zmin=z
      endif
      write(11,*)ii,l,r(ir),realpart(z),imagpart(z),realpart(bessi(ir,fun,l)),imagpart(bessi(ir,fun,l))
      write(12,*)ii,l,r(ir),realpart(z),imagpart(z),realpart(bessk(ir,fun,l)),imagpart(bessk(ir,fun,l))
  
  
    end do
    end do
    end do
  close(11)
  close(12)
  write(*,*)"Minimum argument for Bessel", zmin
  write(*,*)"Maximum argument for Bessel", zmax
  
  endif
  ! Debug end
end subroutine
  


subroutine errfun(n,mu,rsfunC)
!stores the coeficients of erfc expantion in array "rsfunC"
implicit none
integer, intent(in) :: n
real(8), intent(in) :: mu
complex(8),intent(out) :: rsfunC(n,2)
real(8) ::AA1(n),AA2(n),BB1(n),BB2(n)
integer :: i

if (n.eq.9)then
AA1=(/1.1793556928516621d01,-3.3997542606327369d00,-3.1999804759139012d00,&
1.2431945263044286d00,-1.5923036135878112d-02,-2.6330008531846357d-02,&
2.0107048788385695d-03,5.9871920957920612d-06,-1.9014193114960453d-06/)

AA2=(/0.0000000000000000d00,8.5884213008725414d00,-3.1024703144846724d00,&
-5.0090970998006024d-01,2.5412899494655228d-01,-1.5691616874704249d-02,&
-1.2226946313412218d-03,1.2240497422195631d-04,-2.1757399102663054d-06/)

BB1=(/4.5539394422559072d00,4.5539394422559072d00,4.5539394422559072d00,&
4.5539394422559072d00,4.5539394422559072d00,4.5539394422559072d00,&
4.5539394422559072d00,4.5539394422559072d00,4.5539394422559072d00/)

BB2=(/0.0000000000000000d00,9.4528885879869029d-01,1.8906858418787154d00,&
2.8358504737190153d00,3.7810807937266149d00,4.7286024942709624d00,&
5.6705170300130971d00,6.6184139539428530d00,7.7009425117754589d00/)
elseif (n.eq.8)then

AA1=(/1.2860987030345912d01*0.5d0,-4.9051058740247786d00,-1.9357963384430639d00,1.0117925923974511d00,&
        -9.9013582507628040d-02,-3.3813615287899379d-03,1.0330516726180475d-03,-2.2002738770173973d-05/)

AA2=(/0.0000000000000000d00,8.4204906883706165d00,-3.7444370420962709d00,1.5201429359841603d-01,&
        9.5997444318657268d-02,-1.1743959149248155d-02,1.5238534910174336d-04,5.8451000772024675d-06/)

BB1=(/4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00,&
        4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00/)

BB2=(/0.0000000000000000d00,1.0160216352396843d00,2.0413334512076715d00,3.0805139945811777d00,&
        4.1164481828964572d00,5.0646792984847675d00,5.9785623723578114d00,7.2027553450037622d00/)

elseif (n.eq.7)then
AA1=(/1.04388861772302022d01*0.5d0,-4.15230200875813171d00,-1.09637578338378838d00,&
5.85340334546933239d-01,-5.72214605455547975d-02,1.11445947001131431d-03,&
1.37005560127951071d-06/)

AA2=(/0.00000000000000000d00,6.36698813144878617d00,-2.73158403646281123d00,&
1.69412383167246605d-01,3.58625769268040781d-02,-3.56355734717054293d-03,&
7.06617335882476496d-05/)

BB1=(/4.42897206653721742d00,4.42897206653721742d00,4.42897206653721742d00,&
4.42897206653721742d00,4.42897206653721742d00,4.42897206653721742d00,&
4.42897206653721742d00/)

BB2=(/0.00000000000000000d00,1.07811959900154442d00,2.16055514404308591d00,&
3.25051246017242557d00,4.35355433640683831d00,5.50798572653768748d00,&
6.80680166786776208d00/)

elseif (n.eq.6)then

AA1=(/1.16801622501244324d01*0.5d0,-5.58814846386175379d00,4.00500965163793701d-02,&
2.27420591587985083d-01,-1.96862899952504854d-02,2.82940691484453110d-04/)

AA2=(/0.00000000000000000d00,5.91087741873710648d00,-2.71301646794361151d00,&
3.41139126788647828d-01,-1.37953981033003721d-02,2.25664278898155479d-04/)

BB1=(/4.44903651057751848d00,4.44903651057751848d00,4.44903651057751848d00,&
4.44903651057751848d00,4.44903651057751848d00,4.44903651057751848d00/)

BB2=(/0.00000000000000000d00,1.16136522308686563d00,2.34291495991631837d00,&
3.57030720633363519d00,4.88590055700419601d00,6.37901581267626661d00/)

elseif (n.eq.5)then
AA1=(/5.88649164984019890d00*0.5d0,-2.42999364511904048d00,-6.52318163977785026d-02,&
5.35112928141888433d-02,-1.53165601674370706d-03/)

AA2=(/0.00000000000000000d00,2.89879784656921213d00,-9.81481713958860635d-01,&
7.58072684236921479d-02,-1.61854017761636944d-03/)

BB1=(/4.03651233059081704d0,4.03651233059081704d0,4.03651233059081704d0,&
4.03651233059081704d0,4.03651233059081704d0/)

BB2=(/0.00000000000000000d0,1.27567126568311573d0,2.58394634761278397d0,&
3.97319627725742608d0,5.52839978282953215d0/)

elseif (n.eq.4)then

AA1=(/3.34075293484632052d00*0.5d0,-1.14940001089277311d00,-2.63828385336642122d-02,&
5.40643159457097119d-03/)


AA2=(/0.00000000000000000d00,1.50422252025890724d00,-3.29745255232857859d-01,&
1.22265802959771012d-02/)


BB1=(/3.62435558273901570d0,3.62435558273901570d0,3.62435558273901570d0,&
3.62435558273901570d0/)

BB2=(/0.00000000000000000d0,1.43494484832715830d0,2.93170108446797961d0,&
4.59735152907431743d0/)
elseif (n.eq.1)then

AA1=(/1d0/)
AA2=(/0d0/)
BB1=(/1d0/)
BB2=(/0d0/)

endif
   do i = 1,n
      rsfunC(i,1)=cmplx(AA1(i),AA2(i),8)
      rsfunC(i,2)=cmplx(BB1(i),BB2(i),8)
   end do
 rsfunC(:,2)=rsfunC(:,2)*mu
end subroutine

subroutine msbesselic (n, z, rez)
  Implicit None

  Integer, Intent (In) :: n
  Complex (8), Intent (In) :: z
  Complex (8), Intent (Out) :: rez
  Integer :: i
  Complex (8) :: g(-n-1:n)
  Complex (16) :: g16(-n-1:n),y,rez16
  Complex (8) :: rezi,rezk
real(8) :: absz
absz=sqrt( realpart(z)**2+imag(z)**2)


  if (n.gt.18)then
    write(*,*)"msbesselic.f90 Bessel functions with order higher than 18 are not supported"
    stop
    rez=cmplx(0d0,0d0,8)
    return

  endif

  
if ((n.eq.0).and. (abs(realpart(z)).lt.1d-3)) then
      rez=1d0+(1d0/6d0)*z**2 + (1d0/120d0)*z**4 + (1d0/5040d0)*z**6  

elseif ((n.eq.1) .and. (abs(realpart(z)).lt.1d-1).and.(abs(imagpart(z)).lt.1d-1) ) then
      rez=(1d0/3d0)*z +(1d0/30d0)*z**3 + (1d0/840d0)*z**5 + (1d0/45360d0)*z**7+&
              (1d0/3991680d0)*z**9

elseif ((n.eq.2) .and. (abs(realpart(z)).lt.0.5d0).and.(abs(imagpart(z)).lt.0.5d0))  then

      rez=(1d0/15d0)*z**2 +(1d0/210d0)*z**4 + (1d0/7560d0)*z**6 + (1d0/498960d0)*z**8 +&
              (1d0/51891840d0)*z**10+ (1d0/7783776000d0)*z**12+ (1d0/1587890304000d0)*z**14

elseif ((n.eq.3).and. (abs(realpart(z)).lt.1d-1)) then
      rez=(1d0/105d0)*z**3 + (1d0/1890d0)*z**5 + (1d0/83160d0)*z**7&
             + (1d0/6486480d0)*z**9 
elseif ((n.eq.4).and. (abs(realpart(z)).lt.1d-1)) then
      rez=(1d0/945d0)*z**4 + (1d0/20790d0)*z**6 + (1d0/1081080d0)*z**8&
            +(1d0/97297200d0)*z**10  
elseif ((n.eq.5).and. (abs(realpart(z)).lt.1d-1)) then
      rez=(1d0/10395d0)*z**5 + (1d0/270270d0)*z**7 + (1d0/16216200d0)*z**9&
             +(1d0/1654052400d0)*z**11 
elseif ((n.eq.6).and. (abs(realpart(z)).lt.1d-1)) then
      rez=(1d0/135135d0)*z**6 + (1d0/4054050d0)*z**8 + (1d0/275675400d0)*z**10&
              +(1d0/31426995600d0)*z**12
elseif ((n.eq.7).and. (abs(realpart(z)).lt.1d-1)) then
      rez=(1d0/2027025d0)*z**7 + (1d0/68918850d0)*z**9 + (1d0/5237832600d0)*z**11&
              +(1d0/659966907600d0)*z**13
elseif ((n.eq.8).and. (abs(realpart(z)).lt.1d-1)) then
      rez=(1d0/34459425d0)*z**8 + (1d0/1309458150d0)*z**10 + (1d0/109994484600d0 )*z**12&
              +(1d0/15179238874800d0)*z**14
elseif ((n.eq.9).and. (abs(realpart(z)).lt.2d-1)) then
      rez=(1d0/654729075d0)*z**9 + (1d0/27498621150d0)*z**11 + (1d0/2529873145800d0)*z**13&
              +(1d0/379480971870000d0)*z**15
elseif ((n.eq.10).and. (abs(realpart(z)).lt.3d-1)) then
      rez=(1d0/13749310575d0 )*z**10 + (1d0/632468286450d0 )*z**12 + (1d0/63246828645000d0 )*z**14&
              +(1d0/10245986240490000d0 )*z**16
elseif ((n.eq.11).and. (abs(realpart(z)).lt.4d-1)) then
      rez=(1d0/316234143225d0)*z**11 + (1d0/15811707161250d0 )*z**13 + (1d0/1707664373415000d0 )*z**15&
              +(1d0/297133600974210000d0 )*z**17
elseif ((n.eq.12).and. (abs(realpart(z)).lt.6d-1)) then
      rez=(1d0/ 7905853580625d0)*z**12 + (1d0/ 426916093353750d0)*z**14 + (1d0/ 49522266829035000d0)*z**16&
              +(1d0/ 9211141630200510000d0)*z**18
elseif ((n.eq.13).and. (abs(realpart(z)).lt.7d-1)) then
      rez=(1d0/213458046676875d0 )*z**13 + (1d0/12380566707258750d0 )*z**15 + (1d0/1535190271700085000d0 )*z**17&
              +(1d0/303967673796616830000d0 )*z**19
elseif ((n.eq.14).and. (abs(realpart(z)).lt.9d-1)) then
      rez=(1d0/6190283353629375d0 )*z**14 + (1d0/ 383797567925021250d0)*z**16 + (1d0/ 50661278966102805000d0)*z**18&
              +(1d0/ 10638868582881589050000d0)*z**20
elseif ((n.eq.15).and. (abs(realpart(z)).lt.1.1d0)) then
      rez=(1d0/191898783962510625d0)*z**15 + (1d0/12665319741525701250d0)*z**17 + (1d0/1773144763813598175000d0)*z**19&
              +(1d0/393638137566618794850000d0 )*z**21
elseif ((n.eq.16).and. (abs(realpart(z)).lt.1.3d0)) then
      rez=(1d0/6332659870762850625d0 )*z**16 + (1d0/443286190953399543750d0 )*z**18 + (1d0/65606356261103132475000d0 )*z**20&
              +(1d0/15351887365098132999150000d0 )*z**22
elseif ((n.eq.17).and. (absz.lt.1.7d0)) then
      rez=(1d0/221643095476699771875d0 )*z**17 + (1d0/16401589065275783118750d0 )*z**19 + (1d0/2558647894183022166525000d0 )*z**21&
              +(1d0/629427381969023452965150000d0 )*z**23
          
elseif ((n.eq.18).and. (absz.lt.2.1d0)) then
      rez=(1d0/8200794532637891559375d0 )*z**18 + (1d0/639661973545755541631250d0 )*z**20 + (1d0/104904563661503908827525000d0 )*z**22&
              +(1d0/27065377424668008477501450000d0 )*z**24
          


else

if (.false.)then !double


  if (n.eq.0) then
  rez=sinh(z)*z**(-1)
  else
  g(0)=z**(-1)
  g(1)=-z**(-2)


  do i=2,n
    g(i)=g(i-2)-dble(2*i-1)*g(0)*g(i-1)
  enddo
  do i=-1,-n-1,-1
    g(i)=g(i+2)+dble(2*i+3)*g(0)*g(i+1)
  enddo

  rez=g(n)*sinh(z)+g(-n-1)*cosh(z)
  endif

else !quad
  y=cmplx(realpart(z),imagpart(z),16)

  if (n.eq.0) then
  rez16=sinh(y)*y**(-1)
  else
  g16(0)=y**(-1)
  g16(1)=-y**(-2)


  do i=2,n
    g16(i)=g16(i-2)-real(2*i-1,16)*g16(0)*g16(i-1)
  enddo
  do i=-1,-n-1,-1
    g16(i)=g16(i+2)+real(2*i+3,16)*g16(0)*g16(i+1)
  enddo

  rez16=g16(n)*sinh(y)+g16(-n-1)*cosh(y)
  endif
rez=cmplx(realpart(rez16),imagpart(rez16),8)


endif !end quad

endif
end subroutine 

subroutine msbesselkc (n, z, rez)
  Implicit None
 Integer, Intent (In) :: n
 Complex (8), Intent (In) :: z
 Complex (8), Intent (Out) :: rez
 Integer :: i
 Complex (8) :: f(0:n)

 Complex (8) :: rezi,rezk

Complex(16) :: f16(0:n)
complex(16)  :: rez16
complex(16)  :: y



if ((n.eq.1) .and. (abs(realpart(z)).lt.1d-3).and.(abs(imagpart(z)).lt.1d-3))  then
rez=1d0/z**2 - (1d0/2d0)+ (1d0/3d0)*z - (1d0/8d0)*z**2 +  (1d0/30d0)*z**3&
        -(1d0/144d0)*z**4 +(1d0/840d0)*z**5 - (1d0/5760d0)*z**6
   !(1d0/45360d0)*z**7 -(1d0/403200d0)*z**8      
else

if (.false.)then !double


 if (n.eq.0) then
 rez=exp(-z)*z**(-1)
 else
 f(0)=z**(-1)
 f(1)=(z+dble(1))*z**(-2)


 do i=2,n
   f(i)=f(i-2)+dble(2*i-1)*f(0)*f(i-1)
 enddo

 rez=exp(-z)*f(n)
 endif

else !!quad

y=cmplx(realpart(z),imagpart(z),16)


 if (n.eq.0) then
 rez16=exp(-y)*y**(-1)
 else
 f16(0)=y**(-1)
 f16(1)=(y+1q0)*y**(-2)


 do i=2,n
   f16(i)=f16(i-2)+real(2*i-1,16)*f16(0)*f16(i-1)
 enddo

 rez16=exp(-y)*f16(n)
 endif
rez=cmplx(realpart(rez16),imagpart(rez16),8)

endif !!end quad
endif

end subroutine 


end module
