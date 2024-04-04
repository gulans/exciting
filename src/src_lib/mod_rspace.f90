module modrspace 

  use mod_atoms, only: natmtot
implicit none

integer,allocatable    :: rgrid_nmtpoints(:) !(ias)
integer                :: rgrid_max_nmtpoints
real (8),allocatable   :: rgrid_mt_rv(:,:,:) !(3,ias,ig)
real (8),allocatable   :: rgrid_mt_rabs(:,:) !(ias,ig)
integer,allocatable    :: rgrid_mt_map(:,:)  !(ig,ias)
complex(8),allocatable :: rgrid_zylm(:,:,:)  !(lm,ias,ig)


contains

subroutine generate_rgrid_mt_data(lmax)
  use modinput
  use mod_Gvector, only: ngrtot, ngrid, ngvec, cfunir
  use mod_atoms, only: nspecies, natoms, atposc, idxas
  use mod_muffin_tin, only: rmt
  use constants, only: zzero
  implicit none
  integer, intent(in) :: lmax

  real(8) :: a1(3),a2(3),a3(3),a1_abs,a2_abs,a3_abs
  real(8) :: base_mat(3,3),base_mat_inv(3,3),col(1,3)
  real(8) :: ratom(3),ratom_i(3), rv(3), rv_abs
  integer :: is,ia,ias,ir1,ir2,ir3,i1,i2,i3,ig,lm,igr
  real(8) :: tp(2)
  real(8) :: rgrid_mt_rv_tmp(3,natmtot,int(ngrtot/natmtot))
  integer :: rgrid_mt_map_tmp(int(ngrtot/natmtot),natmtot)

integer :: p1,p2,p3, count_by_cfunir
real(8) :: pvec(3), r_dist(3),rabs
count_by_cfunir=0
  a1=input%structure%crystal%basevect(1, :)/ngrid(1)
  a2=input%structure%crystal%basevect(2, :)/ngrid(2)
  a3=input%structure%crystal%basevect(3, :)/ngrid(3)
  a1_abs=sqrt(a1(1)**2+a1(2)**2+a1(3)**2)
  a2_abs=sqrt(a2(1)**2+a2(2)**2+a2(3)**2)
  a3_abs=sqrt(a3(1)**2+a3(2)**2+a3(3)**2)

  base_mat(1,1:3)=a1
  base_mat(2,1:3)=a2
  base_mat(3,1:3)=a3
  
  Call r3minv (base_mat, base_mat_inv)
  allocate(rgrid_nmtpoints(natmtot))
  rgrid_nmtpoints=0

if (.true.) then
do ir1=0,ngrid(1)-1
  do ir2=0,ngrid(2)-1
    do ir3=0,ngrid(3)-1
      ig = ir3*ngrid(2)*ngrid(1) + ir2*ngrid(1) + ir1 + 1 
      if (cfunir(ig).lt.0.8d0)then !inside of MT
        count_by_cfunir=count_by_cfunir+1
        rv = ir1*a1 + ir2*a2 + ir3*a3 !point coordinates
        do is=1, nspecies
          do ia=1, natoms(is)
            ias=idxas(ia,is)
            ratom=atposc (:, ia, is) 
            do p1=-1,1
              do p2=-1,1
                do p3=-1,1
                  pvec=input%structure%crystal%basevect(1, :)*p1 + input%structure%crystal%basevect(2, :)*p2 + input%structure%crystal%basevect(3, :)*p3
                  r_dist = rv + pvec - ratom
                  rabs=sqrt(r_dist(1)**2+r_dist(2)**2+r_dist(3)**2)
                  if (rabs.le.rmt(is)) then !found which MT the point belongs to
                    rgrid_nmtpoints(ias) = rgrid_nmtpoints(ias) + 1
                    igr=rgrid_nmtpoints(ias)
                    rgrid_mt_map_tmp(igr,ias)=ig
                    rgrid_mt_rv_tmp(:,ias,igr)=r_dist
                    !have to exit some loops to save time
                  endif
                enddo !p3
              enddo !p2
            enddo !p1

          enddo!ia
        enddo!is
      endif ! cfunir
    enddo
  enddo
enddo
!write(*,*)"count_by_cfunir",count_by_cfunir
else 
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      ratom=atposc (:, ia, is) 

      col(1,:)=ratom
      col=matmul(col, base_mat_inv)
      ratom_i(1)=col(1,1)
      ratom_i(2)=col(1,2)
      ratom_i(3)=col(1,3)
      !write(*,*)"atoma poz",col(1,1),col(1,2),col(1,3)
      !write(*,*)"mt izmērs",rmt(is)/a1_abs,rmt(is)/a2_abs,rmt(is)/a3_abs

  
      do ir1=floor(ratom_i(1)-rmt(is)/a1_abs-20),ceiling(ratom_i(1)+rmt(is)/a1_abs+20)
        if (ir1.lt.0) then 
          i1=ngrid(1)+ir1+1
        elseif (ir1.gt.(ngrid(1)-1)) then
          i1=-ngrid(1)+ir1+1
        else
          i1=ir1+1
        endif
        do ir2=floor( ratom_i(2)-rmt(is)/a2_abs-20),ceiling( ratom_i(2)+rmt(is)/a2_abs+20)
          if (ir2.lt.0) then
            i2=ngrid(2)+ir2+1 
          elseif (ir2.gt.(ngrid(2)-1)) then
            i2=-ngrid(2)+ir2+1
          else
            i2=ir2+1
          endif
          do ir3=floor(ratom_i(3)-rmt(is)/a3_abs-20),ceiling(ratom_i(3)+rmt(is)/a3_abs+20)
            if (ir3.lt.0) then
              i3=ngrid(3)+ir3+1
            elseif (ir3.gt.(ngrid(3)-1)) then
              i3=-ngrid(3)+ir3+1
            else 
              i3=ir3+1
            endif
            rv=ir1*a1+ir2*a2+ir3*a3-ratom
            rv_abs=dsqrt(rv(1)**2+rv(2)**2+rv(3)**2)

            if (rv_abs.le.rmt(is)) then
              rgrid_nmtpoints(ias) = rgrid_nmtpoints(ias) + 1
              ig = (i3-1)*ngrid(2)*ngrid(1) + (i2-1)*ngrid(1) + i1
              
              if(ig.gt.ngrtot)then
                write(*,*)"ig,ngrtot"
                write(*,*)ig,ngrtot
                write(*,*)"i1,i2,i3"
                write(*,*)i1,i2,i3
                write(*,*)"ir1,ir2,ir3"
                write(*,*)ir1,ir2,ir3
                write(*,*)"ratom"
                write(*,*)ratom
                stop
              endif

              rgrid_mt_rv_tmp(:,ias,rgrid_nmtpoints(ias))=rv
              rgrid_mt_map_tmp(rgrid_nmtpoints(ias),ias)=ig            
            endif
!write(*,*)ir1,ir2,ir3
          enddo !ir1
        enddo !ir2
      enddo !ir3
    enddo !ia
  enddo !is
endif
  rgrid_max_nmtpoints=maxval(rgrid_nmtpoints(:))
  
  allocate( rgrid_mt_rv(3,natmtot,rgrid_max_nmtpoints) )
  allocate( rgrid_mt_map(rgrid_max_nmtpoints,natmtot) )
  allocate( rgrid_mt_rabs(natmtot,rgrid_max_nmtpoints) )
  allocate( rgrid_zylm((lmax+1)**2,natmtot,rgrid_max_nmtpoints) )
  rgrid_mt_rv=0d0
  rgrid_mt_map=0
  rgrid_mt_rabs=0d0
  rgrid_zylm=zzero

  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
  
      rgrid_mt_rv (:,ias,1:rgrid_nmtpoints(ias)) = rgrid_mt_rv_tmp(:,ias,1:rgrid_nmtpoints(ias))
      rgrid_mt_map(1:rgrid_nmtpoints(ias),ias) = rgrid_mt_map_tmp(1:rgrid_nmtpoints(ias),ias)

      do ig=1, rgrid_nmtpoints(ias)
        rv=rgrid_mt_rv(:,ias,ig)
        rv_abs=dsqrt(rv(1)**2+rv(2)**2+rv(3)**2)
        rgrid_mt_rabs(ias,ig)=rv_abs

        !if (rv_abs.lt.1e-12) then 
        if (rv_abs.eq.0d0) then 
          tp(1)=0d0
        else
          tp(1)=dacos(rv(3)/rv_abs) !theta
        endif

        !if ((abs(rv(1)).lt.1e-12).and.(abs(rv(2)).lt.1e-12))then
        if ((abs(rv(1)).eq.0d0).and.(abs(rv(2)).eq.0d0))then
          tp(2)= 0d0 !to avoid NaN
        else
          tp(2)= sign(1d0,rv(2))*dacos(rv(1)/dsqrt(rv(1)**2+rv(2)**2)) !phi
        endif

        call genylm(lmax, tp, rgrid_zylm(:,ias,ig))
      enddo !ig
    enddo !ia
  enddo !is

  open(11,file='rgrid_info.dat',status='replace')
  write(11,*)"ngrtot",ngrtot
  write(11,*)"ngvec",ngvec
  write(11,*)"total number of MT points:",sum(rgrid_nmtpoints)
  write(11,*)"rgrid_max_nmtpoints:",rgrid_max_nmtpoints
  write(11,*)"ias, rgrid_nmtpoints(ias)"
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia,is)
      write(11,*)ias, rgrid_nmtpoints(ias)
    enddo
  enddo
  close(11)


!stop
end subroutine 


end module
