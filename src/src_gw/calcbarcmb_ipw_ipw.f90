!
subroutine calcbarcmb_ipw_ipw(iq)
!
    use modinput
    use constants,             only : zzero, zone, pi, twopi
    use modgw,                 only : Gset, Gamma, kqset, Gqset, Gqbarc
    use mod_product_basis,     only : locmatsiz, mpwipw
    use mod_coulomb_potential, only : barc
    use mod_lattice, only: omega
    Use mod_kpoint, only: nkptnr
    implicit none
    ! input variables
    integer, intent(in) :: iq
    ! local variables 
    integer :: npw, ipw, ipw0  ! PW
    integer :: ngq, igq, jgq   ! IPW
    real(8) :: gqvec(3), gqlen
    real(8) :: vc
    real(8) :: r_c !cutofff radius
    complex(8), allocatable :: tmat1(:,:), tmat2(:,:)
    
    ! external routine 
    external zgemm    
    r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0
    npw = Gqbarc%ngk(1,iq)
    ngq = Gqset%ngk(1,iq)
    
    ! local array
    allocate(tmat1(1:ngq,1:npw))
    tmat1(:,:) = zzero
    
    ipw0 = 1
    !if (Gamma) ipw0 = 2
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,gqvec,gqlen,vc,igq)
!$OMP DO
#endif    
    do ipw = ipw0, npw
      gqvec(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      gqlen = gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2)+gqvec(3)*gqvec(3)
      vc = 4.0d0*pi*(1.0d0-cos(dsqrt(gqlen)*r_c))/gqlen
      if (dsqrt(gqlen) < 1.d-8) then
        write(*,*) 'WARNING(calcbarcmb_ipw_ipw.f90): Zero length vector!', ipw
        vc = twopi*r_c**2 ! cutoff correction for |G+q|=0
        
      endif
    
      
      
      do igq = 1, ngq
        tmat1(igq,ipw) = vc*mpwipw(igq,ipw)
      end do
      
    end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
    
    allocate(tmat2(1:ngq,1:ngq))  
    call zgemm( 'n','c',ngq,ngq,npw, &
    &           zone,tmat1,ngq, &
    &           mpwipw,ngq, &
    &           zzero,tmat2,ngq)
    deallocate(tmat1)
      
    do jgq = 1, ngq
      do igq = 1, ngq
        barc(locmatsiz+igq,locmatsiz+jgq) = tmat2(igq,jgq)
      end do 
    end do
    deallocate(tmat2)
    
    !write(*,*) 'IPW-IPW'
    !do igq = 1, ngq, ngq/10
    !do jgq = 1, ngq, ngq/10
    ! write(*,*) igq, jgq, barc(locmatsiz+igq,locmatsiz+jgq)
    !end do
    !end do
    
end subroutine
    
