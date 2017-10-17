
subroutine calcmwm(nstart, nend, mstart, mend, minm)

    use modmain, only: pi, zone, zzero
    use modgw,   only: vi, kqset, Gamma, singc1, singc2, mbsiz, &
    &                  minmmat, epsilon, epsh, epsw1, epsw2, freq, mwm, fgw
    implicit none

    ! input variables    
    integer(4), intent(in) :: nstart, nend
    integer(4), intent(in) :: mstart, mend
    complex(8), intent(in) :: minm(mbsiz, nstart:nend, mstart:mend)

    ! local variables
    integer(4) :: iom
    integer(4) :: ie1, ie2, nmdim
    real(8)    :: wkq
    real(8)    :: vi4pi, coefs1, coefs2
    complex(8), allocatable :: wm(:)
    complex(8), external    :: zdotu, zdotc
    external zhemm

    vi4pi  = 4.d0*pi*vi
    coefs1 = singc1*sqrt(vi4pi)
    coefs2 = singc2*vi4pi
    wkq    = 1.d0/dble(kqset%nkpt)

    !-------------------------------------------------
    ! calculate \sum_{ij} M^i_{nm}* W^c_{ij} M^j_{nm}
    !-------------------------------------------------
    allocate(wm(mbsiz))
    do iom = 1, freq%nomeg
      do ie2 = mstart, mend
        do ie1 = nstart, nend
          call zhemv( 'u', mbsiz, zone, epsilon(:,:,iom), mbsiz, &
          &           minmmat(:,ie1,ie2), 1, zzero, wm, 1)
          mwm(ie1,ie2,iom) = wkq*zdotc(mbsiz,minmmat(:,ie1,ie2),1,wm,1)
          if ((Gamma).and.(ie1==ie2)) then
            mwm(ie1,ie2,iom) = mwm(ie1,ie2,iom) + &
            &                  coefs2*epsh(iom,1,1) + &
            &                  coefs1*(zdotu(mbsiz,minmmat(:,ie1,ie2),1,epsw2(:,iom,1),1) + &
            &                          zdotc(mbsiz,minmmat(:,ie1,ie2),1,epsw1(:,iom,1),1))
          end if ! singular term
        end do
      end do
    end do ! iom
    deallocate(wm)
      
    return
end subroutine
