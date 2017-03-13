module m_writebevec

  use modinput
  use modmpi
  use mod_constants, only: h2ev
  use m_getunit
  use m_genfilname
  use m_putgetexcitons

  implicit none

  contains 

    subroutine b_writebevec

      integer(4) :: iqmt, alpha, lambda, nvmax, ncmax, icmin, icmax, ivmin, ivmax 
      real(8) :: rbevec, abevec, en1, en2
      real(8), parameter :: epslat = 1.0d-6
      real(8), allocatable :: rbevec_ksum(:,:), abevec_ksum(:,:)
      integer(4) :: i, i1, i2, iv, ic, iknr, ikpnr, ivec(3), rcount, acount
      integer(4) :: un
      real(8) :: abscutoffres(2), abscutoffares(2)
      real(8) :: vklv(3), vklc(3)
      real(8), allocatable :: absvec(:)
      integer(4), allocatable :: idxsort(:), idxsort_desc(:)

      character(256) :: fname, lambdastring
      character(256) :: syscommand, bevecdir, bevecksumdir, excitonevecdir
      character(256) :: tdastring, bsetypestring, tistring, scrtypestring

      if(rank == 0) then 

        bevecdir='BEVEC'
        syscommand = 'test ! -e '//trim(adjustl(bevecdir))//' && mkdir '//trim(adjustl(bevecdir))
        call system(trim(adjustl(syscommand)))
        bevecksumdir='BEVEC_KSUM'
        syscommand = 'test ! -e '//trim(adjustl(bevecksumdir))//' && mkdir '//trim(adjustl(bevecksumdir))
        call system(trim(adjustl(syscommand)))
        excitonevecdir='EXCITON_EVEC'
        syscommand = 'test ! -e '//trim(adjustl(excitonevecdir))//' && mkdir '//trim(adjustl(excitonevecdir))
        call system(trim(adjustl(syscommand)))

        ! Set defaults if writeexcitons is not specified
        if( .not. associated(input%xs%writeexcitons)) then
          input%xs%writeexcitons => getstructwriteexcitons(emptynode)
        end if

        !====================================================!
        ! Read in data to putgetexcitons module              !
        !====================================================!
        ! Requested iqmt
        iqmt = input%xs%bse%iqmt
        if(iqmt == -1) then
          write(*,*) "Warning(b_writebvec): iqmt=-1, setting it to 1"
          iqmt = 1
        end if
        ! Requested excition index range
        if(input%xs%writeexcitons%selectenergy) then 
          en1=input%xs%writeexcitons%minenergyexcitons
          en2=input%xs%writeexcitons%maxenergyexcitons
          if(input%xs%storeexcitons%useev) then 
            en1=en1/h2ev
            en2=en2/h2ev
          end if
          call get_excitons(iqmt=iqmt, e1=en1, e2=en2)
        else
          i1 = input%xs%writeexcitons%minnumberexcitons
          i2 = input%xs%writeexcitons%maxnumberexcitons
          call get_excitons(iqmt=iqmt, a1=i1, a2=i2)
        end if
        !====================================================!

        !====================================================!
        ! Selective write of exciton coefficients to file.   !
        !====================================================!
        if(fcoup_) then
          tdastring=''
        else
          tdastring="-TDA"
        end if
        if(fti_) then 
          tistring="-TI"
        else
          tistring=''
        end if
        bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)//trim(tistring)
        scrtypestring = '-'//trim(input%xs%screening%screentype)

        abscutoffres(1:2) = input%xs%writeexcitons%abscutres
        abscutoffares(1:2) = input%xs%writeexcitons%abscutares

        allocate(absvec(hamsize_))
        allocate(idxsort(hamsize_), idxsort_desc(hamsize_))

        do lambda = iex1_, iex2_

          call genfilname(dirname=trim(excitonevecdir), basename="EXEVEC",&
            & lambda=lambda, iqmt=iq_,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, filnam=fname)

          call getunit(un)

          open(unit=un, file=trim(fname), form='formatted', action='write')

          write(un,'("#",1x,"BSE eigenvector")')
          if(fti_) then 
            write(un,'("#",1x,"Time reversal symmetry was used")') 
          end if
          write(un,'("#")')
          write(un,'("#",1x,"Momentum transfer q_mt:",3f12.7)') vqlmt_
          write(un,'("#",1x,"Size of RR part of the Hamiltonian:",i8)') hamsize_
          write(un,'("#",1x,"Lowest (partially) unoccupied state:",i8)') iuref_
          write(un,'("#")')
          write(un,'("#",1x,"Coefficients are sorted according to thier absolute value")')
          write(un,'("#",2x,"Asolute value cutoff for resonant coefficients:",2E10.3)') abscutoffres
          write(un,'("#",2x,"Asolute value cutoff for anti-resonant coefficients:",2E10.3)') abscutoffares
          write(un,'("#")')
          write(un,'("# Eigenvector number:",i6," with energy/eV: ", f12.7)'), lambda, evals_(lambda)*h2ev
          
          absvec = abs(rvec_(1:hamsize_,lambda))
          call sortidx(hamsize_, absvec, idxsort)
          do i = 1, hamsize_ 
            idxsort_desc(i) = idxsort(hamsize_-i+1)
          end do

          write(un,'("# Resonant contribution")')
          write(un,'("#",1x,a,3x,a,7x,a,7x,a,3x,a,33x,a,33x,a,19x,a,15x,a)')&
            & "alpha","lambda", "ic", "iv", "k_c", "k_v",&
            & "|evec|", "real(evec)", "im(evec)"
          rcount=0
          do i = 1, hamsize_
            alpha = idxsort_desc(i)
            ic = smap_(1,alpha)
            iv = smap_(2,alpha)
            iknr = smap_(3,alpha)
            vklv = vkl0_(1:3,iknr)
            vklc = vkl_(1:3,iknr)
            if(absvec(alpha) > abscutoffres(1) .and. absvec(alpha) < abscutoffres(2)) then
              write(un,'(i7,3(2x,i7),2(3f12.7),3(2x,E23.16))')&
                & alpha, lambda, ic, iv, vklc, vklv,&
                & absvec(alpha), dble(rvec_(alpha,lambda)), aimag(rvec_(alpha, lambda))
              rcount = rcount + 1
            end if
          end do
          write(un,'("# rcount=",i8)'), rcount

          if(fcoup_) then 
            write(un,'("# Anti-resonant contribution")')
            write(un,'("#",1x,a,3x,a,7x,a,7x,a,3x,a,33x,a,33x,a,19x,a,15x,a)')&
              & "alpha","lambda", "ic", "iv", "k_c", "k_v",&
              & "|evec|", "real(evec)", "im(evec)"

            absvec = abs(avec_(1:hamsize_, lambda))
            call sortidx(hamsize_, absvec, idxsort)
            do i = 1, hamsize_ 
              idxsort_desc(i) = idxsort(hamsize_-i+1)
            end do

            acount=0
            do i = 1, hamsize_
              alpha = idxsort_desc(i)
              ic = smap_(1,alpha)
              iv = smap_(2,alpha)
              iknr = smap_(3,alpha)
              vklv = vkl0_(1:3,iknr)
              vklc = vkl_(1:3,iknr)
              if(fti_) then 
                vklv = -vklv
                vklc = -vklc
                call r3frac(epslat, vklv, ivec)
                call r3frac(epslat, vklc, ivec)
              end if
              if(absvec(alpha) > abscutoffares(1) .and. absvec(alpha) < abscutoffares(2)) then
                write(un,'(i7,3(2x,i7),2(3f12.7),3(2x,E23.16))')&
                  & alpha+hamsize_, lambda, ic, iv, vklc, vklv,&
                  & absvec(alpha), dble(avec_(alpha,lambda)), aimag(avec_(alpha, lambda))
                acount=acount+1
              end if
            end do
            write(un,'("# acount=",i8)'), acount

          end if

          close(un)

        end do

        deallocate(idxsort, idxsort_desc)
        deallocate(absvec)
        !====================================================!

        !====================================================!
        ! Write out squared modulus of the resonant          !
        ! (anti-resonant) exciton coefficients to text file. !
        !====================================================!

        nvmax = maxval(koulims_(4,:)-koulims_(3,:)+1)
        ncmax = maxval(koulims_(2,:)-koulims_(1,:)+1)
        icmin = minval(koulims_(1,:))
        icmax = maxval(koulims_(2,:))
        ivmin = minval(koulims_(3,:))
        ivmax = maxval(koulims_(4,:))

        ! Loop over eigenvectors
        lambdaloop: do lambda = iex1_, iex2_

          call getunit(un)

          write(lambdastring, '("_LAMBDA",i4.4)') lambda
          fname=trim('BEVEC'//trim(lambdastring)//'.OUT')
          fname=trim(bevecdir)//'/'//trim(fname)

          open(unit=un, file=trim(fname), form='formatted', action='write')

          ! Loop over transitions
          do alpha = 1, hamsize_

            ! Get absolute indices form combinded index
            ic = smap_(1,alpha)
            iv = smap_(2,alpha)
            iknr = smap_(3,alpha)

            ikpnr = ikmapikq_(iknr)

            ! Resonant
            rbevec = rvec_(alpha, lambda) * conjg(rvec_(alpha, lambda))
            if(fcoup_) then 
              ! Anti-resonant
              abevec = avec_(alpha, lambda) * conjg(avec_(alpha, lambda))
              write(un, '(2(i8, 3E14.7), 2i8, 2E14.7)')&
                & iknr, vkl0_(:, iknr), ikpnr, vkl_(:,ikpnr),&
                & iv, ic, rbevec, abevec
            else
              write(un, '(2(i8, 3E14.7), 2i8, E14.7)')&
                & iknr, vkl0_(:, iknr), ikpnr, vkl_(:,ikpnr),&
                & iv, ic, rbevec
            end if

          end do

          write(un,*)
          if(fcoup_) then 
            write(un, '("k-point, k-point coordinates, kp-point, kp-point coordinates,&
              & valence band, conduction band, abs(X_alpha)^2, abs(Y_alpha)^2 ")')
          else
            write(un, '("k-point, k-point coordinates, kp-point, kp-point coordinates,&
              & valence band, conduction band, abs(X_alpha)^2 ")')
          end if
          write(un,*)
          write(un, '(i8, " : nr. k-points (total)")') nk_max_
          write(un, '(i8, " : nr. k-points (participating)")') nk_bse_
          write(un, '(i8, " : first (partial) unoccupied band")') iuref_
          write(un, '(i8, " : nr. valence states (maximal) ")') nvmax 
          write(un, '(i8, " : nr. conduction states (maximal)")') ncmax
          write(un,*)

          close (un)

        end do lambdaloop
        !====================================================!

        !====================================================!
        ! Write out squared modulus of the resonant          !
        ! (anti-resonant) exciton coefficients summed over   !
        ! all k-points to text file.                         ! 
        !====================================================!

        allocate(rbevec_ksum(icmin:icmax, ivmin:ivmax))
        if(fcoup_) then 
          allocate(abevec_ksum(icmin:icmax, ivmin:ivmax))
        end if

        do lambda = iex1_, iex2_

          rbevec_ksum=0.0d0
          if(fcoup_) abevec_ksum=0.0d0

          do alpha=1, hamsize_
            ic = smap_(1,alpha)
            iv = smap_(2,alpha)
            rbevec_ksum(ic, iv) = rbevec_ksum(ic, iv) +&
              & rvec_(alpha, lambda)*conjg(rvec_(alpha, lambda))
            if(fcoup_) then 
              abevec_ksum(ic, iv) = abevec_ksum(ic, iv) +&
                & avec_(alpha, lambda)*conjg(avec_(alpha, lambda))
            end if
          end do

          call getunit(un)

          write(lambdastring, '("_LAMBDA",i4.4)') lambda
          fname=trim('BEVEC_KSUM'//trim(lambdastring)//'.OUT')
          fname=trim(bevecksumdir)//'/'//trim(fname)

          open(unit=un, file=trim(fname), form='formatted', action='write')

          do iv=ivmin, ivmax
            do ic=icmin, icmax
              if(fcoup_) then 
                write(un, '(2I8, 2E14.7)') iv, ic, rbevec_ksum(ic, iv), abevec_ksum(ic, iv)
              else
                write(un, '(2I8, E14.7)') iv, ic, rbevec_ksum(ic, iv)
              end if
            end do
          end do

          write(un,*)
          if(fcoup_) then 
            write(un, '("valence band, conduction band, sum_k(abs(X)^2), sum_k(abs(Y)^2)")')
          else
            write(un, '("valence band, conduction band, sum_k(abs(X)^2), sum_k(abs(Y)^2)")')
          end if
          write(un,*)
          write(un, '(i8, " : first (partial) unoccupied band")') iuref_ 
          write(un, '(i8, " : nr. valence states (maximal)")') nvmax
          write(un, '(i8, " : nr. conduction states (maximal)")') ncmax
          write(un,*)

          close(un)

        end do

        deallocate(rbevec_ksum)
        if(fcoup_) then 
          deallocate(abevec_ksum)
        end if
        !====================================================!

        !====================================================!
        ! Legacy output for qmt=0, TDA and fixed bands       !
        !====================================================!
        !---------------------------------------------------------------------------
        ! din: new output file for the bandstructure to be able to post-process it
        !---------------------------------------------------------------------------
        if(iq_ ==1 .and. fesel_ == .false.) then

          ! Write out resonant coefficients

          ! Loop over eigenvectors
          do lambda = iex1_, iex2_

            call getunit(un)

            write(fname, '("exciton_evec_res_",i4.4,".dat")') lambda
            fname=trim(excitonevecdir)//'/'//trim(fname)

            open(unit=un, file=trim(fname), form='formatted', action='write')

            ! nkpt total, nv, iv0, nc, ic0
            write(un,*) "# ", nk_bse_, &
            &                 nvmax, ioref_,  &
            &                 ncmax, iuref_

            ! Loop over transitions
            do alpha = 1, hamsize_

              ! Get absolute indices form combinded index
              ic = smap_(1,alpha)
              iv = smap_(2,alpha)
              iknr = smap_(3,alpha)
              vklv(1:3) = vkl0_(:,iknr)

              ! Resonant
              rbevec = rvec_(alpha, lambda) * conjg(rvec_(alpha, lambda))
              write(un, '(i8, 4x, 3f10.6, 2i8, g18.10)')&
                & iknr, vklv, iv, ic, rbevec

            end do

            close (un)

          end do
          
          if(fcoup_) then 

            ! Write out anti-resonant coefficients

            ! Loop over eigenvectors
            do lambda = iex1_, iex2_

              call getunit(un)

              write(fname, '("exciton_evec_ares_",i4.4,".dat")') lambda
              fname=trim(excitonevecdir)//'/'//trim(fname)

              open(unit=un, file=trim(fname), form='formatted', action='write')

              ! nkpt total, nv, iv0, nc, ic0
              write(un,*) "# ", nk_bse_, &
              &                 nvmax, ioref_,  &
              &                 ncmax, iuref_

              ! Loop over transitions
              do alpha = 1, hamsize_

                ! Get absolute indices form combinded index
                ic = smap_(1,alpha)
                iv = smap_(2,alpha)
                iknr = smap_(3,alpha)
                vklv(1:3) = vkl0_(:,iknr)
                ! vklv -> -vklv
                if(fti_) then 
                  vklv = -vklv
                  call r3frac(epslat, vklv, ivec)
                end if

                ! Anti-resonant
                abevec = avec_(alpha, lambda) * conjg(avec_(alpha, lambda))
                write(un, '(i8, 4x, 3f10.6, 2i8, g18.10)')&
                  & iknr, vklv, iv, ic, abevec

              end do

              close (un)

            end do

          end if

        end if
        !====================================================!

        ! Clear read in arrays
        call clear_excitons()

        call barrier

      else

        write(*,*) "m_writebevec: rank", rank, " waiting..."
        call barrier

      end if

    end subroutine b_writebevec

end module m_writebevec