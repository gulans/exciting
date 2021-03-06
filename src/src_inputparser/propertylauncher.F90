
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine propertylauncher
      Use modinput
      Use inputdom
      Use modmain, Only: task
      Use modmpi, Only: rank
      Use spintexture, Only: calculate_spintexture
      Use mod_hdf5, only: fhdf5, hdf5_initialize, hdf5_create_file, hdf5_finalize
      Implicit None
      integer :: l, a, b, c, i
      integer:: k 
#ifdef _HDF5_
  if(rank == 0) then
    call hdf5_initialize()
    fhdf5='property.h5'
    call hdf5_create_file(fhdf5)
  end if
#endif
      call delete_warnings

! properties which depend on the ground state only

       !--------------------------------------------------------
      if( associated( input%properties%wannier)) then
        if( input%properties%wannier%do .ne. "skip") then
          call rereadinput
          call wannierlauncher
        end if
      end if

      !--------------------------------------------------------
      If (associated(input%properties%chargedensityplot)) Then
         call rereadinput
         call rhoplot
      End If

      !--------------------------------------------------------
      if (associated(input%properties%exccplot)) then
         call rereadinput
         call potplot
      end if

      !--------------------------------------------------------
      If (associated(input%properties%wfplot)) Then
         call rereadinput
         ! kstlist should only contain one k-point and state for wave-function plot
         if (size(input%properties%wfplot%kstlist%pointstatepair,2)<1) then
            write(*,*)
            write(*,'("Error(wfplot): /input/properties/wfplot/kstlist must contain")')
            write(*,'(" at least one pointstatepair, but ",i6," were defined")') &
            &  size(input%properties%wfplot%kstlist%pointstatepair,2)
            write(*,*)
            stop
         end if
         select case(input%properties%wfplot%version)
            case('old')
               call wfplot(.false.)
            case('new')
               do i = 1, size(input%properties%wfplot%kstlist%pointstatepair,2)
                  call wfplot_new(input%properties%wfplot%kstlist%pointstatepair(1,i), &
                  &               input%properties%wfplot%kstlist%pointstatepair(2,i))
               end do
            case default
               write(*,*) "Error(propertylauncher): Wrong version! Supported only 'old' and 'new'."
               stop
         end select
      End If

      !--------------------------------------------------------
      If (associated(input%properties%stm).and.(rank==0)) Then
#define STM .true.
         call rereadinput
         Call stm
      End If

      !--------------------------------------------------------
      If (associated(input%properties%LSJ).and.(rank==0)) Then
         call rereadinput
         ! read in input again to reset the magnetic moments for proper symmetry after
         ! a possible run of the groundstate
         task=15
         if (associated(input%properties%LSJ%kstlist)) task=16
         Call writelsj
      End If

      !--------------------------------------------------------
      If (associated(input%properties%TSvdW)) Then
         If ( input%groundstate%do .Eq. "skip" ) Then
            Call init0
            ! read density from file
            Call readstate
         End If
         Call TS_vdW_energy
      End If

      If (associated(input%properties%DFTD2)) Then
         If ( input%groundstate%do .Eq. "skip" ) Then
            Call init0
         End If
         Call DFT_D2_energy
      End If

      If (associated(input%properties%elfplot)) Then
         call rereadinput
         call elfplot
      End If

      !--------------------------------------------------------
      If (associated(input%properties%xcmvecfield)) Then
         If (associated(input%properties%xcmvecfield%plot2d)) Then
            call rereadinput
            task = 82
            Call vecplot
         End If
         If (associated(input%properties%xcmvecfield%plot3d)) Then
            call rereadinput
            task = 83
            Call vecplot
         End If
      End If

      If (associated(input%properties%mvecfield)) Then
         If (associated(input%properties%mvecfield%plot2d)) Then
            call rereadinput
            task = 72
            Call vecplot
         End If
         If (associated(input%properties%mvecfield%plot3d)) Then
            call rereadinput
            task = 73
            Call vecplot
         End If
      End If

      If (associated(input%properties%electricfield)) Then

         If (associated(input%properties%electricfield%plot2d)) Then
            call rereadinput
            task = 142
            Call vecplot
         End If
         If (associated(input%properties%electricfield%plot3d)) Then
            call rereadinput
            task = 143
            Call vecplot
         End If
      End If

      If (associated(input%properties%gradmvecfield)) Then
         call rereadinput
         Call dbxcplot
      End If

      If (associated(input%properties%EFG).and.(rank==0)) Then
         call rereadinput
         Call writeefg
      End If

      If (associated(input%properties%mossbauer).and.(rank==0)) Then
         call rereadinput
         Call mossbauer
      End If

      If (associated(input%properties%expiqr).and.(rank==0)) Then
         call rereadinput
         Call writeexpiqr
      End If

! calculate and print the momentum matrix elements
      If (associated(input%properties%momentummatrix)) Then
         call rereadinput
         Call writepmat
      End If

! IP-RPA dielectric tensor
      If (associated(input%properties%dielmat)) Then
         call rereadinput
         call dielmat
      End If

! BoltzEqu
      If (associated(input%properties%boltzequ)) Then
         call rereadinput
         Call boltzequ
      End If

! Core Overlap
      If (associated(input%properties%coreoverlap)) Then
         call rereadinput
         Call core_overlap
      End If

! MOKE effect
      If (associated(input%properties%moke)) Then
         call rereadinput
         call moke
      End If

! Nonlinear optics: Second Harmonic Generation
      If (associated(input%properties%shg)) Then
         call rereadinput
         do l = 1, size(input%properties%shg%chicomp,2)
           a = input%properties%shg%chicomp(1,l)
           b = input%properties%shg%chicomp(2,l)
           c = input%properties%shg%chicomp(3,l)
           call shg(a,b,c)
         end do
      End If

      If (associated(input%properties%elnes)) Then
         call rereadinput
         Call elnes
      End If

      If (associated(input%properties%eliashberg) .and. (rank==0)) Then
         ! ---- dasabled by DIN ---
         write(*,*)
         write(*,*) " The option is currently disabled!"
         write(*,*)
         stop
         !-------------------
         call rereadinput
         ! this task depends on the results triggered by
         ! "input%properties%phonon"
         ! set the default values if dos element not present
         if (.not.associated(input%properties%dos)) &
           input%properties%dos => getstructdos (emptynode)
         ! electron-phonon coupling
         task=240
         call epcouple
         ! phonon linewidths
         task=245
         call phlwidth
         ! Eliashberg function
         task=250
         Call alpha2f
      End If

      ! Raman scattering
      ! the subroutine raman triggers a phonon calculation, if requested there, and
      ! requires the input of element xs
      If (associated(input%properties%raman)) Then
         call raman
      End If

      If (associated(input%properties%bandstructure)) Then
         call rereadinput
         if( .not. input%properties%bandstructure%wannier) then
           ! tasks are: 20, 21
           task = 20
           Call bandstr
         end if
      End If

      If (associated(input%properties%dos)) Then
         call rereadinput
         if( .not. input%properties%dos%wannier) then
           task = 10
           Call dos
         end if
      End If

      If (associated(input%properties%fermisurfaceplot)) Then
         call rereadinput
         if (associated(input%properties%fermisurfaceplot%plot2d)) then
            task = 101
            Call fermisurf
         Else
            task = 100
            Call fermisurf
         End If
      End If

      If (associated(input%properties%masstensor).and.(rank==0)) Then
         call rereadinput
         task = 25
         Call effmass
      End If

      if( associated( input%properties%polarization)) then
         call rereadinput
         task = 1
         call polarization
      End If

      if (associated(input%properties%ldos)) then
         call rereadinput
         task = 1
         call ldos()
      end if
      !-------------------
      if (associated(input%properties%spintext)) then
         call rereadinput
         task = 42
         call calculate_spintexture(input%properties%spintext%bands(1), input%properties%spintext%bands(2))
      end if

! Finalization of hdf5 output
#ifdef _HDF5_
    if (rank == 0) then
      call hdf5_finalize()
    end if
#endif
End Subroutine propertylauncher
