

Subroutine loaddm
    use mod_potential_and_density, only: mt_dm,dm_copy,DMNullify,DMInitAll
    use mod_atoms, only: natmtot
    implicit none
    complex(8) :: zt1
    integer :: ias,wfi1,wfi2,wfsize


call DMNullify( mt_dm )
call DMInitAll( mt_dm )


open (11, file = "DM.OUT", status = 'old')

read(11,*)wfsize
read(11,*)mt_dm%maxnlo
read(11,*)mt_dm%losize
read(11,*)mt_dm%maxaa
if (.not.allocated(dm_copy)) then
    allocate(dm_copy (wfsize,wfsize,natmtot) )
endif
do ias=1, natmtot
      Do wfi1 = 1, wfsize
            Do wfi2 = 1, wfsize
                  read(11,*)zt1
                  dm_copy(wfi1,wfi2,ias)=zt1
            enddo
      enddo
enddo
close(11)

end Subroutine

! Subroutine loaddm
!     use mod_potential_and_density
!     use mod_kpoint, only: vkl, nkpt
!     use modinput, only: input
!     use mod_eigensystem, only: nmatmax
!     use mod_eigenvalue_occupancy, only: nstfv, nstsv
!     use mod_spin, only: ncmag, nspnfv
!     use precision, only: dp
!     use mod_Gkvector, only: vgkl
!     implicit none
!     integer :: ik
!     complex(dp), allocatable :: evecfv(:, :, :), evecsv(:, :)
!     integer :: is,ia,ias,ir

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !  Generate Density matrix !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     allocate( evecfv(nmatmax, nstfv, nspnfv) )
!     allocate( evecsv(nstsv, nstsv) )

!     if ( input%groundstate%useDensityMatrix ) then
!       call DMNullify( mt_dm )
!       call DMInitAll( mt_dm )
!       mt_dm%alpha%ff = 0._dp
!       if ( associated(input%groundstate%spin) ) then
!         mt_dm%beta%ff = 0._dp
!         if ( ncmag ) then
!           mt_dm%ab%ff = 0._dp
!         end if
!       end if
!       mt_dm%main%ff => mt_dm%alpha%ff
!     end if
!     !call distribute_loop( mpi_env_k, nkpt, firstk, lastk )
!     do ik = 1, nkpt

!       ! get the eigenvectors from file
!       call Getevecfv( vkl(:, ik), vgkl(:, :, :, ik), evecfv )
!       call Getevecsv( vkl(:, ik), evecsv )


!       ! if ( input%groundstate%useDensityMatrix ) then
!       !   ! add to the density and magnetisation
!         call Gendmatmt( ik, evecfv, evecsv )
!       ! else
!       !   call Rhovalk (ik, evecfv, evecsv )
!       ! end if
!       ! call Genrhoir (ik, evecfv, evecsv )
      
!       stop
!     end do ! ik
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! END Generate Density matrix !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
! end Subroutine