subroutine b_bsesurvey
  use modxs
  use modbse
  use mod_kpoint
  use m_genfilname
  use mod_symmetry, only: nsymcrys

  integer(4) :: iqmt
  logical :: fcoup, fti
  integer(4) :: nsymcrys_save

  write(*,*) "b_bsesurvey: Hi!"

  ! General init
  call init0
  ! k-grid init
  call init1
  ! Save variables of the unshifted (apart from xs:vkloff) k grid 
  ! to modxs (vkl0, ngk0, ...)
  call xssave0
  ! q-point and qmt-point setup
  !   Init 2 sets up (task 448):
  !   * A list of momentum transfer vectors form the q-point list 
  !     (modxs::vqmtl and mod_qpoint::vql)
  !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
  !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
  !   * G+qmt quantities (modxs)
  !   * The square root of the Coulomb potential for the G+qmt points
  !   * Reads STATE.OUT
  !   * Generates radial functions (mod_APW_LO)
  call init2

  iqmt = input%xs%bse%iqmt
  fcoup = input%xs%bse%coupling
  fti = input%xs%bse%ti

  if(iqmt == -1) then
    write(*,'("Warning(b_bsesurvey): iqmt=-1, setting it to iqmt = 1")')
    iqmt = 1
  end if
  write(unitout, '("Info(b_bsesurvey): Inspecting transitions for iqmt=",i6)') iqmt

  ! Read Fermi energy from file
  ! Use EFERMI_QMT001.OUT
  call genfilname(iqmt=iqmtgamma, setfilext=.true.)
  call readfermi

  ! Set ist* variables and ksgap in modxs using findocclims
  ! This also reads in 
  ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
  ! modxs:evalsv0, modxs:occsv0
  call setranges_modxs(iqmt, fcoup, fti)

  ! WARNING: ONTOP OF GW STILL IS INCONSISTENT, SINCE OCCUPATION SEARCH IS
  ! NOT DONE ON EVALQP.OUT !?!
  ! If on top of GW
  if(associated(input%gw) .and. iqmt==1) then
    ! Save KS eigenvalues to use them later for renormalizing PMAT
    if(allocated(eval0)) deallocate(eval0)
    allocate(eval0(nstsv, nkptnr))
    eval0=evalsv
    ! If scissor correction is presented, one should nullify it
    input%xs%scissor=0.0d0
    ! Read QP Fermi energies and eigenvalues from file
    ! NOTE: QP evals are shifted by -efermi-eferqp with respect to KS evals
    ! NOTE: getevalqp sets mod_symmetry::nsymcrys to 1
    ! NOTE: getevalqp needs the KS eigenvalues as input
    nsymcrys_save = nsymcrys
    call getevalqp(nkptnr,vkl0,evalsv)
    nsymcrys = nsymcrys_save
    ! Set k and k'=k grid eigenvalues to QP energies
    evalsv0=evalsv
    if(mpiglobal%rank==0) then
      write(unitout,'("Info(b_bsesurvey): Quasi particle energies are read from EVALQP.OUT")')
    end if
  else if(associated(input%gw) .and. iqmt /= 1) then 
    if(mpiglobal%rank==0) then 
      write(*,'("Error(b_bsesurvey): BSE+GW only supported for 0 momentum transfer.")')
    end if
    call terminate
  end if
  call barrier

  ! Select relevant transitions for the construction
  ! of the BSE hamiltonian
  ! Also sets nkkp_bse, nk_bse 
  !   Note: Operates on mpiglobal
  call select_transitions(iqmt, serial=.false., dirname="BSESURVEY")

end subroutine b_bsesurvey