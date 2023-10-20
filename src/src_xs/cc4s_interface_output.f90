module cc4s_interface_output

    use precision, only: dp
    implicit none


contains 


    !> Computes Coulomb potential for a given vector \(\mathbf{q}\)
    !> according to 
    !> \[
    !>    v(\mathbf{q}) = \frac{4\pi}{|\mathbf{q}|^2}
    !>                                                 \].
    !> Todo (Max): Replace v(0) = 0 with expression from
    !> Spencer & Alavi, PHYSICAL REVIEW B 77, 193110 (2008)
    function coulomb_potential(vector) result(coulomb_pot)
      use constants, only: fourpi
      use math_utils, only: all_zero
      !> Recprocal-space vector in cartesian coordinates
      real(dp), intent(in) :: vector(3) 
      real(dp) :: coulomb_pot

      if(.not. all_zero(vector)) then
        coulomb_pot = fourpi / (norm2(vector)**2._dp)
      else 
        coulomb_pot =  0._dp
      end if
    end function

    !> Computes Coulomb potential for a given vector \(\mathbf{q}\)
    !> according to Spencer & Alavi, PHYSICAL REVIEW B 77, 193110 (2008)
    function coulomb_potential_spencer_alavi(unit_cell_volume, vector) result(coulomb_pot)
      use constants, only: fourpi, twopi
      use math_utils, only: all_zero

      !> Recprocal-space vector in cartesian coordinates
      real(dp), intent(in) :: vector(3) 
      real(dp), intent(in) :: unit_cell_volume 
      real(dp) :: coulomb_pot
      real(dp) :: spencer_alavi_factor 
      real(dp) :: r_c 

      r_c = ( unit_cell_volume * 3._dp / fourpi ) **(1._dp / 3._dp) 
      spencer_alavi_factor =  (1._dp - cos( norm2(vector) * r_c))

      if(.not. all_zero(vector)) then
        coulomb_pot = fourpi / (norm2(vector)**2._dp) * spencer_alavi_factor
      else 
        coulomb_pot = twopi * r_c ** 2._dp
      end if

    end function


    subroutine CC4S_out_coulomb_potential_element(grid_vectors)
            
        !> (G+q)-vectors in cartesian coordinates
        real(dp), dimension(:), intent(in) :: grid_vectors(:,:,:)

        !> Number of q-vectors
        integer :: n_qvecs
        !> Number of (G+q)-vectors per q-vector 
        integer :: n_gqvecs
        !> Running indices q- and (G+q)-vectors
        integer :: iq, igq

        n_gqvecs = size(grid_vectors,dim=2)

        n_qvecs= size(grid_vectors,dim=3)

        open (unit=498, file='CoulombPotential.elements')

        do iq = 1, n_qvecs
        do igq = 1, n_gqvecs
            write(498,'(E22.15)') coulomb_potential(grid_vectors(:,igq, iq))
        end do
        end do

        close (498)

    end subroutine CC4S_out_coulomb_potential_element

   ! Works only for Gamma-only calculations!!!
  subroutine CC4S_out_coulomb_vertex_singular_vectors_element(left_singular_vectors)
    implicit none

    !> (G+q)-vectors in cartesian coordinates
    complex(dp), intent(in) :: left_singular_vectors(:, :)
    integer :: n_qvecs, n_gqvecs_per_qvector, iq, igq, q_last



  end subroutine CC4S_out_coulomb_vertex_singular_vectors_element

  ! Copy from cc4s-aims
  subroutine CC4S_out_coulomb_potential_yaml(grid_vectors, do_gamma_only)
    implicit none

    !> (G+q)-vectors
    real(dp), dimension(:), intent(in) :: grid_vectors(:,:,:)
    integer :: n_qvecs, n_gqvecs_per_qvector

    logical, intent(in) :: do_gamma_only

    n_gqvecs_per_qvector = size(grid_vectors,dim=2)
    if(do_gamma_only) then 
      n_qvecs = 1
    else
      n_qvecs= size(grid_vectors,dim=3)
    end if

    !local


    open (unit=498, file='CoulombPotential.yaml')
    write (498, '(A12)') "version: 100"
    write (498, '(A12)') "type: Tensor"
    write (498, '(A18)') "scalarType: Real64"
    
    write (498, '(A11)') "dimensions:"
    write (498, '(2X,A1,1X, A7,1X,I5)') "-", "length:",  n_qvecs*n_gqvecs_per_qvector
    write (498, '(4X, A14)') "type: Momentum"
    write (498, '(A9)') "elements:"
    write (498, '(2X, A14)') "type: TextFile"
    write (498, '(A25)') "unit: 1.0000   # =Hartree"

    close (498)

  end subroutine CC4S_out_coulomb_potential_yaml

  !*******************************************************

  !> Writes the fil 'GridVectors.yaml' for the cc4s interface.
  !> See https://manuals.cc4s.org/user-manual/objects/GridVectors.html#ID-GridVectors
  subroutine CC4S_out_grid_vectors_yaml(reciprocal_lattice, grid_vectors)
    implicit none

    !> (G+q)-vectors in lattice coordinates 
    real(dp), intent(in) :: grid_vectors(:,:,:)
    real(dp), intent(in) :: reciprocal_lattice(3, 3)

    !local
    integer ::iq, igq, n_qvecs, n_gqvecs, ivec, shape_gq_vectors(3)

    n_gqvecs = size(grid_vectors,dim=2)* size(grid_vectors,dim=2)

    shape_gq_vectors = shape(grid_vectors)
    open (unit=498, file='GridVectors.yaml')
    write (498, '(A12)') "version: 100"
    write (498, '(A12)') "type: Tensor"
    write (498, '(A18)') "scalarType: Real64"
    !Leave out "properties"-section as it is not used anywhere yet
    write (498, '(A11)') "dimensions:"
    write (498, '(2X,A1,1X, A7,1X,I5)') "-", "length:",  size(grid_vectors,dim=1)
    write (498, '(4X, A12)') "type: Vector"
    ! Shape N_Gvecs*N_qvecs
    write (498, '(2X,A1,1X, A7,1X,I10)') "-", "length:",    n_gqvecs
    write (498, '(4X, A14)') "type: Momentum"

    write (498, '(A9)') "elements:"
    write (498, '(4X, A14)') "type: TextFile"
    write (498, '(A5,E22.15)') "unit:", 1._dp ! Atomic units

     write (498, '(A9)') "metaData:"

    ! Write reciprocal lattice

      write (498, '(2X,A3,1X,A1,E22.15,A1,1X,E22.15,A1,1X,E22.15,A1)') "Gi:", "[",&
          reciprocal_lattice(1,1), ',', &
          reciprocal_lattice(2,1), ',',&
          reciprocal_lattice(3,1),']'
      write (498, '(2X,A3,1X,A1,E22.15,A1,1X,E22.15,A1,1X,E22.15,A1)') "Gj:", "[",&
          reciprocal_lattice(1,2), ',', &
          reciprocal_lattice(2,2), ',',&
          reciprocal_lattice(3,2),']'
      write (498, '(2X,A3,1X,A1,E22.15,A1,1X,E22.15,A1,1X,E22.15,A1)') "Gk:", "[",&
          reciprocal_lattice(1,3), ',', &
          reciprocal_lattice(2,3), ',',&
          reciprocal_lattice(3,3),']'

    close (498)

  end subroutine CC4S_out_grid_vectors_yaml

  !> Writes the file 'GridVectors.element' for the cc4s interface.
  !> See https://manuals.cc4s.org/user-manual/objects/GridVectors.html#ID-GridVectors
  subroutine CC4S_out_grid_vectors_element(grid_vectors, do_gamma_only)
    implicit none

    logical, intent(in) :: do_gamma_only


    !> (G+q)-vectors in lattice coordinates
    real(dp), intent(in) :: grid_vectors(:,:,:)

    !local
    integer ::iq, igq, n_qvecs, n_gqvecs_per_qvector, ivec, shape_gq_vectors(3)

    n_gqvecs_per_qvector = size(grid_vectors,dim=2)
    if(do_gamma_only) then 
      n_qvecs = 1
    else
      n_qvecs= size(grid_vectors,dim=3)
    end if
    shape_gq_vectors = shape(grid_vectors)
    open (unit=498, file='GridVectors.element')
   
    do iq = 1, n_qvecs
      do igq = 1, n_gqvecs_per_qvector
      
        write (498, '(3E20.10)') grid_vectors(:, igq, iq)
      end do
    end do

    close (498)

  end subroutine CC4S_out_grid_vectors_element

  ! Copied from cc4s-aims
  ! Assume that full vertex for one (k,q)-pair is to be written,
  ! therefore neglect OAF_subsizes and OAF_starts and replace MPI_write_all
  ! by MPI_write_at
  subroutine CC4S_parallel_writing_cmplx(coulomb_vertex, fh, &
                &  curr_kq_pair, CC4S_n_spin, CC4S_nuf_states)


#ifdef MPI
    use mpi
#endif

    implicit none

    !Copied from cc-aims and modified.
    !PURPOSE: Primarily to be used to write the OAF Coulomb vertex to
    !file in parallel. It assumes that a file object has already been
    !opened via MPI_File_open according to:
    !
    !   call MPI_File_open(mpiglobal%comm, coulomb_vertex_filename_bin, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
    !   &MPI_INFO_NULL, coulomb_vertex_filehandle, ierr)
    ! where the parameters should be set according to
    ! 
    ! character(len=200) :: coulomb_vertex_filename_bin = "CoulombVertex.elements"
    ! 
    complex(dp), intent(in)  :: coulomb_vertex(:, :, :)
    integer, intent(in) :: fh !MPI File handle
    ! Number of basis vectors (= number of G-vectors or number of auxiliary field vectors) used for construction of Coulomb-Vertex
    integer :: n_basis_vectors
    ! Index of (k,q)-pair (=1 for Gamma only calculation)
    integer, intent(in) :: curr_kq_pair
    integer, intent(in) :: CC4S_n_spin
    ! Number of unfrozen states (in exciting always equal to total number of states?)
    integer, intent(in) :: CC4S_nuf_states


    !local
    integer :: ierr
    integer :: subarr
    integer :: nbytes
#ifdef MPI
    integer(kind=MPI_OFFSET_KIND) :: offset

    !CC4S always expects a complex vertex, even if the imaginary part is 0
    ! call MPI_Type_create_subarray(3, OAF_sizes, OAF_subsizes, OAF_starts, &
    ! &MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, subarr, ierr)

    ! call MPI_Type_commit(subarr, ierr)

    ! Number of basis vectors (G-vectors, auxiliary field functions, or mixed basis functions )
    n_basis_vectors = size(coulomb_vertex, dim=1)

    ! TODO (Max): Change this to the proper offset for parallel writing
    !offset = 16*n_basis_vectors*(CC4S_nuf_states**2)*CC4S_n_spin*curr_kq_pair
    offset = 0
    ! call MPI_File_set_view(fh, offset, MPI_DOUBLE_COMPLEX, subarr, 'native', MPI_INFO_NULL, ierr)

    nbytes = size(coulomb_vertex)
    write(*, *) "Write to file Coulomb Vertex of size (MB):", sizeof(coulomb_vertex) / (1024._dp)**2
    ! TODO(Max): Replace by MPI_file_write_at for parallel writing
    !call MPI_File_write_at(fh, offset, coulomb_vertex, nbytes, MPI_DOUBLE_COMPLEX, &
    !       &MPI_STATUS_IGNORE, ierr)
    call MPI_File_write(fh, coulomb_vertex, nbytes, MPI_DOUBLE_COMPLEX, &
            &MPI_STATUS_IGNORE, ierr)
#endif
  end subroutine CC4S_parallel_writing_cmplx

  ! Assume that full vertex for one (k,q)-pair is to be written,
  ! therefore neglect OAF_subsizes and OAF_starts and replace MPI_read_all
  ! by MPI_read_at
  subroutine CC4S_parallel_reading_cmplx(coulomb_vertex, fh, &
                & curr_kq_pair, CC4S_n_spin, CC4S_nuf_states)


#ifdef MPI
    use mpi
#endif

    implicit none

    !PURPOSE: Primarily to be used to write the OAF Coulomb vertex to
    !file in parallel. It assumes that a file object has already been
    !opened via MPI_File_open.

    !input
    ! Assume OAF_size = OAF_subsize because no SCALAPACK involved
    ! integer, dimension(3), intent(in) :: OAF_sizes, OAF_subsizes, OAF_starts
    ! double complex, dimension(OAF_subsizes(1), OAF_subsizes(2), OAF_subsizes(3)), &
    complex(dp), intent(out)  :: coulomb_vertex(:, :, :)
    integer, intent(in) :: fh !MPI File handle
    integer :: n_basis_vectors
    integer, intent(in) :: curr_kq_pair, CC4S_n_spin, CC4S_nuf_states

    !local
    integer :: ierr
    integer :: subarr
    integer :: nbytes
#ifdef MPI
    integer(kind=MPI_OFFSET_KIND) :: offset

    !CC4S always expects a complex vertex, even if the imaginary part is 0
    ! call MPI_Type_create_subarray(3, OAF_sizes, OAF_subsizes, OAF_starts, &
    ! &MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, subarr, ierr)

    ! call MPI_Type_commit(subarr, ierr)

    n_basis_vectors = size(coulomb_vertex, dim=1)

    offset = 16*n_basis_vectors*(CC4S_nuf_states**2)*CC4S_n_spin*curr_kq_pair

    ! call MPI_File_set_view(fh, offset, MPI_DOUBLE_COMPLEX, subarr, 'native', MPI_INFO_NULL, ierr)

    nbytes = size(coulomb_vertex)
    call MPI_File_read_at(fh, offset, coulomb_vertex, nbytes, MPI_DOUBLE_COMPLEX, &
            &MPI_STATUS_IGNORE, ierr)
#endif
  end subroutine CC4S_parallel_reading_cmplx


  ! Copied from cc-aims
  subroutine CC4S_out_eigenenergies_yaml(do_gamma_only, version, e_unit, e_fermi, scf_energies, &
                  &n_k_points, n_spin)
  
          character(len=*), intent(in) :: version
          double precision, intent(in) :: e_unit, e_fermi
          integer, intent(in) :: n_k_points
          integer, intent(in) :: n_spin
          ! Flattened and sorted scf-energies
          double precision, dimension(:) :: scf_energies

          ! True if Gamma-only calculation. In this case, nonZeroCondition is not needed. 
          logical, intent(in) :: do_gamma_only

          !local
          integer :: i_k_point
          integer :: i_spin
          integer :: i_element


          
  
          open(unit=501, file='EigenEnergies.yaml')
          write(501,'(A8, 1X, A5)') "version:", version
          write(501,'(A12)') "type: Tensor"
          write(501,'(A18)') "scalarType: Real64"
          !write(501,'(A8)') "indices:"
          !write(501,'(2X, A8)') "orbital:"
          !if(n_spin .eq. 2) then
          !  write(501,'(4X, A10)') "type: spin"
          !else
          !  write(501,'(4X, A13)') "type: spatial"
          !end if 
          write(501,'(A11)') "dimensions:"
          write(501,'(A9, 1X, I0)') "- length:", size(scf_energies)
          write(501,'(2X, A11)') "type: State"
          !write(501,'(4X, A13)') "type: orbital"
          !write(501,'(A5, 1X, A100)') "data:", scf_energies_filename
          write(501,'(A9)') "elements:"
          write(501,'(2X, A14)') "type: TextFile"
          write(501,'(A5, 1X, F12.7)') "unit:", e_unit
          write(501, '(A9)' ) "metaData:"
          write(501,'(2X, A12, 1X, F25.16)') "fermiEnergy:", e_fermi
          write(501,'(2X, A9)') "energies:"

          do i_element=1,size(scf_energies)
            write(501,'(2X, A1, 1X, F25.16)') "-", scf_energies(i_element)
          end do
          
          ! nonZeroCondition block
          if(.not. do_gamma_only) then
            write(501,'(A17)') "nonZeroCondition:"
            write(501, '(2X, A4)') "all:"
            !Spin
            write(501, '(4X, A2)') "0:"
            write(501, '(6X, A22)') "name: SpinConservation"
            write(501, '(6X, A11)') "properties:"
            write(501, '(8X, A2)') "0:"
            write(501, '(10X, A14)') "property: Spin"
            write(501, '(10X, A12)') "dimension: 0"
            write(501, '(6X, A9)') "nonZeros:"

            do i_spin=0,n_spin-1
                write(501,'(8X, A1, 1X, A1, I0, A1)') "-", "[", i_spin, "]"
            end do

            !crystal momentum
            write(501, '(4X, A2)') "1:"
            write(501, '(6X, A26)') "name: MomentumConservation"
            write(501, '(6X, A11)') "properties:"
            write(501, '(8X, A2)') "0:"
            write(501, '(10X, A25)') "property: CrystalMomentum"
            write(501, '(10X, A12)') "dimension: 0"
            write(501, '(6X, A9)') "nonZeros:"

            do i_k_point=0,n_k_points-1
                write(501,'(8X, A1, 1X, A1, I0, A1)') "-", "[", i_k_point, "]"
            end do

          end if
              
          close(unit=501)
  
  end subroutine CC4S_out_eigenenergies_yaml 

  ! Copied from cc-aims
  subroutine CC4S_output_scf_energies(scf_energies, n_states, n_spin, &
                &n_k_points, lowest_state, bin)

        !PURPOSE: write the eigenvalues, intent(in)  of the underlying SCF-method to file called 'file_name'
        !INPUT: 
        !       *scf_energies: A pointer to a three dimensional array of dimension 
        !                      n_states x n_spin x n_k_points
        !                      which contains the scf-eigenvalues for a generally unrestricted
        !                      (n_spin=2) calculation with an arbitrary k-mesh. n_states denotes
        !                      the number of states per k-point (i.e band-index). The array is 
        !                      expected to not be distributed
        !  
        !       *file_name:    Name of the file the scf-energies will be written to
        !
        !       *lowest_state: For frozen-core calculations lowest_state contains the first unfrozen
        !                      band-index. For full-electron calculations specify lowest_state=1
        !
        !       *bin:          Output file formatted (bin=.false.) or as binary (bin=.true.)

        !input variables
        integer, intent(in) :: n_states
        integer, intent(in) :: n_spin
        integer, intent(in) :: n_k_points 
        integer, intent(in) :: lowest_state
        double precision, dimension(n_states, n_spin, n_k_points), intent(in) :: scf_energies
        character(len=22) :: file_name = "EigenEnergies.elements"
        logical, intent(in) :: bin

        !local variables
        integer :: i_k_point, i_state, i_spin

        if(.not. bin) then
          open(unit=334, file=trim(file_name), action='write')
          do i_k_point=1,n_k_points
            do i_spin=1,n_spin
              do i_state=lowest_state,n_states !        1,n_states-lowest_state+1
                write(334,'(F25.16)') scf_energies(i_state, i_spin, i_k_point)
              end do
            end do
          end do
        else
          open(unit=334, file=trim(file_name), action='write', form='unformatted', access='sequential')
          do i_k_point=1,n_k_points
            do i_spin=1,n_spin
              do i_state=lowest_state, n_states !1,n_states-lowest_state+1
                write(334) scf_energies(i_state, i_spin, i_k_point)
              end do
            end do
          end do

        end if 
        close(unit=334)

end subroutine CC4S_output_scf_energies

!*******************************************************
  subroutine CC4S_out_orbital_properties_yaml(ordered_indices, &
                  &n_k_points, n_states, n_spins)
    implicit none

    integer, dimension(:), intent(in) :: ordered_indices
    integer, intent(in) :: n_k_points, n_states, n_spins

    !local
    integer :: i,j,lb
    integer, dimension(:), allocatable :: unordered_k_indices
    integer, dimension(:), allocatable :: unordered_spin_indices
    
    integer :: CC4S_nuf_states ! Added by Max

    ! TODO (Max): Check if number of unfrozen states is really equal to number of states
    CC4S_nuf_states = n_states

    allocate(unordered_k_indices(n_k_points*CC4S_nuf_states*n_spins))
    allocate(unordered_spin_indices(n_k_points*CC4S_nuf_states*n_spins))


    !Momentum indices
    do i=1,n_k_points
      lb = (i-1)*CC4S_nuf_states*n_spins + 1
      unordered_k_indices(lb:lb+CC4S_nuf_states*n_spins-1) = i
    end do

    !Spin indices
    do i=1,n_k_points
      do j=1,n_spins
        lb = ((i-1)*n_spins + (j-1))*CC4S_nuf_states + 1
        unordered_spin_indices(lb:lb+CC4S_nuf_states-1) = j
      end do
    end do


    open(unit=499, file='State.yaml')
    write(499,'(A12)') "version: 100"
    write(499,'(A20)') "dimensionType: State"
    write(499,'(A11)') "properties:"
    write(499,'(2X, A5)') "Spin:"

    if(n_spins == 2) then
      write(499,'(4X, A7)') "0: +0.5"
      write(499,'(4X, A7)') "1: -0.5"
    else
      write(499,'(4X, A7)') "0: +0.5"
    end if 

    write(499,'(A16)') "propertyIndices:"
    write(499,'(2X, A16)') "CrystalMomentum:"
    do i=1,size(ordered_indices)
      !write(499,'(2X, I0, A1, 1X, I0)') i-1, ":", unordered_k_indices(ordered_indices(i)) - 1
      write(499,'(2X, A1, 1X, I0)') "-", unordered_k_indices(ordered_indices(i)) - 1
    end do

    write(499,'(2X, A5)') "Spin:"
    do i=1,size(ordered_indices)
      !write(499,'(2X, I0, A1, 1X, I0)') i-1, ":", unordered_spin_indices(ordered_indices(i)) - 1
      write(499,'(2X, A1, 1X, I0)') "-", unordered_spin_indices(ordered_indices(i)) - 1
    end do

    close(499)
  end subroutine CC4S_out_orbital_properties_yaml
 
  ! copied from cc-aims
  subroutine CC4S_out_spin_yaml(version, n_states, n_spin, n_k_points, e_unit)
        implicit none
  
        !character(len=*), intent(in) :: datafile_name
        character(len=*), intent(in) :: version
        integer, intent(in) :: n_states
        integer, intent(in) :: n_spin
        integer, intent(in) :: n_k_points
        double precision, intent(in) :: e_unit

        !local
        integer :: i_spin, i_k_point
  
        open(unit=500, file='Spins.yaml')
        write(500,'(A8, 1X, A5)') "version:", version
        write(500,'(A8)') "indices:"
        write(500,'(2X, A8)') "orbital:"
        !Assuming this file will only be written out for spin-polarized case
        write(500,'(4X, A10)') "type: spin"
        write(500,'(A11)') "dimensions:"
        write(500,'(2X, A9, 1X, I0)') "- length:", n_states
        write(500,'(4X, A15)') "type: orbital"
        !write(500,'(A5, 1X, A100)') "data:", datafile_name
        write(500,'(A5, 1X, F12.7)') "unit:", e_unit

        write(500,'(A17)') "nonZeroCondition:"
        write(500, '(2X, A4)') "all:"
        !Spin
        write(500, '(4X, A2)') "0:"
        write(500, '(6X, A22)') "name: SpinConservation"
        write(500, '(6X, A11)') "properties:"
        write(500, '(8X, A2)') "0:"
        write(500, '(10X, A14)') "property: Spin"
        write(500, '(10X, A12)') "dimension: 0"
        write(500, '(6X, A9)') "nonZeros:"

        do i_spin=0,n_spin-1
          write(500,'(8X, A1, 1X, A1, I0, A1)') "-", "[", i_spin, "]"
        end do

        !crystal momentum
        write(500, '(4X, A2)') "1:"
        write(500, '(6X, A26)') "name: MomentumConservation"
        write(500, '(6X, A11)') "properties:"
        write(500, '(8X, A2)') "0:"
        write(500, '(10X, A25)') "property: CrystalMomentum"
        write(500, '(10X, A12)') "dimension: 0"
        write(500, '(6X, A9)') "nonZeros:"

        do i_k_point=0,n_k_points-1
          write(500,'(8X, A1, 1X, A1, I0, A1)') "-", "[", i_k_point, "]"
        end do
        close(unit=500)
  end subroutine CC4S_out_spin_yaml

  ! Copy from cc4s-aims
  subroutine CC4S_out_coulvertex_yaml(do_gamma_only, version, n_states, n_spin, n_basbas_OAF, &
    &momentum_triples_list, e_unit, bin, isreal)
    implicit none

    character(len=*), intent(in) :: version
    integer, intent(in) :: n_states
    integer, intent(in) :: n_spin
    integer, intent(in) :: n_basbas_OAF
    integer, dimension(:, :), intent(in) :: momentum_triples_list
    double precision, intent(in) :: e_unit
    logical, intent(in) :: bin
    ! TODO (Max): exciting is always complex, input argument not needed
    logical, intent(in) :: isreal
    ! True if Gamma-only calculation. In this case, nonZeroCondition is not needed. 
    logical, intent(in) :: do_gamma_only

    !local
    character(len=25) :: grid
    integer :: i, i_spin

    !if(isreal) then
    !  grid='halfGrid'
    !else
    !  grid='fullGrid'
    !end if

    open (unit=503, file='CoulombVertex.yaml')
    write (503, '(A8, 1X, A5)') "version:", version
    write (503, '(A12)') "type: Tensor"
    write (503, '(A21)') "scalarType: Complex64"
    !write(503,'(A8)') "indices:"
    !write(503,'(2X, A9)') "momentum:"
    !write(503,'(4X, A5, 1X, A25)') "type:", grid
    !write(503,'(2X, A8)') "orbital:"
    !if(n_spin .eq. 2) then
    !  write(503,'(4X, A10)') "type: spin"
    !else
    !  write(503,'(4X, A13)') "type: spatial"
    !end if
    write (503, '(A11)') "dimensions:"
    write (503, '(A9, 1X, I0)') "- length:", n_basbas_OAF
    write (503, '(2X, A20)') "type: AuxiliaryField"
    write (503, '(A9, 1X, I0)') "- length:", n_states
    !write(503,'(2X, A13)') "type: orbital"
    write (503, '(2X, A11)') "type: State"
    write (503, '(A9, 1X, I0)') "- length:", n_states
    !write(503,'(2X, A13)') "type: orbital"
    write (503, '(2X, A11)') "type: State"

    write (503, '(A9)') "elements:"
    if (bin) then
    write (503, '(2X, A20)') "type: IeeeBinaryFile"
    else
    write (503, '(2X, A14)') "type: TextFile"
    end if

    !if(bin) then
    !  write(503,'(A9)') "binary: 1"
    !  write(503,'(A5, 1X, A50)') "data:", coulomb_vertex_filename_bin
    !else
    !  write(503,'(A5, 1X, A50)') "data:", coulomb_vertex_filename_dat
    !end if

    write (503, '(A5, 1X, F12.7)') "unit:", e_unit

    write (503, '(A9)') "metaData:"
    if (isreal) then
    write (503, '(2X, A11)') "halfGrid: 1"
    else
    write (503, '(2X, A11)') "halfGrid: 0"
    end if

    ! TODO (Max): Do not write following part for Gamma-only calcs
    if(.not. do_gamma_only) then
      write (503, '(A17)') "nonZeroCondition:"
      write (503, '(2X, A4)') "all:"
      !Spins
      write (503, '(4X, A2)') "0:"
      write (503, '(6X,A22)') "name: SpinConservation"
      write (503, '(6X, A11)') "properties:"
      write (503, '(8X, A2)') "0:"
      write (503, '(10X, A14)') "property: Spin"
      write (503, '(10X, A12)') "dimension: 1"
      write (503, '(8X, A2)') "1:"
      write (503, '(10X, A14)') "property: Spin"
      write (503, '(10X, A12)') "dimension: 2"
      write (503, '(6X, A9)') "nonZeros:"
      do i_spin = 0, n_spin - 1
      write (503, '(8X, A1, 1X, A1, I0, A1, I0, A1)') "-", "[", &
        &i_spin, ",", i_spin, "]"
      end do

      !Crystal-/Transfer-momenta
      !Re-arrange/permute the following 3 dimensions
      !so that the order is 0, 1, 2
      !and accordingly permute the indices within the
      !momentum triples.
      !This will be required in future formats of CC4S
      write (503, '(4X, A2)') "1:"
      write (503, '(6X,A26)') "name: MomentumConservation"
      write (503, '(6X, A11)') "properties:"
      write (503, '(8X, A2)') "0:"
      write (503, '(10X, A25)') "property: CrystalMomentum"
      write (503, '(10X, A12)') "dimension: 0" !"dimension: 1"
      write (503, '(8X, A2)') "1:"
      write (503, '(10X, A25)') "property: CrystalMomentum"
      write (503, '(10X, A12)') "dimension: 1" !"dimension: 2"
      write (503, '(8X, A2)') "2:"
      write (503, '(10X, A25)') "property: CrystalMomentum"
      write (503, '(10X, A12)') "dimension: 2" !"dimension: 0"
      write (503, '(6X, A9)') "nonZeros:"

      do i = 1, size(momentum_triples_list, 2)
      write (503, '(8X, A1, 1X, A1, I0, A1, I0, A1, I0, A1)') "-", "[", &
      !        &momentum_triples_list(1,i)-1, ",", &
      !        &momentum_triples_list(2,i)-1, ",", &
      !        &momentum_triples_list(3,i)-1, "]"
      momentum_triples_list(2, i) - 1, ",", &
      momentum_triples_list(1, i) - 1, ",", &
      momentum_triples_list(3, i) - 1, "]"
      end do
    end if

    close (unit=503)

  end subroutine CC4S_out_coulvertex_yaml

  ! Copy from cc4s-aims
  subroutine CC4S_out_momentum_properties_yaml(n_basis_vectors)

    !> Number of basis functions (= number of G-vectors or auxiliary field functions)
    !> used for construction of Coulomb Vertex
    integer, dimension(:), intent(in) :: n_basis_vectors

    !> Running indices q-vectors and basis functions per q-vector
    integer :: iq, ibasis
    !> Number of q-points
    integer :: n_qpoints

    n_qpoints = size(n_basis_vectors)

    open (unit=498, file='AuxiliaryField.yaml')
    write (498, '(A12)') "version: 100"
    write (498, '(A29)') "dimensionType: AuxiliaryField"
    !Leave out "properties"-section as it is not used anywhere yet
    write (498, '(A16)') "propertyIndices:"
    write (498, '(2X, A16)') "CrystalMomentum:"
    do iq = 1, n_qpoints
    do ibasis = 1, n_basis_vectors(iq)
    !write(498,'(2X, I0, A1, 1X, I0)') k-1,":", i-1
    write (498, '(2X, A1, 1X, I0)') "-", iq - 1
    end do
    end do
    close (498)

  end subroutine CC4S_out_momentum_properties_yaml




end module