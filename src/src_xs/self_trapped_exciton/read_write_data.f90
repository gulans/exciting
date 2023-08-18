module read_write_data
  use precision, only: dp



  contains 

  ! phs

  subroutine read_eph_matrix(filename, eph_matrix)
    character(*), intent(in) :: filename
    complex(dp), intent(out) :: eph_matrix(:, :, :, :, :) ! e-bands x e-bands x ph-bands x k-points x q-points
  end subroutine 

  subroutine read_ph_energies(filename, ph_energies)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: ph_energies(:, :) ! ph-bands x q-points
  end subroutine 

  subroutine read_ph_eigenvectors(filename, ph_eigen_vectors)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: ph_eigen_vectors(:, :, :, :) ! cart. directions x atoms x ph-bands x q-points
  end subroutine

  subroutine read_force_constant(filename, force_constant)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: force_constant(:, :, :, :, :, :) ! cart. directions x atoms x unit cells x cart. directions x atoms x unit cells (symmetric)
  end subroutine

  subroutine read_atomic_mass(filename, at_masses)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: at_masses(:) ! atoms
  end subroutine

  ! exs

  subroutine read_ex_energies(filename, ex_energies)
    character(*), intent(in) :: filename
    real(dp), intent(out) :: ex_energies(:, :) ! excitons x q-points
  end subroutine

  subroutine read_ex_eigen_vectors(filename, ex_eigen_vectors)
    character(*), intent(in) :: filename
    complex(dp), intent(out) :: ex_eigen_vectors(:, :, :, :, :) ! exs x valence states x conduction states x k-points x q-points
  end subroutine


  
  subroutine write_phonon_amplitudes(filename, ph_amplitudes)
    character(*), intent(in) :: filename
    complex(dp), intent(in) :: ph_amplitudes(:, :) ! phonons x q-points
  end subroutine

  subroutine write_self_trapped_exciton_energies(filename, stex_energies)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: stex_energies(:) ! (excitons x q-points)
  end subroutine

  subroutine write_self_trapped_exciton_coefficients(filename, stex_coefficients)
    character(*), intent(in) :: filename
    complex(dp), intent(in) :: stex_coefficients(:, :) ! (excitons x q-points) x (excitons x q-points)
  end subroutine
end module read_write_data 





































