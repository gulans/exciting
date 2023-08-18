module coupled_equations
  use precision, only: dp 
  use math_utils, only: all_close
  use xlapack

  contains 

  subroutine calculate_exph_matrix(ex_eigen_vectors, eph_matrix, exph_matrix)
    complex(dp), intent(in) :: ex_eigen_vectors(:, :, :, :, :) ! exs x valence states x conduction states x k-points x q-points
    complex(dp), intent(in) :: eph_matrix(:, :, :, :, :) ! e-bands x e-bands x ph-bands x k-points x q-points
    complex(dp), intent(out) :: exph_matrix(:, :, :, :, :) ! ex-bands x ex-bands x ph-bands x q-points x q-points
  end subroutine 

  subroutine guess_stex_coefficients(ex_energies, stex_coefficients)
    real(dp), intent(in) :: ex_energies(:, :) ! excitons x q-points
    complex(dp), intent(out) :: stex_coefficients(:, :) ! excitons x q-points
  end subroutine

  subroutine calculate_phonon_amplitude(exph_matrix, stex_coefficients, ph_energies, ph_amplitudes)
    complex(dp), intent(in) :: exph_matrix(:, :, :, :, :) ! exs x valence states x conduction states x k-points x q-points
    complex(dp), intent(in) :: stex_coefficients(:, :) ! excitons x q-points
    real(dp), intent(in) :: ph_energies(:, :) ! ph-bands x q-points
    complex(dp), intent(out) :: ph_amplitudes(:, :) ! phonons x q-points
  end subroutine

  subroutine calculate_stex_hamiltonian(ph_amplitudes, exph_matrix, ex_energies, stex_hamiltonian)
    complex(dp), intent(in) :: ph_amplitudes(:, :) ! phonons x q-points
    complex(dp), intent(in) :: exph_matrix(:, :, :, :, :) ! exs x valence states x conduction states x k-points x q-points
    real(dp), intent(in) :: ex_energies(:, :) ! excitons x q-points
    complex(dp), intent(out) :: stex_hamiltonian(:, :) ! (excitons x q-points) x (excitons x q-points)
  end subroutine 

  subroutine solve_coupled_equations(ex_energies, exph_matrix, stex_energies, stex_coefficients)
    real(dp), intent(in) :: ex_energies(:, :) ! excitons x q-points
    complex(dp), intent(in) :: exph_matrix(:, :, :, :, :) ! exs x valence states x conduction states x k-points x q-points
    real(dp), intent(out) :: stex_energies(:) ! (excitons x q-points)
    complex(dp), intent(out) :: stex_coefficients(:, :) ! (excitons x q-points) x (excitons x q-points)

    real(dp), allocatable :: stex_energies_before(:)
    logical :: converged

    converged = .false.
    stex_energies_before = 0._dp

    do while(.not. converged)

      call calculate_phonon_amplitude(exph_matrix, stex_coefficients, ph_energies, ph_amplitudes)
      call calculate_stex_hamiltonian(ph_amplitudes, exph_matrix, ex_energies, stex_coefficients)
      call diagonalize(stex_energies, stex_coefficients)

      converged = all_close(stex_energies, stex_energies_before)
      stex_energies_before = stex_energies_before
    end do 

  end subroutine 
end module coupled_equations