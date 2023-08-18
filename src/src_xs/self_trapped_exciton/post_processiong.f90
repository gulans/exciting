module post_processing
  use precision, only: dp 

  contains

  subroutine calculate_atomic_displacement(at_masses, ph_energies, ph_eigen_vectors, ph_amplitudes, at_displacements)
    real(dp), intent(in) :: at_masses(:) ! atoms
    real(dp), intent(in) :: ph_energies(:, :) ! ph-bands x q-points
    real(dp), intent(in) :: ph_eigen_vectors(:, :, :, :) ! cart. directions x atoms x ph-bands x q-points
    complex(dp), intent(in) :: ph_amplitudes(:, :) ! phonons x q-points
    real(dp), intent(out) :: at_displacements(:, :, :) ! cart. directions x atoms x unit cells 
  end subroutine

  subroutine calculate_stex_real_space_distrubution(stex_coefficients, ex_real_space_distr, stex_real_space_distrubution)
    complex(dp), intent(in) :: stex_coefficients(:) ! (excitons x q-points)
    real(dp), intent(in) :: ex_real_space_distr(:, :, :) ! r-points x excitons x q-points
    real(dp), intent(out) :: stex_real_space_distrubution(:) ! r-points 
  end subroutine 
end module post_processing