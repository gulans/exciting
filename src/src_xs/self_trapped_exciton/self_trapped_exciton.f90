module self_trapped_exciton
  use precision, only: dp 
  use read_write_data
  use coupled_equations
  use post_processing

  contains 


  subroutine main()
    complex(dp), allocatable :: eph_matrix(:, :, :, :, :), &
                                ph_eigen_vectors(:, :, :, :), &
                                ph_amplitudes(:, :), &
                                ex_eigen_vectors(:, :, :, :, :), &
                                stex_coefficients(:, :)
    real(dp), allocatable :: ph_energies(:, :), &
                             force_constant(:, :, :, :, :, :), &
                             at_masses(:), &
                             ex_energies(:, :), &
                             stex_energies(:)

    real(dp), allocatable :: stex_energies_before(:)
    logical :: converged

    ! Read in phonons 
    call read_eph_matrix(filename, eph_matrix)
    call read_ph_energies(filename, ph_energies)
    call read_ph_eigenvectors(filename, ph_eigen_vectors)
    call read_force_constant(filename, force_constant)
    
    ! Read in excitons
    call read_ex_energies(filename, ex_energies)
    call read_ex_eigen_vectors(filename, ex_eigen_vectors)

    call calculate_exph_matrix(ex_eigen_vectors, eph_matrix, exph_matrix)
    call guess_stex_coefficients(ex_energies, stex_coefficients)

    do while(.not. converged)

      call calculate_phonon_amplitude(exph_matrix, stex_coefficients, ph_energies, ph_amplitudes)
      call calculate_stex_hamiltonian(ph_amplitudes, exph_matrix, ex_energies, stex_coefficients)
      call diagonalize(stex_energies, stex_coefficients)

      converged = all_close(stex_energies, stex_energies_before)
      stex_energies_before = stex_energies_before
    end do 

    ! Write data 
    call write_phonon_amplitudes(filename, ph_amplitudes)
    call write_self_trapped_exciton_energies(filename, stex_energies)
    call write_self_trapped_exciton_coefficients(filename, stex_coefficients)

  end subroutine


  subroutine post_processing()
    ! do stuff here 
  end subroutine

end module 

