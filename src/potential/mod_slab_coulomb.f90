!> This module holds the pre-calculated geometric factor for the 2D slab Coulomb potential.
module mod_slab_coulomb
  use precision, only: dp
  use constants, only: fourpi, twopi 
  implicit none
  private

  public :: slab_coulomb_factor, slab_coulomb_initialized, init_slab_coulomb_factor

  real(dp), allocatable :: slab_coulomb_factor(:)
  logical :: slab_coulomb_initialized = .false.

contains

  subroutine init_slab_coulomb_factor(ngvec, ivgp, bvec, r_c, epslat)
    integer, intent(in) :: ngvec
    integer, intent(in) :: ivgp(:,:)
    real(dp), intent(in) :: bvec(:,:)
    real(dp), intent(in) :: r_c
    real(dp), intent(in) :: epslat

    ! local variables
    integer :: ig
    real(dp) :: Gx, Gy, Gz, G_parallel
    integer :: ix, iy, iz
    real(dp) :: factor
    
    ! ---- MODIFICATION START ----
    ! File I/O variables
    !integer :: file_unit = 10
    !character(len=20) :: filename = 'gvectors.dat'
    ! ---- MODIFICATION END ----


    if (slab_coulomb_initialized) return

    allocate(slab_coulomb_factor(ngvec))

    ! ---- MODIFICATION START ----
    ! Open the file to write the G-vectors
    !open(unit=file_unit, file=filename, status='replace', action='write')
    ! Write a header for clarity
    !write(file_unit, '(A)') '# ix iy iz Gx Gy Gz'
    ! ---- MODIFICATION END ----

    write(*,*) "Initializing slab Coulomb factor for ", ngvec, " G-vectors."

    do ig = 1, ngvec
      ix = ivgp(1, ig)
      iy = ivgp(2, ig)
      iz = ivgp(3, ig)

      Gx = ix * bvec(1, 1) + iy * bvec(1, 2) + iz * bvec(1, 3)
      Gy = ix * bvec(2, 1) + iy * bvec(2, 2) + iz * bvec(2, 3)
      Gz = ix * bvec(3, 1) + iy * bvec(3, 2) + iz * bvec(3, 3)

      G_parallel = sqrt(Gx**2 + Gy**2)

      if (G_parallel > epslat) then
        factor = fourpi * (1.0d0 + exp(-G_parallel * r_c)*(Gz * sin(Gz * r_c) / G_parallel - cos(Gz * r_c))) / (G_parallel**2 + Gz**2)
      else if (G_parallel <= epslat .and. abs(Gz) > epslat) then
        factor = fourpi * (1.0d0 - cos(Gz * r_c) - Gz * r_c * sin(Gz * r_c)) / (Gz**2)
      else
        factor = -twopi * r_c**2
      endif

        ! ---- MODIFICATION START ----
      ! Write the integer and cartesian components to the file
      !write(file_unit, '(3I5, 4F18.10)') ix, iy, iz, Gx, Gy, Gz, factor
      ! ---- MODIFICATION END ----

      slab_coulomb_factor(ig) = factor
    end do
    
    ! ---- MODIFICATION START ----
    ! Close the file
    !close(file_unit)
    !write(*,*) "G-vectors written to ", trim(filename)
    ! ---- MODIFICATION END ----

    slab_coulomb_initialized = .true.

  end subroutine init_slab_coulomb_factor

end module mod_slab_coulomb
