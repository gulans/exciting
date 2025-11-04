module mod_cyl_apprx
  use precision, only: dp
  use iso_c_binding, only: c_double
  use constants, only: twopi, pi
  use omp_lib
  use mod_cyl_apprx_weights, only: alpha_arr, An_arr, w_arr, tau, num_weights_real, num_weights_complex, b_max
  use faddeeva_interface_module, only: faddeeva_from_c, erf_real ! THIS IS FADDEEVA C++ implementation
  use faddeeva_test_module, only: run_faddeeva_tests ! THIS IS FOR TEST
  use Faddeeva_module, only: Faddeeva_w
  use modmain, only: intgv
  use gvectors_analysis, only: unique_counts_real

  implicit none
  private

  public :: cyl_apprx_factor, cyl_apprx_initialized, init_cyl_apprx_factor


  interface
    
  end interface

  real(dp), allocatable :: cyl_apprx_factor(:)
  logical :: cyl_apprx_initialized = .false.

contains

  subroutine init_cyl_apprx_factor(ngvec, ivgp, bvec, real_bvec, R_z, epslat)
    integer, intent(in) :: ngvec
    integer, intent(in) :: ivgp(:,:)
    real(dp), intent(in) :: bvec(:,:)        ! Reciprocal space vectors for factor calculation
    real(dp), intent(in) :: real_bvec(:,:)   ! Real vectors for calculating rho
    real(dp), intent(in) :: R_z              ! Half-height of unit cell in z direction
    real(dp), intent(in) :: epslat

    ! Local variables
    ! --------------------------------------------------------------------------------
    ! 1. Pre-calculation & General Variables
    integer :: ig, ix, iy, iz, ia, im
    integer :: nz_min, nz_max, ngz_total                             ! for unique Gz
    integer :: nx_min, nx_max, ny_min, ny_max                        ! for unique Gp
    integer :: ngp_total                                             ! for unique Gp
    real(dp) :: area, rho, cross_x, cross_y, cross_z, scale
    real(dp) :: Gx, Gy, Gz, Gp
    real(dp), allocatable :: Gp_arr(:), Gz_arr(:)                    ! for unique G-vector magnitudes 
    integer, allocatable :: ivgp_map(:,:)                            ! for storing G-vector integer multiples
    real(dp) :: sq_pi
    real(dp), parameter :: ASYMPTOTE_THRESHOLD = -700.0_dp

    ! 2. Variables Shared by Gz and Gp Calculations
    real(dp), allocatable :: alpha(:), An(:)
    complex(dp) :: term0, term1
    complex(dp), parameter :: I = (0.0_dp, 1.0_dp)

    ! 3. Variables for Gz Matrix Calculation
    real(dp), allocatable :: Gz_matrix(:,:)
    real(dp) :: alpha_val, An_val, sq_alpha, arg_exp, prefactor, cmplx_outer
    complex(dp) :: k0, k1, bracket, val

    ! 4. Variables for Gp Matrix Calculation
    real(dp), allocatable :: Gp_matrix(:,:)
    real(dp), allocatable :: A_arr(:), sq_A_arr(:)
    complex(dp), allocatable :: w(:), m_arr(:), B_arr(:)
    complex(dp) :: B_val, z0, z1, w_z0, w_z1, exp_term
    real(dp) :: A_val, sq_A

    ! 5. For final dot product
    integer :: ig_p

    ! 5. Testing and time
    real(dp) :: t0, t1
    ! Testing Gz_matrix
    complex(dp) :: term0_c, term1_c, term0_f, term1_f
    complex(dp) :: bracket_c, bracket_f
    complex(dp) :: val_c, val_f
    real(dp), allocatable :: Gz_row_c(:), Gz_row_f(:)
    real(dp) :: abs_diff
    real(dp) :: rel_diff
    real(dp), parameter :: threshold = 1.0e-10_dp
    real(dp) :: max_abs, max_rel



    if (cyl_apprx_initialized) return


    write(*,*) "Initializing cylinder Coulomb factor for ", ngvec, " G-vectors."


    ! ===================================================================================
    ! THIS ALL IS MEANT ONLY FOR THE CASE WHEN a_3 HAS JUST A VERTICAL (Z-AXIS) COMPONENT
    ! ===================================================================================

    ! ----------------------------------------------------------------------------------------------------
    ! Calculate rho by assuming it the radius of a circel of same area as the unit cell in the ab plane
    cross_x = real_bvec(2,1) * real_bvec(3,2) - real_bvec(3,1) * real_bvec(2,2)
    cross_y = real_bvec(3,1) * real_bvec(1,2) - real_bvec(1,1) * real_bvec(3,2)
    cross_z = real_bvec(1,1) * real_bvec(2,2) - real_bvec(2,1) * real_bvec(1,2)
    area = sqrt(cross_x**2 + cross_y**2 + cross_z**2)

    rho = sqrt(area / pi)

    write(*,*) "Caclculation parameters: rho = ", rho, ", Rz = ", R_z


    ! Get Gz limits
    nz_min = intgv(3, 1)
    nz_max = intgv(3, 2)

    ngz_total = nz_max - nz_min + 1

    ! allocate arrays to hold one value for each integer
    allocate( Gz_arr(nz_min:nz_max) )

    do iz = nz_min, nz_max
      ! Calculate Gz based only on the z-component of the basis
      ! (Assuming bvec(3,1) and bvec(3,2) are zero)
      Gz_arr(iz) = dble(iz) * bvec(3, 3)
      
      ! Store the corresponding integer 'iz' in the map
      !map_Gz(iz) = iz
    end do

    ! I HAVE TO RESCALE THE Gp_matrix AND Gz_matrix SO THEIR INDECES GO FROM - TO +
    
    ! Get Gx and Gy limits
    nx_min = intgv(1, 1)
    nx_max = intgv(1, 2)

    ny_min = intgv(2, 1)
    ny_max = intgv(2, 2)

    ngp_total = (nx_max - nx_min + 1) * (ny_max - ny_min + 1)


    allocate( Gp_arr(ngp_total) )

    allocate( ivgp_map(nx_min : nx_max, ny_min : ny_max) )

    ig = 0 ! Initialize the 1D counter
    do ix = nx_min, nx_max
        do iy = ny_min, ny_max
            ! Increment the 1D index
            ig = ig + 1

            ! Calculate the Cartesian components
            Gx = dble(ix) * bvec(1,1) + dble(iy) * bvec(1,2)
            Gy = dble(ix) * bvec(2,1) + dble(iy) * bvec(2,2)

            ! Store the Gp magnitude in the 1D array
            Gp_arr(ig) = sqrt(Gx*Gx + Gy*Gy)

            ! Store the (ix, iy -> ig) inverse map
            ivgp_map(ix, iy) = ig
        end do
    end do


    ! get local versions of alpha and An
    allocate(alpha(num_weights_real), An(num_weights_real))
    ! Scale alpha and A
    scale = rho*rho + R_z*R_z

    An    = An_arr / sqrt(scale)
    alpha = alpha_arr / scale

    sq_pi = sqrt(pi)

    ! ---------------------------------------------------------------------------------------------------
    ! Calculate Gz matrix

    allocate(Gz_matrix(num_weights_real, nz_min:nz_max))

    ! ------ FOR TEST ------

    allocate(Gz_row_c(num_weights_real), Gz_row_f(num_weights_real))

    max_abs = 0.0_dp
    max_rel = 0.0_dp

    write(*,'(A)') "  ig |      gz        |     max_abs    |     max_rel    "
    write(*,'(A)') "-----+----------------+----------------+----------------"

    ! ----------------------

    write(*,*) "-----------------------------------------"
    write(*,*) "Starting calculation for Gz matrix"
    call timesec(t0)


    !$OMP PARALLEL DO &
    !$OMP PRIVATE(ia, alpha_val, An_val, sq_alpha, arg_exp, ig, gz) &
    !$OMP PRIVATE(prefactor, cmplx_outer, k0, k1, term0, term1, bracket, val) &
    !$OMP DEFAULT(SHARED)

    ! Should be do ig=nz_min, nz_max
    do ig=nz_min, nz_max
      gz = Gz_arr(ig)
      ! Case when Gz=0
      ! Checked with python max abs 10^-15
      if (abs(gz) <= epslat) then
        do ia=1, num_weights_real
          sq_alpha = sqrt(alpha(ia))
          An_val = An(ia) 
          Gz_matrix(ia,ig) = An_val * sq_pi/sq_alpha * erf_real(sq_alpha * R_z)! Need to add error function
        end do
      ! If Gz != 0
      ! Found no relative error biger than 10^-8, absolute was 10^-6 and it started with k0 and k1 somehow and not cmplx_outer
      else

        ! ------ FOR TEST ------
        !write(*,'(A,F12.6,A,I4,A)') "--- Gz = ", gz, " (ig = ", ig, ") ---"
        !write(*,'(A)') "  ia |      Gz_c      |      Gz_f      |    abs_diff    |    rel_diff"
        !write(*,'(A)') "-----+----------------+----------------+----------------+----------------"
        ! ----------------------

        do ia=1, num_weights_real
          alpha_val = alpha(ia)
          An_val = An(ia)
          sq_alpha = sqrt(alpha_val)

          arg_exp = -alpha_val * R_z*R_z

          ! Case when exponent underflows
          ! Checked with python max abs diff 10^-14
          if (arg_exp <= ASYMPTOTE_THRESHOLD) then
            Gz_matrix(ia, ig) = An_val * sq_pi / sq_alpha * exp(-gz*gz / (4.0_dp * alpha_val))
          ! And when it doesnt
          else
            prefactor = An_val * sq_pi / (2.0_dp * sq_alpha) * exp(arg_exp)

            cmplx_outer = gz / (2.0_dp * sq_alpha)
            k0 = cmplx(-cmplx_outer, -sq_alpha * R_z, kind=dp) ! argument of Faddeeva is i*k, so I just multiply the i here
            k1 = cmplx(-cmplx_outer,  sq_alpha * R_z, kind=dp)

            ! ------ FOR TEST ------
            
            ! For C++
            term0_c = faddeeva_from_c(k0)
            term1_c = faddeeva_from_c(k1)

            bracket_c = exp(cmplx(0.0_dp, R_z*gz, kind=dp)) * term0_c - exp(cmplx(0.0_dp, -R_z*gz, kind=dp)) * term1_c
            val_c = prefactor * bracket_c

            Gz_row_c(ia) = real(val_c, kind=dp)

            ! For Fortran
            term0_f = Faddeeva_w(k0, 0.0_dp)
            term1_f = Faddeeva_w(k1, 0.0_dp)

            bracket_f = exp(cmplx(0.0_dp, R_z*gz, kind=dp)) * term0_f - exp(cmplx(0.0_dp, -R_z*gz, kind=dp)) * term1_f
            val_f = prefactor * bracket_f

            Gz_row_f(ia) = real(val_f, kind=dp)

            ! Now compare them
            abs_diff = abs(Gz_row_c(ia) - Gz_row_f(ia))

            if (Gz_row_c(ia) <= threshold) then
              rel_diff = abs_diff
            else 
              rel_diff = abs_diff/ Gz_row_c(ia)
            end if

            ! Check if max and save
            if (abs_diff > max_abs) then
              max_abs = abs_diff
            end if
            if (rel_diff > max_rel) then
              max_rel = rel_diff
            end if 

            ! Check next arg_exp, if yes then write out values
            arg_exp = -alpha(ia+1) * R_z*R_z
            if (arg_exp <= ASYMPTOTE_THRESHOLD) then
              write(*,'(I4," |",ES15.8," |",ES15.8," |",ES15.8)') ig, gz, max_abs, max_rel
              ! Reset max values
              max_abs = 0.0_dp
              max_rel = 0.0_dp
            end if

            ! Write out the values and differences
            !write(*,'(I4," |",ES15.8," |",ES15.8," |",ES15.8," |",ES15.8)') ia, Gz_row_c(ia), Gz_row_f(ia), abs_diff, rel_diff

            ! ----------------------

            !term0 = faddeeva_from_c(k0)
            !term1 = faddeeva_from_c(k1)
            !term0 = Faddeeva_w(k0, 0.0_dp) ! There was no relative error greater than 1e-12 
            !term1 = Faddeeva_w(k1, 0.0_dp) ! If true value was smaller than 1e-10, then I just left the absolute error 
             
            !bracket = exp(cmplx(0.0_dp, R_z*gz, kind=dp)) * term0 - exp(cmplx(0.0_dp, -R_z*gz, kind=dp)) * term1
            !val = prefactor * bracket
            !Gz_matrix(ia,ig) = real(val, kind=dp)
          end if
        end do
      end if
    end do

    !$OMP END PARALLEL DO

    ! ------ FOR TEST ------

    deallocate(Gz_row_c, Gz_row_f)
    write(*,*) "Stopping"
    stop

    ! ----------------------

    call timesec(t1)
    write(*,*) "Calculation finished, total time: ", t1 - t0, " s"
    write(*,*) "-----------------------------------------"

    deallocate(An)

    ! ---------------------------------------------------------------------------------------------------
    ! Calculate Gp matrix

    allocate(Gp_matrix(num_weights_real, ngp_total))

    ! For ease of use, I create these arrays
    allocate(A_arr(num_weights_real), sq_A_arr(num_weights_real))
    allocate(w(num_weights_complex))

    write(*,*) "-----------------------------------------"
    write(*,*) "Starting calculation for Gp matrix"
    call timesec(t0)

    A_arr = alpha * rho*rho
    sq_A_arr = sqrt(A_arr)

    do im = 1, num_weights_complex
        w(im) = w_arr(im) * exp(tau(im) * pi/b_max)
    end do

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(ig, gp, B_arr, ia, A_val, sq_A, im, B_val, z0, z1) &
    !$OMP PRIVATE(w_z0, w_z1, exp_term, term0, term1, m_arr) &
    !$OMP DEFAULT(SHARED)

    do ig=1, ngp_total
      gp = Gp_arr(ig)
      allocate(m_arr(num_weights_complex))

      ! Case when Gp=0
      ! Checked with python, max abs diff 10^-14
      if (gp <= epslat) then
        do ia=1, num_weights_real
          A_val = A_arr(ia)
          Gp_matrix(ia,ig) = (1.0_dp - exp(-A_val)) / (2.0_dp * A_val)
        end do
      ! All other cases
      ! For the ones I checked, max abs diff 10^-13
      else
        B_arr = tau * gp * rho / b_max

        do ia=1, num_weights_real
          A_val = A_arr(ia)
          sq_A = sq_A_arr(ia)

          ! Check if A_val >= 10, if yes then use closed form
          if (A_val >= 10.0_dp) then
           Gp_matrix(ia, ig) = 1 / (2 * A_val) * exp(- (gp*gp * rho*rho) / (4 * A_val))
          !other cases
          else
            do im=1, num_weights_complex
              B_val = B_arr(im)

              z0 = -B_val / (2.0_dp * sq_A)
              z1 = (2.0_dp*A_val - B_val) / (2.0_dp * sq_A)

              !w_z0 = faddeeva_from_c(I * z0) ! wofz(1j * z0)
              !w_z1 = faddeeva_from_c(I * z1 ) ! wofz(1j * z1)
              w_z0 = Faddeeva_w(I * z0, 0.0_dp) ! wofz(1j * z0)
              w_z1 = Faddeeva_w(I * z1, 0.0_dp) ! wofz(1j * z1)

              exp_term = exp(B_val - A_val)
              term0 = (1.0_dp - exp_term) / (2.0_dp * A_val)
              term1 = B_val * sq_pi / (4.0_dp * A_val * sq_A) * (w_z0 - exp_term * w_z1)

              m_arr(im) = w(im) * (term0 + term1)
            end do
            Gp_matrix(ia,ig) = real(sum(m_arr), kind=dp)
          end if
        end do
      end if
      deallocate(m_arr)
    end do

    !$OMP END PARALLEL DO

    call timesec(t1)
    write(*,*) "Calculation finished, total time: ", t1 - t0, " s"
    write(*,*) "-----------------------------------------"

    deallocate(A_arr, sq_A_arr, alpha, w)

    ! ---------------------------------------------------------------------------------------------------
    ! Map Gz and Gp to full G-vectors and take dot product for those rows

    allocate(cyl_apprx_factor(ngvec))

    write(*,*) "Mapping G-vectors and taking dot product..."
    call timesec(t0)

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(ig, ix, iy, iz, ig_p) &
    !$OMP DEFAULT(SHARED)
    do ig = 1, ngvec
        ! Get the integer multiples for this G-vector
        ix = ivgp(1, ig)
        iy = ivgp(2, ig)
        iz = ivgp(3, ig)

        ! Find the corresponding column indices in the matrices
        ig_p = ivgp_map(ix, iy)
        
        ! Take the dot product between corresponding columns
        cyl_apprx_factor(ig) = twopi * rho*rho * dot_product(Gz_matrix(:,iz), Gp_matrix(:,ig_p))
    end do
    !$OMP END PARALLEL DO

    call timesec(t1)
    write(*,*) "Time of mapping and dot product: ", t1-t0, " s"

    !write(*,*) "Everythings finished"
    !write(*,*) "Stopping"
    !stop

    cyl_apprx_initialized = .true.

    deallocate(Gz_matrix, Gp_matrix)

    ! -----------------------------------------------------------------------------------------------------
    ! For Saving data to .txt
    
    ! --- Variables for file writing ---  
    !integer :: debug_unit
    !character(len=256) :: debug_filename

    ! ====================================================================
    ! 1. OPEN THE DEBUG FILE
    ! Open a file to store the intermediate values.
    ! This is done *before* the loop for efficiency.
    ! ====================================================================
    !debug_unit = 10
    !debug_filename = "cmplx_outer.txt"
    !open(unit=debug_unit, file=debug_filename, status='replace', action='write')
    ! Write a header to make the file easier to understand
    !write(debug_unit, '(A)') "# Gz matrix column"
    
    !write(debug_unit, '(8(ES24.15E3, 1X))') real(k0), imag(k0), real(k1), imag(k1), real(term0), imag(term0), real(term1), imag(term1)

    ! ====================================================================
    ! 3. CLOSE THE DEBUG FILE
    ! Always close the file unit when you're done with it.
    ! ====================================================================
    !close(debug_unit)
    !print *, "Debug values written to ", trim(debug_filename)
    !stop

    end subroutine init_cyl_apprx_factor


end module mod_cyl_apprx