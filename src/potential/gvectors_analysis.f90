module gvectors_analysis
  use precision, only: dp
  implicit none
contains

  recursive subroutine quicksort_real(a, left, right)
    ! In-place quicksort for a(:)
    real(dp), intent(inout) :: a(:)
    integer, intent(in)     :: left, right
    integer :: i, j
    real(dp) :: pivot, tmp

    if (left >= right) return

    i = left
    j = right
    pivot = a((left + right)/2)

    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = a(i); a(i) = a(j); a(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do

    if (left < j)  call quicksort_real(a, left, j)
    if (i    < right) call quicksort_real(a, i, right)
  end subroutine quicksort_real

  subroutine unique_counts_real(a, n, eps, n_unique, uniques, counts)
    ! Given real array a(1:n), count distinct values with absolute tolerance eps.
    ! Outputs:
    !   n_unique : number of groups
    ! Optionally:
    !   uniques(1:n) : the distinct representative values (sorted ascending)
    !   counts(1:n)  : how many entries fall into each value
    real(dp), intent(in) :: a(:)
    integer,  intent(in) :: n
    real(dp), intent(in) :: eps
    integer,  intent(out) :: n_unique
    real(dp), intent(out), optional :: uniques(:)
    integer,  intent(out), optional :: counts(:)

    real(dp), allocatable :: b(:)
    integer :: i, k

    if (n <= 0) then
      n_unique = 0
      if (present(uniques)) then
        if (size(uniques) > 0) uniques(1) = 0.0_dp
      end if
      if (present(counts)) then
        if (size(counts) > 0) counts(1) = 0
      end if
      return
    end if

    allocate(b(n))
    b = a(1:n)
    call quicksort_real(b, 1, n)

    k = 1
    if (present(uniques)) uniques(k) = b(1)
    if (present(counts))  counts(k)  = 1

    do i = 2, n
      if (abs(b(i) - b(k)) <= eps) then
        if (present(counts)) counts(k) = counts(k) + 1
      else
        k = k + 1
        b(k) = b(i)
        if (present(uniques)) uniques(k) = b(i)
        if (present(counts))  counts(k)  = 1
      end if
    end do

    n_unique = k

    ! If caller passed arrays bigger than n_unique, leaving tails unspecified is fine.

    deallocate(b)
  end subroutine unique_counts_real

  subroutine analyze_gvectors(ivgp, bvec, ngvec, eps, &
        n_unique_gperp, n_unique_gz, gperp_zero_count, &
        unique_gperp, count_gperp, unique_gz, count_gz)
    ! Compute G_perp and Gz from integer triplets ivgp(:,ig) and reciprocal basis bvec(3,3)
    ! and summarize how many distinct values each takes (within EPS).
    !
    ! Required inputs:
    !   ivgp(3,ngvec)  : integer G-index triplets
    !   bvec(3,3)      : reciprocal lattice vectors as rows b1,b2,b3 (like in your snippet)
    !   ngvec          : number of G-vectors
    !   eps            : absolute tolerance for grouping (use your epslat)
    !
    ! Required outputs:
    !   n_unique_gperp : number of distinct G_perp values
    !   n_unique_gz    : number of distinct Gz values
    !   gperp_zero_count : how many G_perp are ~0 (<= eps)
    !
    ! Optional outputs (provide arrays of size >= ngvec):
    !   unique_gperp(ngvec), count_gperp(ngvec)
    !   unique_gz(ngvec),    count_gz(ngvec)

    integer, intent(in) :: ngvec
    integer, intent(in) :: ivgp(3, ngvec)
    real(dp), intent(in) :: bvec(3,3)
    real(dp), intent(in) :: eps

    integer, intent(out) :: n_unique_gperp, n_unique_gz
    integer, intent(out) :: gperp_zero_count

    real(dp), intent(out), optional :: unique_gperp(:)
    integer,  intent(out), optional :: count_gperp(:)
    real(dp), intent(out), optional :: unique_gz(:)
    integer,  intent(out), optional :: count_gz(:)

    real(dp), allocatable :: arr_gperp(:), arr_gz(:)
    integer :: ig, ix, iy, iz
    real(dp) :: gx, gy, gz, gperp

    if (ngvec <= 0) then
      n_unique_gperp   = 0
      n_unique_gz      = 0
      gperp_zero_count = 0
      return
    end if

    allocate(arr_gperp(ngvec), arr_gz(ngvec))

    gperp_zero_count = 0

    do ig = 1, ngvec
      ix = ivgp(1, ig)
      iy = ivgp(2, ig)
      iz = ivgp(3, ig)

      gx = ix * bvec(1,1) + iy * bvec(1,2) + iz * bvec(1,3)
      gy = ix * bvec(2,1) + iy * bvec(2,2) + iz * bvec(2,3)
      gz = ix * bvec(3,1) + iy * bvec(3,2) + iz * bvec(3,3)

      gperp = sqrt(gx*gx + gy*gy)

      arr_gperp(ig) = gperp
      arr_gz(ig)    = gz

      if (gperp <= eps) gperp_zero_count = gperp_zero_count + 1
    end do

    call unique_counts_real(arr_gperp, ngvec, eps, n_unique_gperp, unique_gperp, count_gperp)
    call unique_counts_real(arr_gz,    ngvec, eps, n_unique_gz,    unique_gz,    count_gz)

    deallocate(arr_gperp, arr_gz)
  end subroutine analyze_gvectors

end module gvectors_analysis

