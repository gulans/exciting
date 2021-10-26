!> Fortran 95 explicit interfaces for lapack f77 routines calculating array products.
!> This ensures that the compiler is able to fully optimise.
module multiplication_interfaces
  use precision, only: dp
  
  interface 

!---------------
! BLAS Level 1 | 
!---------------
    ! dot product
    real(dp) function ddot(n, dx, incx, dy, incy)
      import :: dp 
      implicit none
      integer,  intent(in) :: n
      real(dp), intent(in) :: dx(*)
      integer,  intent(in) :: incx
      real(dp), intent(in) :: dy(*)
      integer,  intent(in) :: incy
    end function ddot

    complex(dp) function zdotc(n, zx, incx, zy, incy)
      import :: dp 
      implicit none
      integer,     intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer,     intent(in) :: incx
      complex(dp), intent(in) :: zy(*)
      integer,     intent(in) :: incy
    end function zdotc

    complex(dp) function zdotu(n, zx, incx, zy, incy)
      import :: dp 
      implicit none
      integer,     intent(in) :: n
      complex(dp), intent(in) :: zx(*)
      integer,     intent(in) :: incx
      complex(dp), intent(in) :: zy(*)
      integer,     intent(in) :: incy
    end function zdotu

  
!---------------
! BLAS Level 2 | 
!---------------
    ! general matrix vector multiplication
    subroutine dgemv(trans, m, n, alpha, A, lda, X, incx, beta, Y, incy)
      import :: dp 
      implicit none
      character, intent(in)    :: trans 
      integer,   intent(in)    :: m
      integer,   intent(in)    :: n
      real(dp),  intent(in)    :: alpha
      integer,   intent(in)    :: lda
      real(dp),  intent(in)    :: A(lda, *)
      real(dp),  intent(in)    :: X(*)
      integer,   intent(in)    :: incx
      real(dp),  intent(in)    :: beta
      real(dp),  intent(inout) :: Y(*)
      integer,   intent(in)    :: incy
    end subroutine dgemv

    subroutine zgemv(trans, m, n, alpha, A, lda, X, incx, beta, Y, incy)
      import :: dp 
      implicit none
      character,   intent(in)    :: trans
      integer,     intent(in)    :: m
      integer,     intent(in)    :: n
      complex(dp), intent(in)    :: alpha
      complex(dp), intent(in)    :: A(lda, *)
      integer,     intent(in)    :: lda
      complex(dp), intent(in)    :: X(*)
      integer,     intent(in)    :: incx
      complex(dp), intent(in)    :: beta
      complex(dp), intent(inout) :: Y(*)
      integer,     intent(in)    :: incy
    end subroutine zgemv

    ! hermitian matrix vector multiplication
    subroutine dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
      import :: dp 
      implicit none 
      character, intent(in)    :: uplo
      integer,   intent(in)    :: n
      real(dp),  intent(in)    :: alpha
      integer,   intent(in)    :: LDA
      real(dp),  intent(in)    :: A(lda, *)
      real(dp),  intent(in)    :: x(*)
      integer,   intent(in)    :: incx
      real(dp),  intent(in)    :: beta
      real(dp),  intent(inout) :: Y(*)
      integer,   intent(in)    :: incy
    end subroutine dsymv

    subroutine zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
      import :: dp 
      implicit none 
      character,   intent(in)    :: uplo
      integer,     intent(in)    :: n
      complex(dp), intent(in)    :: alpha
      integer,     intent(in)    :: LDA
      complex(dp), intent(in)    :: A(lda, *)
      complex(dp), intent(in)    :: x(*)
      integer,     intent(in)    :: incx
      complex(dp), intent(in)    :: beta
      complex(dp), intent(inout) :: Y(*)
      integer,     intent(in)    :: incy
    end subroutine zhemv

!---------------
! BLAS Level 3 | 
!---------------
    ! general matrix matrix multiplication
    subroutine dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      import :: dp 
      implicit none
      character, intent(in)    :: transa
      character, intent(in)    :: transb
      integer,   intent(in)    :: m
      integer,   intent(in)    :: n 
      integer,   intent(in)    :: k
      real(dp),  intent(in)    :: alpha
      integer,   intent(in)    :: lda
      real(dp),  intent(in)    :: A(lda, *)
      integer,   intent(in)    :: ldb
      real(dp),  intent(in)    :: B(ldb, *)
      real(dp),  intent(in)    :: beta
      integer,   intent(in)    :: ldc
      real(dp),  intent(inout) :: C(ldc, *)
    end subroutine dgemm

    subroutine zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      import :: dp 
      implicit none
      character,   intent(in)    :: transa
      character,   intent(in)    :: transb
      integer,     intent(in)    :: m
      integer,     intent(in)    :: n 
      integer,     intent(in)    :: k
      complex(dp), intent(in)    :: alpha
      integer,     intent(in)    :: lda
      complex(dp), intent(in)    :: A(lda, *)
      integer,     intent(in)    :: ldb
      complex(dp), intent(in)    :: B(ldb, *)
      complex(dp), intent(in)    :: beta
      integer,     intent(in)    :: ldc
      complex(dp), intent(inout) :: C(ldc, *)
    end subroutine zgemm

    ! hermitian matrix matrix multiplication
    subroutine dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
      import :: dp 
      implicit none
      character, intent(in)    :: side
      character, intent(in)    :: uplo
      integer,   intent(in)    :: m
      integer,   intent(in)    :: n
      real(dp),  intent(in)    :: alpha
      integer,   intent(in)    :: lda
      real(dp),  intent(in)    :: A(lda, *)
      integer,   intent(in)    :: ldb
      real(dp),  intent(in)    :: B(ldb, *)
      real(dp),  intent(in)    :: beta
      integer,   intent(in)    :: ldc
      real(dp),  intent(inout) :: C(ldc, *)
    end subroutine dsymm

    subroutine zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
      import :: dp 
      implicit none
      character,    intent(in)    :: side
      character,    intent(in)    :: uplo
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      complex(dp),  intent(in)    :: alpha
      integer,      intent(in)    :: lda
      complex(dp),  intent(in)    :: A(lda, *)
      integer,      intent(in)    :: ldb
      complex(dp),  intent(in)    :: B(ldb, *)
      complex(dp),  intent(in)    :: beta
      integer,      intent(in)    :: ldc
      complex(dp),  intent(inout) :: C(ldc, *)
    end subroutine zhemm
  end interface
  
end module multiplication_interfaces