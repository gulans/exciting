subroutine msbesselic1 (x, rez )
       Implicit None
! arguments
      complex (8), Intent (In) :: x
      complex (8), Intent (Out) :: rez

rez=cosh(x)/x




end subroutine

