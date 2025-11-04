! ==============================================================================
! FILE: faddeeva_module.f90
!
! DESCRIPTION:
! This module provides a reusable Fortran interface to the C++ Faddeeva library.
! It now includes the Faddeeva w(z) function and both real and complex
! versions of the error function (erf).
! ==============================================================================
MODULE faddeeva_interface_module
    USE iso_c_binding
    USE precision, only: dp
    IMPLICIT NONE

    ! Make all functions available to other programs that USE this module
    PUBLIC :: faddeeva_from_c, erf_real, erf_cmplx

CONTAINS

    ! --- EXISTING FUNCTION ---
    FUNCTION faddeeva_from_c(z)
        COMPLEX(KIND=dp), INTENT(IN) :: z
        COMPLEX(KIND=dp)             :: faddeeva_from_c
        REAL(KIND=c_double)          :: w_real, w_imag
        INTERFACE
            SUBROUTINE faddeeva_w_c(z_real, z_imag, w_real, w_imag) BIND(C, name='faddeeva_w_c')
                USE iso_c_binding
                IMPLICIT NONE
                REAL(c_double), VALUE, INTENT(IN) :: z_real, z_imag
                REAL(c_double), INTENT(OUT)       :: w_real, w_imag
            END SUBROUTINE faddeeva_w_c
        END INTERFACE
        CALL faddeeva_w_c(REAL(z), AIMAG(z), w_real, w_imag)
        faddeeva_from_c = CMPLX(w_real, w_imag, kind=dp)
    END FUNCTION faddeeva_from_c


    ! --- NEW FUNCTIONS ADDED ---

    ! 1. Fortran wrapper for the REAL error function, erf(x)
    FUNCTION erf_real(x)
        REAL(KIND=dp), INTENT(IN) :: x
        REAL(KIND=dp)             :: erf_real
        INTERFACE
            ! The C wrapper 'erf_real_c' is a function that directly returns a double.
            FUNCTION erf_real_c(x_c) BIND(C, name='erf_real_c')
                USE iso_c_binding
                IMPLICIT NONE
                REAL(c_double), VALUE, INTENT(IN) :: x_c
                REAL(c_double)                    :: erf_real_c
            END FUNCTION erf_real_c
        END INTERFACE
        erf_real = erf_real_c(x)
    END FUNCTION erf_real

    ! 2. Fortran wrapper for the COMPLEX error function, erf(z)
    FUNCTION erf_cmplx(z)
        COMPLEX(KIND=dp), INTENT(IN) :: z
        COMPLEX(KIND=dp)             :: erf_cmplx
        REAL(KIND=c_double)          :: result_real, result_imag
        INTERFACE
            ! The interface is a subroutine, just like for w(z),
            ! returning the result via output arguments.
            SUBROUTINE erf_cmplx_c(z_real, z_imag, res_real, res_imag) BIND(C, name='erf_cmplx_c')
                USE iso_c_binding
                IMPLICIT NONE
                REAL(c_double), VALUE, INTENT(IN) :: z_real, z_imag
                REAL(c_double), INTENT(OUT)       :: res_real, res_imag
            END SUBROUTINE erf_cmplx_c
        END INTERFACE
        CALL erf_cmplx_c(REAL(z), AIMAG(z), result_real, result_imag)
        erf_cmplx = CMPLX(result_real, result_imag, kind=dp)
    END FUNCTION erf_cmplx

END MODULE faddeeva_interface_module
