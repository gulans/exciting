#include "Faddeeva.hh"
#include <complex>

// This file provides simple C-style "wrapper" functions that can be
// easily called from Fortran. They handle the conversion between basic
// C types (double, pointers) and C++ types (std::complex).

extern "C" {

    // --- EXISTING FUNCTION ---
    // Wrapper for the Faddeeva w(z) function
    void faddeeva_w_c(double z_real, double z_imag, double* w_real_ptr, double* w_imag_ptr) {
        // Create a C++ complex number from the Fortran inputs
        std::complex<double> z_cpp(z_real, z_imag);

        // Call the main C++ Faddeeva::w function
        std::complex<double> result = Faddeeva::w(z_cpp);

        // Use the output pointers to send the result back to Fortran
        *w_real_ptr = result.real();
        *w_imag_ptr = result.imag();
    }

    // --- NEW FUNCTIONS ADDED ---

    // 1. Wrapper for the REAL error function erf(x)
    // This case is simple: it takes a double and returns a double.
    double erf_real_c(double x) {
        return Faddeeva::erf(x);
    }

    // 2. Wrapper for the COMPLEX error function erf(z)
    // This follows the same pattern as w(z): inputs are passed by value,
    // and the complex result is returned via output pointers.
    void erf_cmplx_c(double z_real, double z_imag, double* result_real_ptr, double* result_imag_ptr) {
        std::complex<double> z_cpp(z_real, z_imag);
        std::complex<double> result = Faddeeva::erf(z_cpp);
        *result_real_ptr = result.real();
        *result_imag_ptr = result.imag();
    }

} // extern "C"

