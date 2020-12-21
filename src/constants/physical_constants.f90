! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! For all new units, please:
! a) Define using CODATA 2018 (see reference below)
! b) Define using double precision (dp): elec_mass = 9.10938370e-31_dp;
! c) Use descriptive naming convention

! See the 20th May 2019 redefinition of the SI units and CODATA 2018.
! Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor (2019)
! "The 2018 CODATA Recommended Values of the Fundamental Physical Constants"
! (Web Version 8.0). Database developed by J. Baker, M. Douma, and S. Kotochigova.
! Available at http://physics.nist.gov/constants,
! National Institute of Standards and Technology, Gaithersburg, MD 20899.

!> Physical constants 
module physical_constants
  use precision, only: dp 
  implicit none
  private 

  ! TODO(Alex) Issue #20. Update to CODATA 2018 physical units if required
  !> Boltzmann constant in Hartree/kelvin (CODATA 2006)
  Real(8), Public, Parameter :: kboltz = 3.166815343d-6

  !> Electron mass in kg (CODATA 2018)
  Real(dp), Public, Parameter :: elec_mass   = 9.10938370e-31_dp;

end module physical_constants