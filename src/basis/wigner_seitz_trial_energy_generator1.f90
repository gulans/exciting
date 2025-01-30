module wigner_seitz_trial_energy_generator1
   use modmpi, only: terminate
   use precision, only: dp

   implicit none
   
   private
   public :: generate_wigner_seitz_trial_energies1

contains

!> Generates the trial energy for a given local orbital from
!> its number of nodes by using the Wigner-Seitz rules.
   subroutine generate_wigner_seitz_trial_energies1(is,ia,l, principal_n, spr, nr, vr, wf_tolerance, energy_tolerance, en, en_m, e_p1zerro_m, e_trial_m) 

      !> angular momentum
      integer, intent(in) :: l,is,ia
      ! prinicipal quantum number
      integer, intent(in) :: principal_n
      !> number of muffin-tin radial points for the species
      integer, intent(in) :: nr
      !> species radial mesh
      real(dp), intent(in) :: spr(nr)
      !> radial component of the muffin-tin effective potential
      real(dp), intent(in) :: vr(nr)
      !> tolerance for considering the slope of a wave function as zero
      real(dp), intent(in) :: wf_tolerance
      !> target accuracy for the trial energy in the bisection 
      real(dp), intent(in) :: energy_tolerance
      real(dp), intent(out):: en, en_m, e_p1zerro_m, e_trial_m
      !> computed trial energy
      real(dp) :: e_trial
      
      ! local variables
      !> number of nodes
      integer :: nodes
      ! number of nodes
      integer :: nn
      ! major component of the radial wavefunction
      real(dp) :: p0(nr)
      ! radial derivative of p0
      real(dp) :: p1(nr)
      ! minor component of the radial wavefunction
      real(dp) :: q0(nr)
      ! radial derivative of q0
      real(dp) :: q1(nr)
      ! major component of the radial wavefunction that gets zero at muffin-tin boundary
      real(dp) :: p0_zero_at_mt(nr)
      ! energy for which the wavefunction is zero at muffin-tin boundary
      real(dp) :: energy_for_zero_wf_at_mt
      ! lower bound energy
      real(dp) :: e_lower_bound
      ! upper bound energy
      real(dp) :: e_upper_bound
      ! radial derivative of the wavefunction for lower bound energy at muffin-tin boundary
      real(dp) :: p1_lower_bound
      ! radial derivative of the wavefunction for upper bound energy at muffin-tin boundary
      real(dp) :: p1_upper_bound
      ! mean of lower and upper bound energy
      real(dp) :: e_mean
      ! radial derivative of the wavefunction for mean energy at muffin-tin boundary
      real(dp) :: p1_mean
      ! flag to pick the equation (Dirac or Schroedinger) used in rdirac
      Logical  :: dirac_eq
      ! flag to pick a quick-and-dirty algorithm for integrating the Dirac equation in rdirac
      Logical  :: sloppy
      real(8) :: e_step, e_toler, e_hi,e_lo,e_try,e_lo_all,e_hi_all
      ! Error message
      character(1024) :: message
      integer :: ie, nn_lo, nn_hi,ir

en=0d0
en_m=0d0 
e_p1zerro_m=0d0 
e_trial_m=0d0 

      if (principal_n < 1) Then
         call terminate("Error(wigner_seitz_trial_energy_generator): principal quantum number < 1")
      end if

! Schroedinger equation is used in rdirac
      dirac_eq = .false.
      sloppy = .false.

! Compute number of nodes from principal quantum number
      nodes = principal_n - l - 1

write(*,'("nodes=", I2 ," l=", I2 ," searching e where u_mt=0")')nodes,l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! u_mt=0 for nodes=nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
e_toler=energy_tolerance
e_lo_all=-15d0
e_hi_all=200d0
! ! Compute energy for which the wave function becomes 0 at the muffin-tin boundary
      e_hi=e_hi_all
      e_lo=e_lo_all
      e_try=e_lo
      Call rschroddme2(is,ia,0, l, 0, e_lo, nr, spr, vr, nn_lo, p0, p1, q0, q1)
      Call rschroddme2(is,ia,0, l, 0, e_hi, nr, spr, vr, nn_hi, p0, p1, q0, q1)
      write(*,*)e_lo,nn_lo
      write(*,*)e_hi,nn_hi
      if (.not.((nn_lo.le.nodes).and.(nn_hi.gt.nodes)))then
         write(*,*)"error needed nodes ",nodes," are not in the range"
         return
      endif
      do while (e_hi-e_lo.gt.e_toler) 
         !rschroddme2 (is,ia,m, l, k, e, nr, r, vr, nn, p0, p1, q0, q1)
         e_try=0.5d0*(e_hi + e_lo)
         Call rschroddme2(is,ia,0, l, 0, e_try, nr, spr, vr, nn, p0, p1, q0, q1)
         write(*,*)e_try,nn
         if (nn.le.nodes) then
            e_lo=e_try
         else
            e_hi=e_try
         endif
      enddo
      Call rschroddme2(is,ia,0, l, 0, e_lo, nr, spr, vr, nn_lo, p0_zero_at_mt, p1, q0, q1)
      Call rschroddme2(is,ia,0, l, 0, e_hi, nr, spr, vr, nn_hi, p0, p1, q0, q1)
      write(*,*)"e_lo:",e_lo,nn_lo
      write(*,*)"e_hi:",e_hi,nn_hi
      write(*,'("result = ", F18.12 )')e_lo
      en=e_lo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! u_mt=0 for nodes=nodes-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute energy for which the wave function with one node less becomes 0 at the muffin-tin boundary
      en_m = 0._dp
      if (nodes == 0) then
         en_m = 0d0
         e_lo = en - 5d0
      Else
         write(*,'("searching e where u_mt=0 for one node less, nodes=", I2)')nodes-1
         e_hi=e_hi_all
         e_lo=e_lo_all
         Call rschroddme2(is,ia,0, l, 0, e_lo, nr, spr, vr, nn_lo, p0, p1, q0, q1)
         Call rschroddme2(is,ia,0, l, 0, e_hi, nr, spr, vr, nn_hi, p0, p1, q0, q1)
         write(*,*)e_lo,nn_lo
         write(*,*)e_hi,nn_hi
         if (.not.((nn_lo.le.nodes-1).and.(nn_hi.gt.nodes-1)))then
            write(*,*)"error needed nodes ",nodes-1," are not in the range"
            return
         endif
         do while (e_hi-e_lo.gt.e_toler) 
            e_try=0.5d0*(e_hi + e_lo)
            Call rschroddme2(is,ia,0, l, 0, e_try, nr, spr, vr, nn, p0, p1, q0, q1)
            write(*,*)e_try,nn
            if (nn.le.nodes-1) then
               e_lo=e_try
            else
               e_hi=e_try
            endif
         enddo
         Call rschroddme2(is,ia,0, l, 0, e_lo, nr, spr, vr, nn_lo, p0, p1, q0, q1)
         Call rschroddme2(is,ia,0, l, 0, e_hi, nr, spr, vr, nn_hi, p0, p1, q0, q1)
         write(*,*)"e_lo:",e_lo,nn_lo
         write(*,*)"e_hi:",e_hi,nn_hi
         en_m=e_lo
         write(*,'("result = ", F18.12 )')e_lo
      endif
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Bisection to find du/dr=0 between e(n) and e(n-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      e_upper_bound=en
      e_lower_bound=e_lo
      Call rschroddme2(is,ia,0, l, 0, e_upper_bound, nr, spr, vr, nn, p0, p1, q0, q1)
      p1_upper_bound = p1(nr)
      Call rschroddme2(is,ia,0, l, 0, e_lower_bound, nr, spr, vr, nn, p0, p1, q0, q1)
      p1_lower_bound = p1(nr)
      if (e_upper_bound < e_lower_bound) then
         call terminate("Error(gentrialenergy): e_upper_bound < e_lower_bound")
      end if
      write(*,'("searching e where du/dr_mt=0 between e=", F18.12, " and e=", F18.12)')e_lower_bound,e_upper_bound
      do while (e_upper_bound - e_lower_bound > energy_tolerance)
         e_mean = 0.5_dp*(e_upper_bound + e_lower_bound)
         Call rschroddme2(is,ia,0, l, 0, e_mean, nr, spr, vr, nn, p0, p1, q0, q1)
         p1_mean = p1(nr)
         if (p1_mean*p1_upper_bound < 0) then
            p1_lower_bound = p1_mean
            e_lower_bound = e_mean
         else
            p1_upper_bound = p1_mean
            e_upper_bound = e_mean
         end if
      end do
      write(*,'("result du/dr_mt=0 pie e=", F18.12 )')e_lower_bound
      e_p1zerro_m=e_lower_bound

!!!!!!!!!!
      e_trial_m=0.5d0*(en + e_p1zerro_m)
   
      write(*,'("Trial energy for n=", I2," l=",I2, " e=", F18.12 )'), principal_n, l, e_trial
      write(*,*)""

   end subroutine generate_wigner_seitz_trial_energies1

end module wigner_seitz_trial_energy_generator1
