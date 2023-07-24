module trial_energy_selection

   use constants, only: y00
   use wigner_seitz_trial_energy_generator, only: generate_wigner_seitz_trial_energies
   use precision, only: dp

   implicit none

   private
   public :: select_trial_energies

contains

   !> Replace trial energies if the number of nodes of the local orbital is given in the species file.
   subroutine select_trial_energies(nlorb, lorbord, lorbl, lorbn, nspecies, idxas, nrmt, spr, veffmt, lorbe0)

      !> number of local orbitals
      integer, intent(in) :: nlorb(:)
      !> local orbital order
      integer, intent(in) :: lorbord(:, :)
      !> local-orbital angular momentum
      integer, intent(in) :: lorbl(:, :)
      !> local-orbital number of nodes
      integer, intent(in) :: lorbn(:, :, :)
      !> number of species
      integer, intent(in) :: nspecies
      !> index to atoms and species
      integer, intent(in) :: idxas(:, :)
      !> species radial mesh
      real(dp), intent(in) :: spr(:, :)
      !> number of muffin-tin radial points for each species
      integer, intent(in) :: nrmt(:)
      !> muffin-tin effective potential
      real(dp), intent(in) :: veffmt(:, :)
      !> local orbital energy
      real(dp), intent(inout) :: lorbe0(:, :, :)

      ! local variables
      Integer :: is, ilo, io
      Real(dp) :: v(maxval(nrmt))

      do is = 1, nspecies
         do ilo = 1, nlorb(is)
            do io = 1, lorbord(ilo, is) ! Looping over all local orbitals.
               ! If the number of nodes of a local orbital is given, the trial energy is computed automatically.
               if (lorbn(io, ilo, is) /= -1) then
                  v = veffmt(1:nrmt(is), idxas(1, is))*y00
                  ! Computing trial energy with Wigner-Seitz algorithm.
                  lorbe0(io, ilo, is) = generate_wigner_seitz_trial_energies(lorbl(ilo, is), lorbn(io, ilo, is), &
                                                                             spr(:, is), nrmt(is), v, 1.e-4_dp, 1.e-6_dp)
               end if
            end do
         end do
      end do

   end subroutine select_trial_energies

end module trial_energy_selection
