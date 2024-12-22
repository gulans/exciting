module lo_recommendation1

   use constants, only: y00
   use wigner_seitz_trial_energy_generator1, only: generate_wigner_seitz_trial_energies1
   use precision, only: dp
   use modmpi, only: mpiglobal
   use precision, only: dp


   implicit none

   private
   public :: recommend_local_orbital_trial_energies1

contains

   !> Recommend trial energies for local orbitals.
   !> The trial energies are computed from the number of nodes by using the 
   !> Wigner-Seitz rules [Andersen, O.: Solid State Communications. 1973, vol. 13, no. 2, pp. 133–136.].
   !> According to these rules an appropriate trial energy can be found between the upper bound
   !> energy 
   !>   \begin{align*}
   !>    u_{l\alpha}(r_{\alpha}, E_{\text{max}})\bigg|_{r_{\alpha}=R_{{MT}}} = 0
   !>   \end{align*}
   !> and the lower bound energy 
   !>   \begin{align*}
   !>    \frac{\partial u_{l\alpha}(r_{\alpha}, E_{\text{min}})}{\partial r_{\alpha}}\bigg|_{r_{\alpha}=R_{{MT}}} = 0.
   !>   \end{align*}
   !> We then choose the linearization energy to be the mean of these two boundary values
   !>   \begin{align*}
   !>    \epsilon_{l \alpha} = \frac{E_{\text{max}} + E_{\text{min}}}{2}.
   !>   \end{align*}
   subroutine recommend_local_orbital_trial_energies1(nodesmax, lmax, nspecies, spsymb, idxas, nrmt, spr, veffmt)

      !> maximal number of nodes for which the trial energies are computed
      integer, intent(in) :: nodesmax
      !> maximal angular momentum for which the trial energies are computed
      integer, intent(in) :: lmax
      !> number of species
      integer, intent(in) :: nspecies
      !> species symbols
      character(64), intent(in) :: spsymb(:)
      !> index to atoms and species
      integer, intent(in) :: idxas(:, :)
      !> species radial mesh
      real(dp), intent(in) :: spr(:, :)
      !> number of muffin-tin radial points for each species
      integer, intent(in) :: nrmt(:)
      !> muffin-tin effective potential
      real(dp), intent(in) :: veffmt(:, :)

      ! local variables
      Integer :: is, l, nodes, fid, principal_n
      Real(dp) :: v(maxval(nrmt)), e_trial
      real(8) :: en, en_m, en_p, e_p1zerro_m, e_p1zerro_p, e_trial_m, e_trial_p
      if(mpiglobal%is_root) then

      ! Writing recommended trial energies to file
      open (newunit=fid, File='LO_RECOMMENDATION.OUT', Action='WRITE', Form='FORMATTED')
      write (fid, *) '# Recommended linearization energies computet with Wigner-Seitz rules.'
      write (fid, *) '--------------------------------------------------------------------'
      write (fid, '(" #  n_species: ", I2)') nspecies
      write (fid, '(" # n_l-channels: ", I2)') lmax + 1
      write (fid, '(" # n_nodes: ", I2)') nodesmax + 1

      do is = 1, nspecies
         ! write (fid, '(" # species: ", I3, " (", A2,")")') is, spsymb(is)
         do l = 0, lmax
            write (fid,*) 
            write (fid, '(" # species: ", A2, ", l : ", I2)') spsymb(is), l
            write (fid, *) "nodes : n : e_n(u_mt=0) : e_n-1(u_mt=0) : e_n+1(u_mt=0) : e_n-1(du/dr_mt=0) : e_n(du/dr_mt=0) : e_trial(n-1) : e_trial(n+1)"
            ! looping over numbers of nodes up to maximum number of nodes set in input
            do nodes = 0, nodesmax 
               principal_n = nodes + 1 + l
               v = veffmt(1:nrmt(is), idxas(1, is))*y00
               call generate_wigner_seitz_trial_energies1(is,1,l, principal_n, spr(:, is), nrmt(is), v, 1.e-4_dp, 1.e-6_dp, en, en_m, en_p, e_p1zerro_m, e_p1zerro_p, e_trial_m, e_trial_p)
               write (fid, '(I2," : ", I2, " : ", F18.12, " : ", F18.12," : ", F18.12," : ", F18.12," : ", F18.12," : ", F18.12," : ", F18.12)') nodes, principal_n, en, en_m, en_p, e_p1zerro_m, e_p1zerro_p, e_trial_m, e_trial_p
            end do
         end do
      end do
      close (fid)

   end if

   end subroutine recommend_local_orbital_trial_energies1

end module lo_recommendation1
