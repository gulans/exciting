module mp2_energy
  use precision, only: dp
  use xlapack, only: dot_multiply

  contains


  !> gamma_ph = gamma(:, N_occ + 1:, :N_occ)
  !> gamma_hp = gamma(:, :N_occ, N_occ + 1:)

  subroutine calculate_mp2_energy(gamma_ph, gamma_hp, energies_occ, energies_virt) 
    implicit none
    !> Coulomb vertices
    complex(dp), intent(in) :: gamma_ph(:, :, :) ! N_basis x N_virt x N_occ
    complex(dp), intent(in) :: gamma_hp(:, :, :) ! N_basis x N_occ x N_virt
    !> Eigen energies of the orbitals
    real(dp), intent(in) :: energies_occ(:), energies_virt(:) ! N_occ, N_virt

    real(dp) :: E_mp2

    integer :: i, j, a, b, i_basis, N_virt, N_occ
    real(dp) :: denominator, counter
    complex(dp) int_ijab, int_abji
    integer, allocatable :: pair_a(:),pair_b(:)
    integer :: ipair,maxvirt

write(*,*)"------------ start of mp2 procedure ---------------"
    E_mp2 = 0._dp

    N_occ = size(energies_occ)
    N_virt = size(energies_virt)
if (.false.) then
    allocate(pair_a(N_virt**2))
    allocate(pair_b(N_virt**2))
    ipair=0
   ! !$omp parallel default(shared) private(ipair, a, b)
  !  ipair=0
    do maxvirt=1,N_virt
   ! !$omp do
      do a=1,maxvirt-1
        ipair=ipair+1
        pair_a(ipair)=a
        pair_b(ipair)=maxvirt
      enddo
   ! !$omp end do
   ! !$omp do
      do b=1,maxvirt-1
        ipair=ipair+1
        pair_a(ipair)=maxvirt
        pair_b(ipair)=b
      enddo
   ! !$omp end do
   ! !$omp barrier
      ipair=ipair+1
      pair_a(ipair)=maxvirt
      pair_b(ipair)=maxvirt
    enddo
!!$omp end parallel 
!!$omp parallel default(shared)
!!$omp do
    do ipair=1,N_virt**2
      a=pair_a(ipair)
      b=pair_b(ipair)
      do i=1, N_occ
        do j=1, N_occ
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
        end do
      end do
      if (a.eq.b) write(*,*) a, E_mp2
    end do

    
    write(*,*)"-------------------core-core---------------------"

   E_mp2 = 0._dp

    do ipair=1,N_virt**2
      a=pair_a(ipair)
      b=pair_b(ipair)
      do i=1, 1
        do j=1, 1
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
        end do
      end do
      if (a.eq.b) write(*,*) a, E_mp2
    end do

write(*,*)"-------------------core-valence---------------------"

    E_mp2 = 0._dp
    do ipair=1,N_virt**2
      a=pair_a(ipair)
      b=pair_b(ipair)
      do i=1, 1
        do j=2, 2
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
        end do
      end do
      if (a.eq.b) write(*,*) a, E_mp2
    end do


    write(*,*)"------------------valence-valence---------------------"
    E_mp2 = 0._dp

    do ipair=1,N_virt**2
      a=pair_a(ipair)
      b=pair_b(ipair)
      do i=2, N_occ
        do j=2, N_occ
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
        end do
      end do
      if (a.eq.b) write(*,*) a, E_mp2
    end do





    deallocate(pair_a,pair_b)
else
    do i=1, N_occ
      do j=1, N_occ
        do a=1, N_virt
          do b=1, N_virt
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
          end do
        end do
      end do
    end do
endif
write(*,*)"MP22", E_mp2
  end subroutine
 subroutine calculate_mp2_cc(gamma_ph, gamma_hp, energies_occ, energies_virt)
    implicit none
    !> Coulomb vertices
    complex(dp), intent(in) :: gamma_ph(:, :, :) ! N_basis x N_virt x N_occ
    complex(dp), intent(in) :: gamma_hp(:, :, :) ! N_basis x N_occ x N_virt
    !> Eigen energies of the orbitals
    real(dp), intent(in) :: energies_occ(:), energies_virt(:) ! N_occ, N_virt

    real(dp) :: E_mp2

    integer :: i, j, a, b, i_basis, N_virt, N_occ
    real(dp) :: denominator, counter
    complex(dp) int_ijab, int_abji
    integer, allocatable :: pair_a(:),pair_b(:)
    integer :: ipair,maxvirt

write(*,*)"------------core-core ---------------"
    E_mp2 = 0._dp

    N_occ = size(energies_occ)
    N_virt = size(energies_virt)
    do i=1, 1
      do j=1, 1
        do a=1, N_virt
          do b=1, N_virt
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
          end do
        end do
      end do
    end do
write(*,*)"MP22, c-c", E_mp2
  end subroutine

   subroutine calculate_mp2_cv(gamma_ph, gamma_hp, energies_occ, energies_virt)
    implicit none
    !> Coulomb vertices
    complex(dp), intent(in) :: gamma_ph(:, :, :) ! N_basis x N_virt x N_occ
    complex(dp), intent(in) :: gamma_hp(:, :, :) ! N_basis x N_occ x N_virt
    !> Eigen energies of the orbitals
    real(dp), intent(in) :: energies_occ(:), energies_virt(:) ! N_occ, N_virt

    real(dp) :: E_mp2

    integer :: i, j, a, b, i_basis, N_virt, N_occ
    real(dp) :: denominator, counter
    complex(dp) int_ijab, int_abji
    integer, allocatable :: pair_a(:),pair_b(:)
    integer :: ipair,maxvirt

write(*,*)"------------core-valence ---------------"
    E_mp2 = 0._dp

    N_occ = size(energies_occ)
    N_virt = size(energies_virt)
    do i=1, 1
      do j=2, 2
        do a=1, N_virt
          do b=1, N_virt
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
          end do
        end do
      end do
    end do
write(*,*)"MP22, c-v", E_mp2
  end subroutine

     subroutine calculate_mp2_vv(gamma_ph, gamma_hp, energies_occ, energies_virt)
    implicit none
    !> Coulomb vertices
    complex(dp), intent(in) :: gamma_ph(:, :, :) ! N_basis x N_virt x N_occ
    complex(dp), intent(in) :: gamma_hp(:, :, :) ! N_basis x N_occ x N_virt
    !> Eigen energies of the orbitals
    real(dp), intent(in) :: energies_occ(:), energies_virt(:) ! N_occ, N_virt

    real(dp) :: E_mp2

    integer :: i, j, a, b, i_basis, N_virt, N_occ
    real(dp) :: denominator, counter
    complex(dp) int_ijab, int_abji
    integer, allocatable :: pair_a(:),pair_b(:)
    integer :: ipair,maxvirt

write(*,*)"------------valence-valence ---------------"
    E_mp2 = 0._dp

    N_occ = size(energies_occ)
    N_virt = size(energies_virt)
    do i=2, N_occ
      do j=2, N_occ
        do a=1, N_virt
          do b=1, N_virt
            denominator = energies_occ(i) + energies_occ(j) - energies_virt(a) - energies_virt(b)
            int_ijab = dot_multiply(gamma_ph(:, a, i), gamma_hp(:, j, b))
            int_abji = dot_multiply( gamma_hp(:, i, b), gamma_ph(:, a, j))
            counter = real(int_ijab * (2 * conjg(int_ijab) - int_abji), dp)
            E_mp2 = E_mp2 + counter / denominator
          end do
        end do
      end do
    end do
write(*,*)"MP2, v-v", E_mp2
  end subroutine


end module mp2_energy

