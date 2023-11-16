!> This module provides functionalities to compute the electrostatic potential
!> from a given charge density distribution, i.e., solving Poisson's equation,
!> using Weinert's method [1].
!>
!> [1]: [M. Weinert. Solution of Poisson's equation: Beyond Ewald-type methods. 
!> *J. Math. Phys.* **22**, 2433(1981)](https://doi.org/10.1063/1.524800)
module weinert
  use precision, only: dp

  implicit none
  private

  public :: surface_ir, multipoles_ir, poisson_ir
  public :: poisson_and_multipoles_mt, match_bound_mt
  public :: poisson_and_multipoles_mt_yukawa, multipoles_ir_yukawa, pseudocharge_gspace_yukawa, poisson_ir_yukawa, pseudocharge_rspace

  contains

    !> This subroutine evaluates an interstitial function given by the Fourier series
    !> \[ f({\bf r}) = \sum_{\bf G} \hat{f}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
    !> on the surface of the muffin-tin spheres.
    !>
    !> The function evaluated on the surface of the muffin-tin sphere \(\alpha\) is given by
    !> \[ f^{{\rm SF},\alpha}(\theta,\phi) = f({\bf \tau}_\alpha + R_\alpha \hat{\bf r}(\theta,\phi)) = 
    !> \sum_{l,m} f^{{\rm SF},\alpha}_{lm} \, Y_{lm}(\theta,\phi) \;,\]
    !> with
    !> \[ f^{{\rm SF},\alpha}_{lm} = 4\pi \, {\rm i}^l \sum_{\bf G} {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}
    !> \hat{f}({\bf G+p}) \, j_l(|{\bf G+p}| R_\alpha) \, Y_{lm}^\ast({\bf \widehat{G+p}}) \;,\]
    !> where \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) which is centered at \({\bf \tau}_\alpha\)
    !> and \(j_l(x)\) are the spherical Bessel functions.
    subroutine surface_ir( lmax, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, fsf)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: rmt
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function
      complex(dp), intent(in) :: fig(:)
      !> function values on the muffin-tin sphere surfaces given by \(f^{{\rm SF},\alpha}_{lm}\)
      complex(dp), intent(out) :: fsf(:,:)

      integer :: is, ia, ias
    
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call surface_ir_single_mt( lmax, rmt(is), ngp, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
                                     fig, fsf(:,ias))
        end do
      end do
    end subroutine
    !> Same as [[surface_ir(subroutine)]] but for a single muffin-tin sphere.
    subroutine surface_ir_single_mt( lmax, rmt, ngp, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, fsf)
      use constants, only: zzero, fourpi
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function
      complex(dp), intent(in) :: fig(:)
      !> function values on the muffin-tin sphere surfaces given by \(f^{{\rm SF},\alpha}_{lm}\)
      complex(dp), intent(out) :: fsf(:)

      integer :: l, m, lm, igp, ifg, ig(3)
      complex(dp) :: z1, z2, zil

      real(dp), allocatable :: rl3(:)
    
      fsf = zzero
      
      allocate( rl3(0:lmax))
    
      rl3(0) = rmt**3
      do l = 1, lmax
        rl3(l) = rl3(l-1)*rmt
      end do
!$omp parallel default(shared) private(igp,ig,ifg,z1,z2,zil,l,m,lm) reduction(+:fsf)
!$omp do
      do igp = 1, ngp
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        z1 = fourpi*sfacgp(igp)*fig(ifg)
        zil = cmplx( 1._dp, 0._dp, dp)
        do l = 0, lmax
          z2 = z1*zil*jlgpr( l, igp)
          do m = -l, l
            lm = l*(l+1) + m + 1
            fsf(lm) = fsf(lm) + z2*conjg( ylmgp( lm, igp))
          end do
          zil = cmplx( -aimag(zil), dble(zil), dp)
        end do
      end do
!$omp end do
!$omp end parallel

      deallocate( rl3)
    end subroutine

    !> This subroutine calculates the multipole moments corresponding to the extension of
    !> an interstitial function given by the Fourier series
    !> \[ f({\bf r}) = \sum_{\bf G} \hat{f}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
    !> into the interior of the muffin-tin spheres.
    !>
    !> The multipole moments corresponding to the extension of the function inside the the muffin-tin sphere \(\alpha\)
    !> are given by
    !> \[ q^{{\rm IR},\alpha}_{lm} = 4\pi \, {\rm i}^l \, R_\alpha^{l+3} \sum_{\bf G} {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}
    !> \hat{f}({\bf G+p}) \, \frac{j_{l+1}(|{\bf G+p}| R_\alpha)}{|{\bf G+p}| R_\alpha} \, Y_{lm}^\ast({\bf \widehat{G+p}}) \;,\]
    !> where \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) which is centered at \({\bf \tau}_\alpha\)
    !> and \(j_l(x)\) are the spherical Bessel functions.
    subroutine multipoles_ir( lmax, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, qlm)
      use modinput
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: rmt
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
      complex(dp), intent(in) :: fig(:)
      !> multipole moments of the function's extension inside the muffin-tin spheres
      complex(dp), intent(out) :: qlm(:,:)

      integer :: is, ia, ias
    
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call multipoles_ir_single_mt( lmax, rmt(is), ngp, gpc, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
                                        fig, qlm(:,ias), epslat=input%structure%epslat)
        end do
      end do
    end subroutine
    !> Same as [[multipoles_ir(subroutine)]] but for a single muffin-tin sphere.
    subroutine multipoles_ir_single_mt( lmax, rmt, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, fig, qlm, epslat)
      use constants, only: zzero, fourpi, y00
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
      complex(dp), intent(in) :: fig(:)
      !> multipole moments of the function's extension inside the muffin-tin spheres
      complex(dp), intent(out) :: qlm(:)
      !> threshold below which \(|{\bf G+p}|\) is considered zero
      real(dp), optional, intent(in) :: epslat

      integer :: i, l, m, lm, igp, ngpf, ngpz, ifg, ig(3)
      real(dp) :: eps
      complex(dp) :: z1, z2, zil, figzero

      integer, allocatable :: igp_finite(:), igp_zero(:)
      real(dp), allocatable :: rl3(:)
    
      eps = 1.d-12
      if( present( epslat)) eps = epslat

      qlm = zzero
      figzero = zzero
      
      allocate( rl3(0:lmax))
    
      rl3(0) = rmt**3
      do l = 1, lmax
        rl3(l) = rl3( l-1)*rmt
      end do

      igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > eps, i=1, ngp)])
      ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg,z1,z2,zil,l,m,lm) reduction(+:qlm)
!$omp do
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        z1 = fourpi*sfacgp(igp)*fig(ifg)/(gpc(igp)*rmt)
        zil = cmplx( 1._dp, 0._dp, dp)
        do l = 0, lmax
          z2 = z1*zil*rl3(l)*jlgpr( l+1, igp)
          do m = -l, l
            lm = l*(l+1) + m + 1
            qlm(lm) = qlm(lm) + z2*conjg( ylmgp( lm, igp))
          end do
          zil = cmplx( -aimag(zil), dble(zil), dp)
        end do
      end do
!$omp end do
!$omp end parallel
      igp_zero = pack( [(i, i=1, ngp)], [(gpc(i) <= eps, i=1, ngp)])
      ngpz = size( igp_zero)
      do i = 1, ngpz
        igp = igp_zero(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        qlm(1) = qlm(1) + fourpi/3._dp*rl3(0)*y00*fig(ifg)
      end do

      deallocate( rl3)
    end subroutine

    !> This subroutine solves Poisson's equation for a given complex charge density
    !> contained in an isolated muffin-tin sphere using the Green's function approach
    !> and calculates the multipole moments of the charge distribution.
    !>
    !> The general electrostatic potential arising from an arbitrary charge distribution
    !> \[ n({\bf r}) = n({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} n^\alpha_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> inside the muffin-tin sphere \(\alpha\) is given by
    !> \[ v_{\rm sph}({\bf r}) = v_{\rm sph}({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} v_{\rm sph}[n^\alpha_{lm}](r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> with
    !> \[ v_{\rm sph}[n^\alpha_{lm}](r) = \frac{4\pi}{2l+1} \left[ 
    !> \frac{1}{r^{l+1}} \int_0^r s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s + 
    !> r^l \int_r^{R_\alpha} \frac{n^\alpha_{lm}(s)}{s^{l-1}} \, {\rm d}s - 
    !> \frac{r^l}{R_\alpha^{2l+1}} \int_0^{R_\alpha} s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s \right] \;.\]
    !>
    !> In addition, the multipole moments of the charge distribution in the sphere are calculated as
    !> \[ q^{{\rm MT},\alpha}_{lm} = \int_0^{R_\alpha} s^{l+2} \, n^\alpha_{lm}(s) \, {\rm d}s \,. \]
    !>
    !> @note \(v_{\rm sph}({\rm r})\) as defined above goes to zero on the muffin-tin surface.@endnote
    subroutine poisson_and_multipoles_mt( lmax, nr, r, zrhomt, zvclmt, qlm)
      use constants, only: fourpi
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> number of radial grid points
      integer, intent(in) :: nr
      !> radial grid
      real(dp), intent(in) :: r(:)
      !> complex charge distribution \(n^\alpha_{lm}(r)\)
      complex(dp), intent(in) :: zrhomt(:,:)
      !> complex electrostatic potential \(v_{\rm sph}[n^\alpha_{lm}](r)\)
      complex(dp), intent(out) :: zvclmt(:,:)
      !> multipole moments of the charge distribution \(q^{{\rm MT},\alpha}_{lm}\)
      complex(dp), intent(out) :: qlm(:)

      integer :: l, m, lm, ir
      real(dp) :: t1, t2, t3, t4, t5

      real(dp), allocatable :: ri(:), rl(:), ril1(:), cf(:,:)
      real(dp), allocatable :: fr(:,:), gr(:,:)

      allocate( ri(nr), rl(nr), ril1(nr), cf(3,nr))
      allocate( fr(nr,4), gr(nr,4))
    
      ! initialise r^l and r^(-l-1)
      do ir = 1, nr
        rl(ir) = 1._dp
        ri(ir) = 1._dp/r(ir)
        ril1(ir) = ri(ir)
      end do
      lm = 0
      do l = 0, lmax
        t1 = fourpi/dble( 2*l+1)
        do m = -l, l
          lm = lm + 1
          do ir = 1, nr
            t2 = rl(ir)*r(ir)*r(ir)    ! r^(2+l)
            t3 = ril1(ir)*r(ir)*r(ir)  ! r^(1-l)
            t4 = dble( zrhomt(lm, ir))
            t5 = aimag( zrhomt(lm, ir))
            fr(ir,1) = t2*t4
            fr(ir,2) = t2*t5
            fr(ir,3) = t3*t4
            fr(ir,4) = t3*t5
          end do
          call fderiv( -1, nr, r, fr(:,1), gr(:,1), cf)
          call fderiv( -1, nr, r, fr(:,2), gr(:,2), cf)
          call fderiv( -1, nr, r, fr(:,3), gr(:,3), cf)
          call fderiv( -1, nr, r, fr(:,4), gr(:,4), cf)
          qlm(lm) = cmplx( gr(nr,1), gr(nr,2), 8)
          do ir = 1, nr
            t2 = ril1(ir)*gr(ir,1) + rl(ir)*(gr(nr,3) - gr(ir,3) - gr(nr,1)*ril1(nr)/rl(nr))
            t3 = ril1(ir)*gr(ir,2) + rl(ir)*(gr(nr,4) - gr(ir,4) - gr(nr,2)*ril1(nr)/rl(nr))
            zvclmt( lm, ir) = t1*cmplx( t2, t3, 8)
          end do
        end do
        ! update r^l and r^(-l-1)
        if( l < lmax) then
          do ir = 1, nr
            rl(ir) = rl(ir)*r(ir)
            ril1(ir) = ril1(ir)*ri(ir)
          end do
        end if
      end do
    
      deallocate( ri, rl, ril1, cf, fr, gr)
      return
    end subroutine

    !> This subroutine solves Poisson's equation for a complex charge density given in the
    !> interstitial region by 
    !> \[ n({\bf r}) = \sum_{\bf G} \hat{n}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \]
    !> and having multipole moments \(q^\alpha_{lm}\) inside the muffin-tin spheres.
    !>
    !> This is done by constructing a pseudodensity with the given multipole moments and having a 
    !> rapidly converging Fourier series. The Fourier components of the pseudodensity are given by
    !> \[ \hat{n}^{\rm ps}({\bf G+p}) = \hat{n}({\bf G+p}) + \sum_\alpha \hat{n}^{{\rm ps},\alpha}({\bf G+p}) \;.\]
    !> See [[pseudodensity_ir_single_mt(subroutine)]] for further details.
    !>
    !> From the pseudodensity the electrostatic potential in the interstitial region is obtained
    !> by solving Poisson's equation in reciprocal space, i.e., 
    !> \[ V({\bf r}) = \sum_{\bf G} \hat{V}({\bf G+p}) \, {\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf r}} \;, \]
    !> with
    !> \[ \hat{V}({\bf G+p}) = 4\pi \frac{\hat{n}^{\rm ps}({\bf G+p})}{|{\bf G+p}|^2} \;.\]
    subroutine poisson_ir( lmax, npsden, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, zrhoig, qlm, zvclig, cutoff_in,hybrid_in)
      use modinput
      use mod_lattice, only: omega
      use mod_atoms, only: nspecies, natoms, idxas
      Use mod_kpoint, only: nkptnr
      use mod_muffin_tin, only: rmt
      use constants, only: zzero, fourpi, twopi
      use mod_Gvector, only: ngrid
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> pseudodensity expansion order
      integer, intent(in) :: npsden
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> on entry: Fourier components \(\hat{n}({\bf G+p})\) of the interstitial charge density on the FFT grid
      !> on exit: Fourier components \(\hat{n}^{\rm ps}({\bf G+p})\) of the pseudodensity on the FFT grid
      complex(dp), intent(inout) :: zrhoig(:)
      !> muffin-tin multipole moments \(q^\alpha_{lm}\) of the charge density
      complex(dp), intent(in) :: qlm(:,:)
      !> option for using coulomb cutoff for solving Poisson's equation
      logical, optional, intent(in) :: cutoff_in

      logical, optional, intent(in) :: hybrid_in 
      !> Fourier components \(\hat{V}({\bf G+p})\) of the interstitial electrostatic potential on the FFT grid
      complex(dp), intent(out) :: zvclig(:)

      integer :: i, is, ia, ias, igp, ngpf, ifg, ig(3)
      real :: r_c
      integer, allocatable :: igp_finite(:)
      
      logical :: cutoff, hybrid

  if (present(hybrid_in)) then 
        hybrid=hybrid_in
  else
        hybrid=.false.
  endif

  if (present(cutoff_in)) then 
    cutoff=cutoff_in
else
    cutoff=.false.
endif

    !if (hybrid) then
    if(.false.) then      
          call pseudocharge_rspace(lmax,npsden,qlm,zrhoig)
          
    else

      ! add Fourier components of pseudodensity from multipole moments
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          call pseudodensity_ir_single_mt( lmax, rmt(is), omega, npsden, ngp, gpc, ivgp, jlgpr(:,:,is), ylmgp, sfacgp(:,ias), intgv, ivgig, igfft, &
                                           zrhoig, qlm(:,ias), epslat=input%structure%epslat)
        end do
      end do
  
    endif

      ! solve Poisson's equation in reciprocal space
!      zvclig = zzero
      igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > input%structure%epslat, i=1, ngp)])
      ngpf = size( igp_finite)



      if (cutoff) then
        r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0
        zvclig = zrhoig*twopi*r_c**2
!$omp parallel default(shared) private(i,igp,ig,ifg)
!$omp do
      ! cutoff correction for > epslat
        do i = 1, ngpf
          igp = igp_finite(i)
          ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
          ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
          zvclig(ifg) = fourpi*zrhoig(ifg)*(1d0-cos(gpc(igp)*r_c))/(gpc(igp)**2)
        end do
!$omp end do
!$omp end parallel
      else
        zvclig = zzero
      ! without cutoff correction
!$omp parallel default(shared) private(i,igp,ig,ifg)
!$omp do
        do i = 1, ngpf
          igp = igp_finite(i)
          ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
          ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
          zvclig(ifg) = fourpi*zrhoig(ifg)/(gpc(igp)**2)
        end do
!$omp end do
!$omp end parallel
      end if
    end subroutine

    !> This subroutine computes the Fourier components of a quickly converging pseudodensity with multipole moments \(q_{lm}\)
    !> in muffin-tin \(\alpha\) and adds them to the input argument `zrhoig`.
    !>
    !> The Fourier components to be added are given by
    !> \[ \hat{n}^{{\rm ps},\alpha}({\bf G+p}) = \frac{4\pi}{\Omega} {\rm e}^{-{\rm i}({\bf G+p}) \cdot {\bf \tau}_\alpha}
    !> \sum_{l,m} \left( -\frac{\rm i}{R_\alpha} \right)^l \frac{(2l + 2N + 3)!!}{(2l + 1)!!} 
    !> \frac{j_{l+N+1}(|{\bf G+p}|R_\alpha)}{(|{\bf G+p}|R_\alpha)^{N+1}} \, q^\alpha_{lm} \, Y_{lm}(\widehat{\bf G+p}) \;,\]
    !> where \(\Omega\) is the unit cell volume, \(R_\alpha\) is the radius of the muffin-tin sphere \(\alpha\) 
    !> which is centered at \({\bf \tau}_\alpha\) and \(j_l(x)\) are the spherical Bessel functions.
    !> \(N\) is the expansion order of the pseudodensity.
    subroutine pseudodensity_ir_single_mt( lmax, rmt, omega, npsden, ngp, gpc, ivgp, jlgpr, ylmgp, sfacgp, intgv, ivgig, igfft, zrhoig, qlm, epslat)
      use constants, only: zzero, fourpi, zi, y00
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> unit cell volume \(\Omega\)
      real(dp), intent(in) :: omega
      !> pseudodensity expansion order
      integer, intent(in) :: npsden
      !> total number of \({\bf G+p}\) vectors
      integer, intent(in) :: ngp
      !> lengths of \({\bf G+p}\) vectors
      real(dp), intent(in) :: gpc(:)
      !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
      integer, intent(in) :: ivgp(:,:)
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:)
      !> bounds for integer components of \({\bf G}\)
      integer, intent(in) :: intgv(3,2)
      !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
      integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
      !> map from \({\bf G}\) vector index to point in FFT grid
      integer, intent(in) :: igfft(:)
      !> Fourier components \(\hat{n}^{\rm ps}({\bf G+p})\) of the pseudodensity on the FFT grid
      complex(dp), intent(inout) :: zrhoig(:)
      !> muffin-tin multipole moments \(q^\alpha_{lm}\) of the charge density
      complex(dp), intent(in) :: qlm(:)
      !> threshold below which \(|{\bf G+p}|\) is considered zero
      real(dp), optional, intent(in) :: epslat

      integer :: i, l, m, lm, igp, ngpf, ngpz, ifg, ig(3)
      real(dp) :: eps, fpo, t1, t2
      complex(dp) :: z1, z2
    
      integer, allocatable :: igp_finite(:), igp_zero(:)
      complex(dp), allocatable :: zrp(:)

      real(dp), external :: factnm
      
      eps = 1.d-12
      if( present( epslat)) eps = epslat

      fpo = fourpi/omega
      allocate( zrp( (lmax+1)**2))
    
      ! add Fourier components of pseudodensity from multipole moments
      t1 = 1.d0
      do l = 0, lmax
        t2 = t1*factnm( 2*l + 2*npsden + 3, 2)/factnm( 2*l+1, 2)
        t1 = t1/rmt
        do m = -l, l
          lm = l*(l+1) + m + 1
          zrp(lm) = conjg( zi**l)*t2*qlm(lm)
        end do
      end do
    
      igp_finite = pack( [(i, i=1, ngp)], [(gpc(i) > eps, i=1, ngp)])
      ngpf = size( igp_finite)
!$omp parallel default(shared) private(i,igp,ig,ifg,l,m,lm,z1,z2,t1,t2)
!$omp do
      do i = 1, ngpf
        igp = igp_finite(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        z1 = fpo*conjg( sfacgp(igp))
        t1 = (gpc(igp)*rmt)**(npsden + 1)
        do l = 0, lmax
          t2 = jlgpr( l+npsden+1, igp)/t1
          z2 = zzero
          do m = -l, l
            lm = l*(l+1) + m + 1
            z2 = z2 + zrp(lm)*ylmgp( lm, igp)
          end do
          zrhoig(ifg) = zrhoig(ifg) + z1*t2*z2
        end do
      end do
!$omp end do
!$omp end parallel
      igp_zero = pack( [(i, i=1, ngp)], [(gpc(i) <= eps, i=1, ngp)])
      ngpz = size( igp_zero)
      do i = 1, ngpz
        igp = igp_zero(i)
        ig = modulo( ivgp(:,igp)-intgv(:,1), intgv(:,2)-intgv(:,1)+1) + intgv(:,1)
        ifg = igfft( ivgig( ig(1), ig(2), ig(3)))
        zrhoig(ifg) = zrhoig(ifg) + fpo/factnm( 2*npsden + 3, 2)*zrp(1)*y00
      end do
      deallocate( zrp)
    end subroutine

    !> Given a complex muffin-tin function
    !> \[ f({\bf r}) = f({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} f^\alpha_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> and an interstitial function given on the muffin-tin sphere surface by
    !> \(f^{{\rm SF},\alpha}_{lm}\) according to subroutine [[surface_ir(subroutine)]],
    !> this subroutine matches the muffin-tin function to the interstial function by
    !> adding the homogeneous function
    !> \[ f^{\rm hom}({\bf r}) = f^{\rm hom}({\bf \tau}_\alpha + r_\alpha \, \hat{\bf r}_\alpha) = 
    !> \sum_{l,m} \left[ f^{{\rm SF},\alpha}_{lm} - f^\alpha_{lm}(R_\alpha) \right] \,
    !> \left( \frac{r_\alpha}{R_\alpha} \right)^l \, Y_{lm}(\hat{\bf r}_\alpha) \]
    !> to the muffin-tin function.
    subroutine match_bound_mt( lmax, nr, r, rmt, firsf, fmt, yukawa, zbessi)
      !> maximum angular momentum \(l\)
      integer, intent(in) :: lmax
      !> number of radial grid points
      integer, intent(in) :: nr
      !> radial grid
      real(dp), intent(in) :: r(:)
      !> muffin-tin radius \(R_\alpha\)
      real(dp), intent(in) :: rmt
      !> interstitial function on muffin-tin sphere surface \(f^{{\rm SF},\alpha}_{lm}\)
      complex(dp), intent(in) :: firsf(:)
      !> muffin-tin function \(f^\alpha_{lm}(r)\)
      complex(dp), intent(inout) :: fmt(:,:)

      logical, optional, intent(in) :: yukawa
      complex(dp),optional, intent(in) :: zbessi(:,0:) !nrmtmax, 0:input%groundstate%lmaxvr+input%groundstate%npsden+1
      integer :: l, m, lm, ir
      complex(dp) :: df,zt1

      real(dp), allocatable :: rr(:,:)

      

      if ((present(yukawa)).and.yukawa) then

        lm = 0
        do l = 0, lmax
          do m = -l, l
            lm = lm + 1
            zt1 = firsf(lm) - fmt (lm, nr)
            Do ir = 1, nr
               fmt (lm, ir) = fmt (lm, ir) + zt1 &
               & * zbessi(ir,l)/zbessi(nr,l)
  

            End Do !ir
          enddo
        enddo
      
      
      
      else


        allocate( rr(nr,2), source=1._dp)
        do ir = 1, nr
          rr(ir,2) = r(ir)/rmt
        end do
     
        lm = 0
        do l = 0, lmax
          do m = -l, l
            lm = lm + 1
            df = firsf(lm) - fmt(lm,nr)
            do ir = 1, nr
              fmt(lm,ir) = fmt(lm,ir) + df*rr(ir,1)
            end do
          end do
          do ir = 1, nr
            rr(ir,1) = rr(ir,1)*rr(ir,2)
          end do
        end do

        deallocate( rr)
    endif
    end subroutine









subroutine poisson_and_multipoles_mt_yukawa( lmax, nr, r, zrhomt, zvclmt, qlm, zlambda, il, kl, is)
use modinteg
use constants, only: fourpi
!> maximum angular momentum \(l\)
integer, intent(in) :: lmax
!> number of radial grid points
integer, intent(in) :: nr
!> radial grid
real(dp), intent(in) :: r(:)
!> complex charge distribution \(n^\alpha_{lm}(r)\)
complex(dp), intent(in) :: zrhomt(:,:)
!> complex electrostatic potential \(v_{\rm sph}[n^\alpha_{lm}](r)\)
complex(dp), intent(out) :: zvclmt(:,:)
!> multipole moments of the charge distribution \(q^{{\rm MT},\alpha}_{lm}\)
complex(dp), intent(out) :: qlm(:)
Complex (8),Intent (In) :: il(nr,0:lmax), kl(nr,0:lmax),zlambda

integer, intent(in) :: is
integer :: l, m, lm, ir


complex(dp) :: f1(nr),f2(nr),g1(nr),g2(nr),zt1,zt2
!external functions
Real (8) :: factnm
External factnm



lm = 0
Do l = 0, lmax
  Do m= -l, l
      lm = lm + 1

      f1 = il(:,l) * r**2 * zrhomt(lm, :)

      call integ_cf (nr, is, f1, g1, mt_integw)
      f1 = kl(:,l) * g1

      f2 = kl(:,l) * r**2 * zrhomt(lm, :)
      call integ_cf (nr, is, f2, g2, mt_integw)
      f2= il(:,l) * (g2(nr)-g2)
      zvclmt (lm, :)=fourpi * zlambda * (f1+f2)
   Enddo
enddo

lm=0
Do l = 0, lmax
  zt1 = factnm (2*l+1, 2) / (zlambda ** l)
  zt2 = 1d0 / ( kl(nr,l)* fourpi * zlambda)
  Do m = - l, l
    lm = lm + 1
    qlm (lm) = zt1 * zt2 * zvclmt (lm,nr)
  End Do
enddo


!write(*,*)"mans"
!do ir=1, nr
  !write(*,*)dble(il(ir, 0)),",",imag(il(ir, 0))
  !write(*,*)dble(zrhomt(1, ir)),",",imag(zrhomt(1, ir))
  !write(*,*)dble(zvclmt (1, ir)),",", imag(zvclmt (1, ir))
!enddo
!stop
!write(*,*)qlm(1:3)


end subroutine


subroutine multipoles_ir_yukawa( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qi, zlambda,zilmt)
  use modinput
  use mod_atoms, only: nspecies, natoms, idxas
  use mod_muffin_tin, only: rmt,nrmtmax,nrmt
  use constants, only: fourpi,y00,zil,zzero
  !> maximum angular momentum \(l\)
  integer, intent(in) :: lmax
  !> total number of \({\bf G+p}\) vectors
  integer, intent(in) :: ngvec
  !> lengths of \({\bf G+p}\) vectors
  real(dp), intent(in) :: gpc(:)
  !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
  !integer, intent(in) :: ivgp(:,:)
  !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
  real(dp), intent(in) :: jlgpr(0:,:,:)
  !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
  complex(dp), intent(in) :: ylmgp(:,:)
  !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
  complex(dp), intent(in) :: sfacgp(:,:)
  !> bounds for integer components of \({\bf G}\)
  !integer, intent(in) :: intgv(3,2)
  !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
  !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
  !> map from \({\bf G}\) vector index to point in FFT grid
  integer, intent(in) :: igfft(:)
  !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
  complex(dp), intent(in) :: zvclir(:)
  !> multipole moments of the function's extension inside the muffin-tin spheres
  complex(dp), intent(out) :: qi(:,:)
  complex(dp), intent(in) :: zlambda
  complex(dp), intent(in) :: zilmt(0:,:) 

  integer :: is, ia, ias, ig, ifg, l, lm, m
  complex(dp) :: zt1,zt2
  real(dp) :: t1
!external functions
  Real (8) :: factnm
  External factnm

  qi (:, :) = zzero


  do is = 1, nspecies
     do ia = 1, natoms(is)
      ias = idxas( ia, is)

        Do ig = 1, ngvec
          ifg = igfft (ig)
          
          If (gpc(ig) .Gt. input%structure%epslat) Then
             zt1 = zvclir (ifg) * sfacgp (ig, ias)/((gpc(ig) ** 2)+(zlambda ** 2))
             
             lm = 0
             Do l = 0, input%groundstate%lmaxvr
                if (l .eq. 0) then 
                zt2 = zt1 *  fourpi * (rmt(is) ** 2) &
               &* ((zlambda * jlgpr (0, ig, is) * (zilmt(1,is)+ (zilmt(0,is)/(zlambda*rmt(is)))))&
               &- (gpc(ig) * ((jlgpr(0, ig, is)/(gpc(ig)*rmt(is)))-jlgpr(1, ig, is)) * zilmt(0,is)))
                else
                zt2 = zt1 *  fourpi * zil(l) * (rmt(is) ** 2) * factnm (2*l+1, 2) /(zlambda ** (l))&
               &* ((zlambda * jlgpr (l, ig, is) * zilmt(l-1,is))&
               &- (gpc(ig) * jlgpr (l-1, ig, is) * zilmt(l,is)))
                endif

                Do m = - l, l
                   lm = lm + 1
                   qi (lm, ias) = qi (lm, ias) + zt2 * conjg (ylmgp(lm, ig))
                End Do
             End Do
          Else
             t1 = fourpi * y00 * (rmt (is) ** 2) * zilmt(1,is) / zlambda
             qi (1,ias) = qi (1,ias) + t1 * zvclir (ifg)
          End If
        enddo                        
    end do
  end do

end subroutine


subroutine pseudocharge_gspace_yukawa( lmax, ngvec, gpc, jlgpr, ylmgp, sfacgp, igfft, zvclir, qlm, zlambda, zilmt)
  use mod_lattice, only: omega
  use modinput
  use mod_atoms, only: nspecies, natoms, idxas
  use mod_muffin_tin, only: rmt,nrmtmax,nrmt
  use constants, only: fourpi,y00,zil
  !> maximum angular momentum \(l\)
  integer, intent(in) :: lmax
  !> total number of \({\bf G+p}\) vectors
  integer, intent(in) :: ngvec
  !> lengths of \({\bf G+p}\) vectors
  real(dp), intent(in) :: gpc(:)
  !> integer components \({\bf G}\) of \({\bf G+p}\) vectors
  !integer, intent(in) :: ivgp(:,:)
  !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
  real(dp), intent(in) :: jlgpr(0:,:,:)
  !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
  complex(dp), intent(in) :: ylmgp(:,:)
  !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
  complex(dp), intent(in) :: sfacgp(:,:)
  !> bounds for integer components of \({\bf G}\)
  !integer, intent(in) :: intgv(3,2)
  !> map from integer components of \({\bf G}\) vectors to \({\bf G}\) vector index
  !integer, intent(in) :: ivgig(intgv(1,1):,intgv(2,1):,intgv(3,1):)
  !> map from \({\bf G}\) vector index to point in FFT grid
  integer, intent(in) :: igfft(:)
  !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
  complex(dp), intent(inout) :: zvclir(:)
  !> multipole moments of the function's extension inside the muffin-tin spheres
  complex(dp), intent(in) :: qlm(:,:)
  complex(dp), intent(in) :: zlambda

  complex(dp), intent(in) :: zilmt(0:,:)

  complex(dp) :: zrp((lmax+1)**2)
  
  
  integer :: is, ia, ias, ig, ifg, l, lm, m,npsd
  complex(dp) :: zt1,zsum1, zsum2
  real(dp) :: t1,z1,fpo
!external functions
  Real (8) :: factnm
  External factnm


  fpo = fourpi/omega
npsd=input%groundstate%npsden

  do is = 1, nspecies
      do ia = 1, natoms(is)
      ias = idxas( ia, is)
      lm=0
      Do l = 0, lmax
        t1 = 1d0 / (factnm (2*l+1, 2))
        Do m = - l, l
           lm = lm + 1
           zrp (lm) = qlm(lm, ias) * t1
           !write(*,*)"qlm",qlm(lm, ias),t1
        End Do
     End Do
                 
 

  Do ig = 1, ngvec
    ifg = igfft (ig)
    If (gpc(ig) .Gt. input%structure%epslat) Then
       z1 = gpc (ig)
       zt1 = fpo * conjg (sfacgp(ig, ias))/(z1 ** (npsd+1))
       lm = 0
       Do l = 0, lmax
          zsum1 = (jlgpr (npsd+l+1, ig, is)*(zlambda ** (l+npsd+1))/ zilmt(npsd+l+1,is))* conjg(zil(l))
          lm = lm + 1
          zsum2 = zrp (lm) * ylmgp (lm, ig)
          Do m = - l + 1, l
             lm = lm + 1
             zsum2 = zsum2 + zrp (lm) * ylmgp (lm, ig)
          End Do
          zvclir (ifg) = zvclir (ifg) +  zt1 * zsum1 * zsum2 
       End Do

    Else
       !Sie divi labojumi no oriģinālā szpotcoul:
       zt1 = (fpo * y00 * (rmt(is) ** (npsd+1))&
      &* (zlambda ** (npsd+1))) / ((factnm (2*npsd+3, 2))* zilmt(npsd+1,is))
      
       zvclir (ifg) = zvclir (ifg) + zt1 * qlm(1, ias) 
     
    End If
 End Do
End Do
End Do

!do ig=1,20
!  ifg = igfft (ig)
!  write(*,*)zvclir (ifg)
!enddo
end subroutine

subroutine poisson_ir_yukawa( lmax, ngvec, gpc, igfft, zrhoir, zlambda,zvclir,cutoff)
  use mod_lattice, only: omega
  Use mod_kpoint, only: nkptnr
  use modinput
  use constants, only: zzero
 ! use mod_atoms, only: nspecies, natoms, idxas
 ! use mod_muffin_tin, only: rmt,nrmtmax,nrmt
  use constants, only: fourpi 
  !> maximum angular momentum \(l\)
  integer, intent(in) :: lmax
  !> total number of \({\bf G+p}\) vectors
  integer, intent(in) :: ngvec
  !> lengths of \({\bf G+p}\) vectors
  real(dp), intent(in) :: gpc(:)
  !> map from \({\bf G}\) vector index to point in FFT grid
  integer, intent(in) :: igfft(:)
  !> Fourier components \(\hat{f}({\bf G+p})\) of the function on the FFT grid
  complex(dp), intent(in) :: zrhoir(:)

  complex(dp), intent(in) :: zlambda
  complex(dp), intent(out) :: zvclir(:)
  logical, intent(in) :: cutoff
  
  
  integer :: ig, ifg
  real(dp) :: r_c
!external functions
  Real (8) :: factnm
  External factnm

  zvclir=zzero
r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0
Do ig = 1, ngvec
  ifg = igfft (ig)

if (cutoff) then

  If (gpc(ig) .Gt. input%structure%epslat) Then
     zvclir (ifg) = fourpi * zrhoir (ifg) *&
        (1d0 - exp(-zlambda*r_c)* ( zlambda * sin(gpc(ig)*r_c) / gpc(ig) + cos(gpc(ig)*r_c) ) )/&
        ((gpc(ig) ** 2) + (zlambda ** 2))
  Else
     zvclir (ifg) = fourpi * zrhoir(ifg)* (1d0 - exp(-zlambda*r_c)*(zlambda*r_c + 1))/(zlambda ** 2)
  End If 
  !write(*,*)ig,zvclir (ifg)
else

  zvclir (ifg) = fourpi * zrhoir (ifg) / ((gpc(ig) ** 2) + (zlambda ** 2))

  !write(*,*)ig,zvclir (ifg)
endif

End Do
!stop
end subroutine


subroutine pseudocharge_rspace(lmax,npsd,qlm,zvclir,yukawa_in,zlambda,zilmt,zbessi)

  use modinput
use modinteg
  use mod_atoms, only: nspecies, natoms, idxas, natmtot,atposc,spr
  use mod_Gvector, only: ngrid,ngrtot
  use mod_muffin_tin, only: rmt, nrmtmax,nrmt
  use constants, only: zzero, fourpi, y00
  implicit none
  
  integer, intent(in) :: lmax
  integer, intent(in) :: npsd
  complex(8), intent(in) :: qlm((lmax+1)**2,natmtot)
  Complex (8), Intent (InOut) :: zvclir(ngrtot)
  logical, optional, intent(in) :: yukawa_in
  complex(dp),optional, intent(in) :: zlambda
  complex(dp),optional, intent(in) :: zilmt(0:,:)
  complex(dp),optional, intent(in) :: zbessi(:,0:,:)

  complex(dp) :: cor(ngrtot)
logical :: yukawa

  complex (8) :: zylm((lmax+1)**2), zrp((lmax+1)**2),alm((lmax+1)**2),zf1(nrmtmax),zf2(nrmtmax)
  complex (8) :: zt1,zt2, rr

  Real (8) :: a1(3),a2(3),a3(3),rv(3),tp(2),a1_abs,a2_abs,a3_abs
  Real (8) :: t1,t2,eps,rv_abs,ratom(3)
  real (8) :: bmat(3,3),binv(3,3), col(1,3)
  integer :: is, ia, ias, lm, l,ig,m
  integer :: i1,i2,i3,ir1,ir2,ir3,ir
  !external functions
  Real (8) :: factnm
  External factnm
  

if(present(yukawa_in))then 
  yukawa=yukawa_in
else
  yukawa=.false.
endif

  cor=zzero

  call zfftifc( 3, ngrid, 1, zvclir)
  
  ! write(*,*)"nspecies",nspecies, shape(nspecies)
  ! write(*,*)"natoms",natoms, shape(natoms)
  do is=1, nspecies
  do ia=1, natoms(is)
  ias=idxas(ia,is)
  
  ! write(*,*)"ngrtot",ngrtot
  ! write(*,*)"ias", ias
  ! write(*,*)"qlm",shape(qlm)
  ! write(*,*)"ngrtot?",shape(zvclir)
  ratom=atposc (:, ia, is)
  
  ! write(*,*)"ratom",ratom


  
  
  a1=input%structure%crystal%basevect(1, :)/ngrid(1)
  a2=input%structure%crystal%basevect(2, :)/ngrid(2)
  a3=input%structure%crystal%basevect(3, :)/ngrid(3)
  
  bmat(1,1:3)=a1
  bmat(2,1:3)=a2
  bmat(3,1:3)=a3
  Call r3minv (bmat, binv)
  
  col(1,:)=ratom
  col=matmul(col,binv)
  
  
  
  
  a1_abs=sqrt(a1(1)**2+a1(2)**2+a1(3)**2)
  a2_abs=sqrt(a2(1)**2+a2(2)**2+a2(3)**2)
  a3_abs=sqrt(a3(1)**2+a3(2)**2+a3(3)**2)
  
  ! write(*,*)"atoma poz",col(1,1),col(1,2),col(1,3)
  ! write(*,*)"mt izmērs",rmt(is)/a1_abs,rmt(is)/a2_abs,rmt(is)/a3_abs

    do l = 0, lmax
      t1 = factnm( 2*l + 2*npsd + 3, 2)/factnm( 2*l+1, 2)
      do m = -l, l
        lm = l*(l+1) + m + 1
        zrp (lm) = (qlm(lm, ias)) * t1
      end do
    end do

    do ir1=-floor(rmt(is)/a1_abs+col(1,1)),floor(rmt(is)/a1_abs+col(1,1))
        if (ir1.lt.0) then 
            i1=ngrid(1)+ir1+1
        else
            i1=ir1+1
        endif
        do ir2=-floor(rmt(is)/a2_abs+col(1,2)),floor(rmt(is)/a2_abs+col(1,2))
            if (ir2.lt.0) then
                i2=ngrid(1)+ir2+1 
            else 
                i2=ir2+1
            endif
            do ir3=-floor(rmt(is)/a3_abs+col(1,3)),floor(rmt(is)/a3_abs+col(1,3))
                if (ir3.lt.0) then
                    i3=ngrid(1)+ir3+1 
                else 
                    i3=ir3+1
                endif
                rv=ir1*a1+ir2*a2+ir3*a3-ratom
    
                rv_abs=dsqrt(rv(1)**2+rv(2)**2+rv(3)**2)
              
    
    
                if (rv_abs.le.rmt(is)) then
                    ig = (i3-1)*ngrid(2)*ngrid(1) + (i2-1)*ngrid(1) + i1
                    lm=0
                    
                    if (rv_abs.eq.0d0) then 
                      tp(1)=0d0
                    else
                      tp(1)=dacos(rv(3)/rv_abs) !theta
                    endif


                    if ((rv(1).eq.0d0).and.(rv(2).eq.0d0))then
                      tp(2)= 0d0 !to avoid NaN
                    else
                      tp(2)= sign(1d0,rv(2))*dacos(rv(1)/dsqrt(rv(1)**2+rv(2)**2)) !phi
                    endif
                    call genylm(lmax, tp, zylm)
  
            if(yukawa) then
              
              
if(.true.)then !skaitkiski
 
  lm=0
  Do l = 0, lmax
  zf1=zzero
  do ir=1, nrmt(is)
    rr=spr(ir,is)/rmt(is)
    zf1(ir)=(spr(ir,is)/rmt(is))**l*(1d0-(spr(ir,is)/rmt(is))**2)**npsd*zbessi(ir,l,is)*spr(ir,is)**2
  enddo

  call integ_cf (nrmt(is), is, zf1(1:nrmt(is)), zf2(1:nrmt(is)), mt_integw)
  zt1=zlambda**l/(factnm(2*l+1, 2)*zf2(nrmt(is)))

  Do m = - l, l
    lm = lm + 1
    alm(lm)=qlm(lm,ias)*zt1
  enddo
  enddo


  lm=0
  Do l = 0, lmax
    
    zt1=(rv_abs/rmt(is))**l*(1d0-(rv_abs/rmt(is))**2)**npsd
    if (rv_abs .Gt. input%structure%epslat) then
      Do m = - l, l
        lm = lm + 1
        cor(ig)= cor(ig) + zt1*alm(lm)*zylm(lm)
        zvclir(ig) = zvclir(ig) + zt1*alm(lm)*zylm(lm)
      enddo
    else
      cor(ig)= cor(ig) + zt1*alm(1)*y00
      zvclir(ig) = zvclir(ig) + zt1*alm(1)*y00
      exit
    endif
  enddo


else! analītiski

              Do l = 0, lmax
                t1=(rv_abs**2-rmt(is)**2)**npsd*rv_abs**l
                zt1=zlambda**(l+npsd+1)*t1
                t2=(-2d0)**npsd*factnm (npsd, 1)*rmt(is)**(l+npsd+2)* factnm(2*l+1, 2)
                zt2=zilmt(l+npsd+1,is)*t2
                zt1=zt1/zt2

                if (rv_abs .Gt. input%structure%epslat) then
                  Do m = - l, l
                    lm = lm + 1
                    cor(ig)= cor(ig)+zt1*zylm(lm)
                    zvclir(ig) = zvclir(ig)+zt1*zylm(lm)
                  enddo ! m
                else
                  cor(ig)=cor(ig)+zt1*y00
                  zvclir(ig) = zvclir(ig)+zt1*y00
                  exit
                endif
              enddo !l 
endif
            
            else !if not Yukawa
                    
                    Do l = 0,lmax
                      t1=(rv_abs/rmt(is))**l*(1d0-rv_abs**2/rmt(is)**2)**npsd/(rmt(is)**(l+3)*2**(npsd)*factnm (npsd, 1))
                      !if (rv_abs .Gt. 1d-6) then
                      if (rv_abs .Gt. input%structure%epslat) then
                          Do m = - l, l
                            lm = lm + 1
                            cor(ig)=cor(ig)+zrp(lm)*t1*zylm(lm)
                            zvclir(ig) = zvclir(ig)+zrp(lm)*t1*zylm(lm)
                          enddo ! m
                      else
                          cor(ig)=cor(ig)+zrp(1)*t1*y00
                          zvclir(ig) = zvclir(ig)+zrp(1)*t1*y00
                          exit   
                      endif
                    enddo ! l

            endif !if yukawa or not

                   
                endif 
            enddo !ir3 z
        enddo !ir2 y
    enddo !ir1 x

    enddo !ia
    enddo !is

    ! open(11,file='is_pseudo_r_cor_a.dat',status='replace')
    ! Do ig = 1, ngrtot
    !     write(11,*)(ig-1)*a1(1),",",dble(cor(ig))
    ! End Do
    ! close(11)
    



  !write(*,*)"pseudodensity_ir"
  
  
  call zfftifc( 3, ngrid, -1, zvclir)
  
  
  
  
  end subroutine

end module weinert
