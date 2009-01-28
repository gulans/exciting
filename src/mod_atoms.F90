#include "maxdefinitions.inc"
module mod_atoms

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species

! maximum allowed atoms per species

! number of species
integer nspecies
! number of atoms for each species
integer natoms(_MAXSPECIES_)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot
! index to atoms and species
integer idxas(_MAXATOMS_,_MAXSPECIES_)
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! atomic positions in lattice coordinates
real(8) atposl(3,_MAXATOMS_,_MAXSPECIES_)
! atomic positions in Cartesian coordinates
real(8) atposc(3,_MAXATOMS_,_MAXSPECIES_)

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
character(256) sppath
! species filenames
character(256) spfname(_MAXSPECIES_)
! species name
character(256) spname(_MAXSPECIES_)
! species symbol
character(256) spsymb(_MAXSPECIES_)
! species nuclear charge
real(8) spzn(_MAXSPECIES_)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
logical ptnucl
! species electronic charge
real(8) spze(_MAXSPECIES_)
! species mass
real(8) spmass(_MAXSPECIES_)
! smallest radial point for each species
real(8) sprmin(_MAXSPECIES_)
! effective infinity for species
real(8) sprmax(_MAXSPECIES_)
! number of radial points to effective infinity for each species
integer spnr(_MAXSPECIES_)
! maximum spnr over all the species
integer spnrmax
! maximum allowed states for each species
integer, parameter :: maxspst=40
! number of states for each species
integer spnst(_MAXSPECIES_)
! maximum spnst over all the species
integer spnstmax
! state principle quantum number for each species
integer spn(maxspst,_MAXSPECIES_)
! state l value for each species
integer spl(maxspst,_MAXSPECIES_)
! state k value for each species
integer spk(maxspst,_MAXSPECIES_)
! spcore is .true. if species state is core
logical spcore(maxspst,_MAXSPECIES_)
! state eigenvalue for each species
real(8) speval(maxspst,_MAXSPECIES_)
! state occupancy for each species
real(8) spocc(maxspst,_MAXSPECIES_)
! species radial mesh
real(8), allocatable :: spr(:,:)
! species charge density
real(8), allocatable :: sprho(:,:)
! species self-consistent potential
real(8), allocatable :: spvr(:,:)

end module