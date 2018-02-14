!################################################################################
module grid_mod

  implicit none

  integer :: ngrid
  real(kind(0d0)) :: dgrid
  real(kind(0d0)) :: x0
  real(kind(0d0)) :: rmax
  real(kind(0d0)) :: pmask
  real(kind(0d0)) :: xmask
  real(kind(0d0)) :: intmask
  real(kind(0d0)) :: yukawa
  integer :: llrmax, ulrmax

  logical :: domask
  logical :: docap
  character(len=16) :: mask_type
  character(len=16) :: norm_type
  real(kind(0d0)), allocatable :: x(:)
  real(kind(0d0)), allocatable :: mask(:)
  real(kind(0d0)), allocatable :: v1(:)
  real(kind(0d0)), allocatable :: gv1(:)
  real(kind(0d0)), allocatable :: v2(:,:)
  real(kind(0d0)), allocatable :: p(:)
  real(kind(0d0)), allocatable :: p2(:,:)
  real(kind(0d0)), allocatable :: efree(:)
  real(kind(0d0)), allocatable :: efree2(:,:)

  character(len=16) :: tc_type
  character(len=16) :: f12_type
  real(kind(0d0)) :: zeta
  real(kind(0d0)), allocatable :: tcv2(:,:,:)
  real(kind(0d0)), allocatable :: dfdx(:,:)

  real(kind(0d0)) :: rmax_atom
  integer, allocatable :: grid_inner(:)

  integer :: gll, gul
  integer :: fd_order
  integer :: fd_ohalf
  real(kind(0d0)) :: fd_fac
  real(kind(0d0)), allocatable :: fd_coeff1_in(:)
  real(kind(0d0)), allocatable :: fd_coeff1_l(:,:)
  real(kind(0d0)), allocatable :: fd_coeff1_r(:,:)
  real(kind(0d0)), allocatable :: fd_coeff_in(:)
  real(kind(0d0)), allocatable :: fd_coeff_l(:,:)
  real(kind(0d0)), allocatable :: fd_coeff_r(:,:)

  logical :: gbasis
  character(len=256) :: g09out
  integer :: ngbas,ngb2,ngeris
  real(kind(0d0)), allocatable :: govlp(:,:)
  real(kind(0d0)), allocatable :: goinv(:,:)
  real(kind(0d0)), allocatable :: gcore(:,:)
  real(kind(0d0)), allocatable :: geris(:)
  integer, allocatable :: gmap2(:,:)
  integer, allocatable :: gmap4(:,:)

end module grid_mod
!################################################################################
