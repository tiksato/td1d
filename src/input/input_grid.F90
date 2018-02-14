!################################################################################
subroutine input_grid(ioin)

  use root_mod, only : ical
  use const_mod, only : zero, one, half
  use grid_mod, only : ngrid, dgrid, x0, rmax, pmask, xmask, intmask, domask, yukawa, &
       & docap, mask_type, norm_type, fd_order, fd_ohalf, zeta, tc_type, f12_type, &
       & gll, gul, rmax_atom, gbasis, ngbas, g09out

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr
  namelist /grid/ ngrid, dgrid, x0, rmax, pmask, xmask, intmask, domask, yukawa, docap, &
       & mask_type, norm_type, fd_order, zeta, tc_type, f12_type, rmax_atom, gbasis, ngbas, g09out

  ngrid = -1               ! number of grid points
  dgrid = -one             ! grid spacing
  rmax = -one              ! radius defining core region
  domask = .true.          ! mask after each propagation
  docap = .false.          ! cap included in 1e potential
  yukawa = 0.25d+0         ! exponential of yukawa potential
  pmask = 0.85             ! mask boundary (percentage to x0)
  intmask = -one           ! intensity of mask (for QUDRATIC mask)
  mask_type = 'COS4'  ! QUADRATIC, COS4, COS8
  norm_type = 'FORMAL'     ! FORMAL, DIRECT
  fd_order = 8             ! finite difference order for kinetic operator
  zeta = one               ! exponent of weight gaussian for tcham
  tc_type = 'BIORTHOGONAL' ! NOTC, ORTHOGONAL
  f12_type = 'GAUSS'       ! BARE, GAUSS2
  rmax_atom = 10.d+0
  gbasis = .false.         ! Use gaussian basis instead of 1D finite difference
  ngbas = -1               ! Number of gaussian basis function
  g09out = ""              ! Name of the gaussian output file.

  rewind(ioin)
  read(unit=ioin, nml=grid, iostat=ioerr)
  if(ioerr /= 0) stop "error in namelist grid."

  call util_transchar(mask_type)
  call util_transchar(norm_type)
  call util_transchar(tc_type)
  call util_transchar(f12_type)

  if(ngrid < 0) stop "error in grid.ngrid."
  if(dgrid < zero) stop "error in grid.dgrid."
  if(ical >= 0 .and. rmax < zero)  stop "error in grid.rmax."
  if(ical >= 0 .and. domask .and. docap) stop "error in grid.domask or grid.docap."
  if(ical >= 0 .and. intmask < zero .and. trim(mask_type) == 'QUADRATIC') stop "error in grid.intmask and grid.mask_type."
  if(fd_order < 2 .or. mod(fd_order,2) /= 0) stop "error in grid.fd_order."
  if(gbasis .and. ngbas < 0) stop "error in grid.ngbas."
  if(gbasis .and. trim(g09out)=="") stop "error in grid.g09out."
  if(gbasis) ngrid = ngbas-1
  if(gbasis) dgrid = 1

  x0 = half * dgrid * ngrid
  xmask = x0 * pmask
  fd_ohalf = fd_order / 2
!bug  gll = fd_ohalf
!bug  gul = ngrid - fd_ohalf
  gll = 0
  gul = ngrid

  write(6, nml=grid)

end subroutine input_grid
!################################################################################
