!################################################################################
subroutine grid_mkgrid

  use grid_mod, only : ngrid, x0, dgrid, x

  implicit none
  integer :: igrid

  do igrid = 0, ngrid
     x(igrid) = -x0 + dgrid * igrid
  end do

end subroutine grid_mkgrid
!################################################################################
