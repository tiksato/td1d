!################################################################################
subroutine grid_print

  use io_mod, only : iostdo
  use grid_mod, only : ngrid, x, mask

  implicit none
  integer :: igrid

  do igrid = 0, ngrid
     write(iostdo, "(i10, 2f20.10)") igrid, x(igrid), mask(igrid)
  end do

end subroutine grid_print
!################################################################################
