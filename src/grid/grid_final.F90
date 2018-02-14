!################################################################################
subroutine grid_final

  use grid_mod, only : gbasis
  implicit none

  call grid_dealloc
  if (gbasis) call gbasis_final

end subroutine grid_final
!################################################################################
