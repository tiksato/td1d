!################################################################################
subroutine grid_init

  use wfn_mod, only : tcham
  use grid_mod, only : gbasis, g09out

  implicit none

  call grid_alloc
  call grid_mkgrid
  call grid_mkmask
  call grid_v1
  call grid_v2
  if (tcham) call grid_f12
  call grid_efree
  call grid_fdcoeff1
  call grid_fdcoeff2
  call grid_mkinner
  if (gbasis) call gbasis_init(g09out)
  !stop "FOR DEBUG @ grid_init!"
end subroutine grid_init
!################################################################################
