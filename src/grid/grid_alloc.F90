!################################################################################
subroutine grid_alloc

  use grid_mod, only : ngrid, x, mask, v1, gv1, v2, p, p2, efree, efree2, tcv2, dfdx, &
                     & fd_order, fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r, &
                     & fd_coeff_in, fd_coeff_l, fd_coeff_r, grid_inner

  implicit none

  allocate(x     (0:ngrid))
  allocate(mask  (0:ngrid))
  allocate(v1    (0:ngrid))
  allocate(gv1   (0:ngrid))
  allocate(v2    (0:ngrid, 0:ngrid))
  allocate(p     (0:ngrid-1))
  allocate(p2    (0:ngrid-1, 0:ngrid-1))
  allocate(efree (0:ngrid-1))
  allocate(efree2(0:ngrid-1, 0:ngrid-1))
  allocate(tcv2  (0:ngrid, 0:ngrid, 1:2))
  allocate(dfdx  (0:ngrid, 0:ngrid))

  allocate(fd_coeff1_in(-fd_ohalf:fd_ohalf))
  allocate(fd_coeff1_l(0:(fd_order + 1), 0:fd_ohalf))
  allocate(fd_coeff1_r(0:(fd_order + 1), 0:fd_ohalf))
  allocate(fd_coeff_in(-fd_ohalf:fd_ohalf))
  allocate(fd_coeff_l(0:(fd_order + 1), 0:fd_ohalf))
  allocate(fd_coeff_r(0:(fd_order + 1), 0:fd_ohalf))
  allocate(grid_inner(0:ngrid))

end subroutine grid_alloc
!################################################################################
