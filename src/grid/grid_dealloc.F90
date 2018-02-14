!################################################################################
subroutine grid_dealloc

  use grid_mod, only : x, mask, v1, gv1, v2, p, p2, efree, efree2, tcv2, dfdx, &
       & fd_coeff1_in, fd_coeff1_l, fd_coeff1_r, fd_coeff_in, fd_coeff_l, fd_coeff_r, &
       & grid_inner

  implicit none

  deallocate(grid_inner)
  deallocate(fd_coeff_r)
  deallocate(fd_coeff_l)
  deallocate(fd_coeff_in)
  deallocate(fd_coeff1_r)
  deallocate(fd_coeff1_l)
  deallocate(fd_coeff1_in)
  deallocate(dfdx)
  deallocate(tcv2)
  deallocate(efree2)
  deallocate(efree)
  deallocate(p2)
  deallocate(p)
  deallocate(v2)
  deallocate(gv1)
  deallocate(v1)
  deallocate(mask)
  deallocate(x)

end subroutine grid_dealloc
!################################################################################
