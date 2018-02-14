!######################################################################
integer function x2e_size()

  use grid_mod, only : ngrid

  implicit none

  x2e_size = (ngrid + 1) * (ngrid + 1)

end function x2e_size
!######################################################################
