!######################################################################
subroutine x2e_write(iunit, wfn)

  use grid_mod, only : ngrid, x

  implicit none
  integer, intent(in) :: iunit
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)

  write(iunit, "( I25)") ngrid
  write(iunit, "( E25.15)") x(0:ngrid)
  write(iunit, "(2E25.15)") wfn(0:ngrid, 0:ngrid)

end subroutine x2e_write
!######################################################################
