!######################################################################
subroutine x2e_print(iunit, wfn)

  use grid_mod, only : ngrid, x

  implicit none
  integer, intent(in) :: iunit
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid

  write(iunit, "('x2e wavefunction:', i10)") ngrid
  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        write(iunit, "(4F20.10)") x(igrid), x(jgrid), wfn(igrid, jgrid)
     end do
     ! for gnuplot
     write(iunit, *)
  end do

end subroutine x2e_print
!######################################################################
