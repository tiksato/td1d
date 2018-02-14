!######################################################################
subroutine x2e_den2e(wfn, den2e)

  use grid_mod, only : ngrid

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: den2e(0:ngrid, 0:ngrid)
  integer :: igrid, jgrid

  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        den2e(jgrid, igrid) = conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do

end subroutine x2e_den2e
!######################################################################
