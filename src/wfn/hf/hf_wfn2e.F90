!######################################################################
subroutine hf_wfn2e(wfn, wfn2e)

  use grid_mod, only : ngrid
  use wfn_mod, only : nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nspin)
  complex(kind(0d0)), intent(out) :: wfn2e(0:ngrid, 0:ngrid)
  integer :: igrid, jgrid

  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        wfn2e(jgrid, igrid) = wfn(jgrid, 1) * wfn(igrid, nspin)
     end do
  end do

end subroutine hf_wfn2e
!######################################################################
