!######################################################################
subroutine hf_den2e(wfn, den2e)

  use grid_mod, only : ngrid
  use wfn_mod, only : nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nspin)
  complex(kind(0d0)), intent(out) :: den2e(0:ngrid, 0:ngrid)
  integer :: igrid, jgrid
  complex(kind(0d0)) :: tmp

  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        tmp = wfn(jgrid, 1) * wfn(igrid, nspin)
        den2e(jgrid, igrid) = conjg(tmp) * tmp
     end do
  end do

end subroutine hf_den2e
!######################################################################
