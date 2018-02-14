!######################################################################
subroutine x2e_readhf(wfnin, wfn)

  use const_mod, only : czero
  use grid_mod, only : ngrid
  use wfn_mod, only : nspin
  use mol_mod, only : mult

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: wfnin(0:ngrid, 1:*)
  complex(kind(0d0)), intent(out) :: wfn(0:ngrid, 0:ngrid)
  !--------------------------------------------------------------------
  integer :: igrid, jgrid
  complex(kind(0d0)) :: tmp

  call x2e_clear(wfn)

  if (mult == 0) then
     do igrid = 0, ngrid
        do jgrid = 0, igrid - 1
           tmp = wfnin(jgrid, 1) * wfnin(igrid, nspin)
           wfn(jgrid, igrid) = tmp
           wfn(igrid, jgrid) = tmp
        end do
        wfn(igrid, igrid) = wfnin(igrid, 1) * wfnin(igrid, nspin)
     end do
  else
     do igrid = 0, ngrid
        do jgrid = 0, igrid - 1
           tmp = wfnin(jgrid, 1) * wfnin(igrid, 2) - wfnin(jgrid, 2) * wfnin(igrid, 1)
           wfn(jgrid, igrid) = + tmp
           wfn(igrid, jgrid) = - tmp
        end do
     end do
  end if

  call x2e_ort(wfn)

end subroutine x2e_readhf
!######################################################################
