!######################################################################
real(kind(0d0)) function x2e_ene(wfnl, wfnr, ng0, ng1)

  use const_mod, only : two, czero
  use grid_mod, only : ngrid, dgrid

  implicit none
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid, 0:ngrid)
  integer, intent(in) :: ng0, ng1
  integer :: igrid, jgrid
  complex(kind(0d0)) :: val, tmp

  val = czero

  tmp = czero
  do igrid = ng0, ng1
     tmp = tmp + conjg(wfnl(igrid, igrid)) * wfnr(igrid, igrid)
  end do
  val = val + tmp

  tmp = czero
  do igrid = ng0, ng1
     do jgrid = 1, igrid - 1
        tmp = tmp + conjg(wfnl(jgrid, igrid)) * wfnr(jgrid, igrid)
     end do
  end do
  val = val + tmp * two

  val = val * dgrid ** two
  x2e_ene = real(val)

  return

end function x2e_ene
!######################################################################
