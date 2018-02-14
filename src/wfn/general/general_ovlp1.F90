!######################################################################
complex(kind(0d0)) function general_ovlp1(rmax, wfnl, wfnr)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid)

  complex(kind(0d0)) :: tmp
  integer :: igrid, llgrid, ulgrid

!old  call get_irmax(rmax, llgrid, ulgrid, 1, ngrid-1)
  call get_irmax(rmax, llgrid, ulgrid, 0, ngrid)

  tmp = czero
  do igrid = llgrid, ulgrid
     tmp = tmp + conjg(wfnl(igrid)) * wfnr(igrid)
  end do
  general_ovlp1 = tmp * dgrid

  return

end function general_ovlp1
!######################################################################
