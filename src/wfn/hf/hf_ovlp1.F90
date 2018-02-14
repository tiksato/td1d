!######################################################################
complex(kind(0d0)) function hf_ovlp1(rmax, wfnl, wfnr)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid, gll, gul

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid)

  complex(kind(0d0)) :: tmp
  integer :: igrid, llgrid, ulgrid

  call get_irmax(rmax, llgrid, ulgrid, gll, gul)

  tmp = czero
  do igrid = llgrid, ulgrid
     tmp = tmp + conjg(wfnl(igrid)) * wfnr(igrid)
  end do
  hf_ovlp1 = tmp * dgrid

  return

end function hf_ovlp1
!######################################################################
complex(kind(0d0)) function hf_ovlp1_gbasis(wfnl, wfnr)

  use const_mod, only : czero
  use grid_mod, only : ngbas

  implicit none
  complex(kind(0d0)), intent(in) :: wfnl(1:ngbas)
  complex(kind(0d0)), intent(in) :: wfnr(1:ngbas)

  complex(kind(0d0)) :: tmp
  integer :: mu

  tmp = czero
  do mu = 1, ngbas
     tmp = tmp + conjg(wfnl(mu)) * wfnr(mu)
  end do
  hf_ovlp1_gbasis = tmp

  return

end function hf_ovlp1_gbasis
!######################################################################
