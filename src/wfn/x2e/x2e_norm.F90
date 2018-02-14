!######################################################################
real(kind(0d0)) function x2e_norm(rmax, wfn, p0, p1, p2)

  use const_mod, only : zero, two, czero
  use grid_mod, only : ngrid, dgrid

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)
  real(kind(0d0)), intent(out) :: p0, p1, p2

  integer :: igrid, jgrid, llgrid, ulgrid
  complex(kind(0d0)) :: tmp

  p0 = zero
  p1 = zero
  p2 = zero
!old  call get_irmax(rmax, llgrid, ulgrid, 1, ngrid-1)
  call get_irmax(rmax, llgrid, ulgrid, 0, ngrid)

  tmp = czero
  do igrid = llgrid, ulgrid
     tmp = tmp + conjg(wfn(igrid, igrid)) * wfn(igrid, igrid)
  end do
  p0 = p0 + real(tmp)

  tmp = czero
  do igrid = llgrid, ulgrid
     do jgrid = llgrid, igrid - 1
        tmp = tmp + conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do
  p0 = p0 + real(tmp) * two

  tmp = czero
  do igrid = llgrid, ulgrid
     do jgrid = 1, llgrid - 1
        tmp = tmp + conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do
  do igrid = ulgrid + 1, ngrid - 1
     do jgrid = llgrid, ulgrid
        tmp = tmp + conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do
  p1 = p1 + real(tmp) * two

  tmp = czero
  do igrid = 1, llgrid - 1
     tmp = tmp + conjg(wfn(igrid, igrid)) * wfn(igrid, igrid)
  end do
  do igrid = ulgrid + 1, ngrid - 1
     tmp = tmp + conjg(wfn(igrid, igrid)) * wfn(igrid, igrid)
  end do
  p2 = p2 + real(tmp)

  tmp = czero
  do igrid = 1, llgrid - 1
     do jgrid = 1, igrid - 1
        tmp = tmp + conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do
  do igrid = ulgrid + 1, ngrid - 1
     do jgrid = ulgrid + 1, igrid - 1
        tmp = tmp + conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do
  do igrid = ulgrid + 1, ngrid - 1
     do jgrid = 1, llgrid - 1
        tmp = tmp + conjg(wfn(jgrid, igrid)) * wfn(jgrid, igrid)
     end do
  end do
  p2 = p2 + real(tmp) * two

  p0 = p0 * dgrid ** two
  p1 = p1 * dgrid ** two
  p2 = p2 * dgrid ** two

  x2e_norm = p0 + p1 + p2

end function x2e_norm
!######################################################################
