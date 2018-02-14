!######################################################################
complex(kind(0d0)) function x2e_iprod(rmax, wfnl, wfnr)

  use omp_mod
  use const_mod, only : two, czero
  use grid_mod, only : ngrid, dgrid

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid, llgrid, ulgrid
  complex(kind(0d0)) :: val, tmp

  x2e_iprod = czero
!old  call get_irmax(rmax, llgrid, ulgrid, 1, ngrid-1)
  call get_irmax(rmax, llgrid, ulgrid, 0, ngrid)

!$omp parallel default(shared) private(igrid, jgrid, val, tmp) reduction(+:x2e_iprod)

  call omp_mod_thread(llgrid, ulgrid)

  val = czero
  tmp = czero
  do igrid = ng0, ng1
     tmp = tmp + conjg(wfnl(igrid, igrid)) * wfnr(igrid, igrid)
  end do
  val = val + tmp

  tmp = czero
  do igrid = ng0, ng1
     do jgrid = llgrid, igrid - 1
        tmp = tmp + conjg(wfnl(jgrid, igrid)) * wfnr(jgrid, igrid)
     end do
  end do
  val = val + tmp * two

  val = val * dgrid ** two
  x2e_iprod = x2e_iprod + val

!$omp end parallel
  
end function x2e_iprod
!######################################################################
