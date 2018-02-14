!######################################################################
subroutine general_ovlp(inner, rmax, wfnl, wfnr, slr)
  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid, gll, gul
  use wfn_mod, only : nfun

  implicit none
  logical, intent(in) :: inner
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(out) :: slr(1:nfun, 1:nfun)
  complex(kind(0d0)) :: tmp
  integer :: igrid, ifun, jfun, llgrid, ulgrid

  slr(1:nfun, 1:nfun) = czero
  call get_irmax(rmax, llgrid, ulgrid, gll, gul)

  if (inner) then
     do ifun = 1, nfun
        do jfun = 1, nfun
           tmp = czero
           do igrid = llgrid, ulgrid
              tmp = tmp + conjg(wfnl(igrid, ifun)) * wfnr(igrid, jfun)
           end do
           slr(ifun, jfun) = tmp * dgrid
        end do
     end do
  else
     do ifun = 1, nfun
        do jfun = 1, nfun
           tmp = czero
           do igrid = 1, llgrid - 1
              tmp = tmp + conjg(wfnl(igrid, ifun)) * wfnr(igrid, jfun)
           end do
           do igrid = ulgrid + 1, ngrid - 1
              tmp = tmp + conjg(wfnl(igrid, ifun)) * wfnr(igrid, jfun)
           end do
           slr(ifun, jfun) = tmp * dgrid
        end do
     end do
  end if

end subroutine general_ovlp
!######################################################################
