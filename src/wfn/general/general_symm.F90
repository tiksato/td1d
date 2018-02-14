!######################################################################
subroutine general_symm(nfun, wfn)
!
! assumptions:
! 1) g, u, g, u, ... sequence
! 2) even ngrid
! 3) no population near the boundary [0,fd_ohalf) and (ngrid-fd_ohalf,ngrid]
!
  use const_mod, only : czero, runit
  use grid_mod, only : ngrid, fd_ohalf

  implicit none
  integer, intent(in) :: nfun
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 1:nfun)

  integer :: igrid, igl, igr, ifun
  complex(kind(0d0)) :: fac

  do ifun = 1, nfun

     fac = (-runit) ** (ifun - 1)

!     do igrid = 0, fd_ohalf - 1
!        igl = igrid
!        igr = ngrid - igrid
!        wfn(igl, ifun) = czero
!        wfn(igr, ifun) = czero
!     end do
!
!     do igrid = fd_ohalf, ngrid / 2 - 1
!        igl = igrid
!        igr = ngrid - igrid
!        wfn(igr, ifun) = wfn(igl, ifun) * fac
!     end do

     do igrid = 0, ngrid / 2 - 1
        igl = igrid
        igr = ngrid - igrid
        wfn(igr, ifun) = wfn(igl, ifun) * fac
     end do

  end do

end subroutine general_symm
!######################################################################
