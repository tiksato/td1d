!################################################################################
subroutine general_hprod1(time, wfn, hwfn)

  use omp_mod
  use root_mod, only : icomp
  use field_mod, only : gauge
  use const_mod, only : zero, czero, iunit, one
  use grid_mod, only : ngrid, gll, gul, x, v1
  use grid_mod, only : fd_order, fd_ohalf, fd_coeff_in, fd_coeff1_in

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid)
  complex(kind(0d0)), intent(out) :: hwfn(0:ngrid)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: lfield
  real(kind(0d0)), external :: field
  integer :: gridll, gridul, k, kmax, igrid
  complex(kind(0d0)) :: tmp, coeff_in(-fd_ohalf:fd_ohalf)

  hwfn = czero
  lfield = field(time)
  kmax = 2 * fd_ohalf + 1
!  gridll = max(gll, fd_ohalf)
!  gridul = min(gul, ngrid-fd_ohalf)
  gridll = gll
  gridul = gul

  if (icomp == 0 .or. trim(gauge) == 'L') then
     coeff_in(-fd_ohalf:fd_ohalf) = fd_coeff_in(-fd_ohalf:fd_ohalf)

    !$omp parallel default(shared) private(tmp)
    !$omp do
     do igrid = gridll, gridul
        tmp = wfn(igrid) * (coeff_in(0) + v1(igrid) - lfield * x(igrid))
!!!        do k = 1, fd_ohalf
!!!           tmp = tmp + wfn(igrid - k) * coeff_in(-k) &
!!!                   & + wfn(igrid + k) * coeff_in(+k)
!!!        end do
        do k = 1, fd_ohalf
           if (igrid - k >= 0) then
              tmp = tmp + wfn(igrid - k) * coeff_in(-k)
           end if
           if (igrid + k <= ngrid) then
              tmp = tmp + wfn(igrid + k) * coeff_in(+k)
           end if
        end do
        hwfn(igrid) = hwfn(igrid) + tmp
     end do
     !$omp end do
     !$omp end parallel
  else
     coeff_in(-fd_ohalf:fd_ohalf) = fd_coeff_in (-fd_ohalf:fd_ohalf) &
               & - iunit * lfield * fd_coeff1_in(-fd_ohalf:fd_ohalf)

     !$omp parallel default(shared) private(tmp)
     !$omp do
     do igrid = gridll, gridul
        tmp = wfn(igrid) * (coeff_in(0) + v1(igrid))
!!!        do k = 1, fd_ohalf
!!!           tmp = tmp + wfn(igrid - k) * coeff_in(-k) &
!!!                   & + wfn(igrid + k) * coeff_in(+k)
!!!        end do
        do k = 1, fd_ohalf
           if (igrid - k >= 0) then
              tmp = tmp + wfn(igrid - k) * coeff_in(-k)
           end if
           if (igrid + k <= ngrid) then
              tmp = tmp + wfn(igrid + k) * coeff_in(+k)
           end if
        end do
        hwfn(igrid) = hwfn(igrid) + tmp
     end do
     !$omp end do
     !$omp end parallel
  end if

end subroutine general_hprod1
!################################################################################
