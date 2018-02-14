!################################################################################
subroutine util_proj(npt0, npt1, nval, y0, y1)

  use const_mod, only : zero, two, czero, runit, half

  implicit none
  integer, intent(in) :: npt0, npt1, nval
  complex(kind(0d0)), intent(in)  :: y0(0:npt0, 1:nval)
  complex(kind(0d0)), intent(out) :: y1(0:npt1, 1:nval)
  integer :: ipt0, ipt1, ii, ival, ratio
  complex(kind(0d0)) :: facl, facr, facd

  ratio = npt1 / npt0
  write(6, "('util_proj: ratio = ', i5)") ratio

  do ipt0 = 0, npt0 - 1
     do ii = 0, ratio - 1
        ipt1 = ipt0 * ratio + ii
        facl = runit * (ratio - ii)
        facr = runit * ii
        facd = facl + facr
        do ival = 1, nval
           y1(ipt1, ival) = (y0(ipt0, ival) * facl + y0(ipt0 + 1, ival) * facr) / facd
        end do
     end do
  end do

  do ival = 1, nval
     y1(npt1, ival) = y0(npt0, ival)
  end do

end subroutine util_proj
!################################################################################
!################################################################################
subroutine util_proj_old(npt0, npt1, nval, x0, y0, x1, y1)

  use const_mod, only : zero, two, czero, half

  implicit none
  integer, intent(in) :: npt0, npt1, nval
  real(kind(0d0)), intent(in) :: x0(1:npt0)
  complex(kind(0d0)), intent(in) :: y0(1:npt0, 1:nval)
  real(kind(0d0)), intent(in) :: x1(1:npt1)
  complex(kind(0d0)), intent(out) :: y1(1:npt1, 1:nval)
  integer :: ipt0, ipt1
  real(kind(0d0)) :: diff_l, diff_r, denom
  complex(kind(0d0)), allocatable :: numer(:)
  real(kind(0d0)), parameter :: thresh = 1.D-6

  allocate(numer(1:nval))
  y1(1:npt1, 1:nval) = czero

  do ipt1 = 1, npt1
     do ipt0 = 1, npt0 - 1
        diff_l = x1(ipt1) - x0(ipt0)
        diff_r = x1(ipt1) - x0(ipt0 + 1)
        if (abs(diff_l) < thresh) then
           y1(ipt1, 1:nval) = y0(ipt0, 1:nval)
           exit
        else if (abs(diff_r) < thresh) then
           y1(ipt1, 1:nval) = y0(ipt0 + 1, 1:nval)
           exit
        else if (diff_l * diff_r < zero) then
           numer(1:nval) = y0(ipt0, 1:nval) * abs(diff_r) + y0(ipt0 + 1, 1:nval) * abs(diff_l)
           denom = abs(diff_l) + abs(diff_r)
           y1(ipt1, 1:nval) = numer(1:nval) / denom
           exit
        end if
     end do
  end do

  deallocate(numer)

end subroutine util_proj_old
!################################################################################
!################################################################################
subroutine util_proj_real(npt0, npt1, nval, x0, y0, x1, y1)

  use const_mod, only : zero, two, half

  implicit none
  integer, intent(in) :: npt0, npt1, nval
  real(kind(0d0)), intent(in) :: x0(1:npt0)
  real(kind(0d0)), intent(in) :: y0(1:npt0, 1:nval)
  real(kind(0d0)), intent(in) :: x1(1:npt1)
  real(kind(0d0)), intent(out) :: y1(1:npt1, 1:nval)
  integer :: ipt0, ipt1
  real(kind(0d0)) :: diff_l, diff_r, denom
  real(kind(0d0)), allocatable :: numer(:)
  real(kind(0d0)), parameter :: thresh = 1.D-6

  allocate(numer(1:nval))
  y1(1:npt1, 1:nval) = zero

  do ipt1 = 1, npt1
     do ipt0 = 1, npt0 - 1
        diff_l = x1(ipt1) - x0(ipt0)
        diff_r = x1(ipt1) - x0(ipt0 + 1)
        if (abs(diff_l) < thresh) then
           y1(ipt1, 1:nval) = y0(ipt0, 1:nval)
           exit
        else if (abs(diff_r) < thresh) then
           y1(ipt1, 1:nval) = y0(ipt0 + 1, 1:nval)
           exit
        else if (diff_l * diff_r < zero) then
           numer(1:nval) = y0(ipt0, 1:nval) * abs(diff_r) + y0(ipt0 + 1, 1:nval) * abs(diff_l)
           denom = abs(diff_l) + abs(diff_r)
           y1(ipt1, 1:nval) = numer(1:nval) / denom
           exit
        end if
     end do
  end do

  deallocate(numer)

end subroutine util_proj_real
!################################################################################
