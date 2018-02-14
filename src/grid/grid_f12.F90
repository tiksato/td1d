!################################################################################
subroutine grid_f12
  use root_mod, only : iprint, isoftr12
  use const_mod, only : zero, two
  use grid_mod, only : ngrid, x, v2, tcv2, dfdx, tc_type

  implicit none
  integer :: igrid, jgrid
  real(kind(0d0)), allocatable :: f12(:,:)

  allocate(f12(0:ngrid,2))

  call grid_f12_cal(f12)

  if(isoftr12 == 1) then
     stop 'f12 with isoftr12 = 1: nyi'
  else 
     if (trim(tc_type) == 'NOTC') then
        do igrid = 0, ngrid
           do jgrid = 0, ngrid
              ! scalar part
              tcv2(jgrid, igrid, 1) = v2(jgrid, igrid)
              tcv2(jgrid, igrid, 2) = v2(jgrid, igrid)
              ! vector part
              dfdx(jgrid, igrid) = zero
           end do
        end do
     else if (trim(tc_type) == 'ORTHOGONAL') then
        do igrid = 0, ngrid
           do jgrid = 0, ngrid
              ! scalar part
              tcv2(jgrid, igrid, 1) = v2(jgrid, igrid) - f12(abs(jgrid-igrid), 1) ** two - f12(abs(jgrid-igrid), 2)
              tcv2(jgrid, igrid, 2) = v2(jgrid, igrid) - f12(abs(jgrid-igrid), 1) ** two - f12(abs(jgrid-igrid), 2)
              ! vector part
              if (jgrid >= igrid) then
                 dfdx(jgrid, igrid) =  f12(abs(jgrid-igrid), 1)
              else
                 dfdx(jgrid, igrid) = -f12(abs(jgrid-igrid), 1)
              end if
           end do
        end do
     else if (trim(tc_type) == 'BIORTHOGONAL') then
        do igrid = 0, ngrid
           do jgrid = 0, ngrid
              tcv2(jgrid, igrid, 1) = v2(jgrid, igrid) - f12(abs(jgrid-igrid), 1) ** two - f12(abs(jgrid-igrid), 2)
              tcv2(jgrid, igrid, 2) = v2(jgrid, igrid) - f12(abs(jgrid-igrid), 1) ** two + f12(abs(jgrid-igrid), 2)
              ! vector part
              if (jgrid >= igrid) then
                 dfdx(jgrid, igrid) =  f12(abs(jgrid-igrid), 1)
              else
                 dfdx(jgrid, igrid) = -f12(abs(jgrid-igrid), 1)
              end if
           end do
        end do
     end if
  end if

  if (iprint > 1) then
     write(6, "('F12:')")
     do igrid = 0, ngrid
        write(6, "(i10,5f20.10)") igrid, x(igrid), v2(igrid,0), tcv2(igrid,0,1), tcv2(igrid,0,2), dfdx(igrid,0)
     end do
  end if

  deallocate(f12)

!debug
!stop 'debug: grid_f12'
!debug

end subroutine grid_f12
!################################################################################
!################################################################################
subroutine grid_f12_cal(f12)
  use root_mod, only : iprint
  use const_mod, only : one, two
  use grid_mod, only : ngrid, x, dgrid

  implicit none
  real(kind(0d0)), intent(out) :: f12(0:ngrid,2)

  integer :: ixx, nxx, igrid
  real(kind(0d0)) :: dxx, xx, fun, dfun

  nxx = ngrid * 1000
  dxx = dgrid / 1000

  do ixx = 0, nxx

     xx = dxx * ixx

     if (ixx == 0) then
        call grid_f12_bc(fun, dfun)
     else
        call grid_f12_fun(xx, fun, dfun)
     end if

!debug     if (iprint > 1) call grid_f12_print(ixx, xx, fun, dfun)
     if (mod(ixx,1000) == 0) then
        igrid = ixx / 1000
        f12(igrid, 1) = fun
        f12(igrid, 2) = dfun
     end if

     call grid_f12_rk4(xx, dxx, fun, dfun)

  end do


end subroutine grid_f12_cal
!################################################################################
!################################################################################
subroutine grid_f12_bc(fun, dfun)

  implicit none

  real(kind(0d0)), intent(out) :: fun
  real(kind(0d0)), intent(out) :: dfun

  real(kind(0d0)), parameter :: zero = 0.d+0, one = 1.d+0

  fun = zero
  dfun = one

!debug
!  fun = one
!  dfun = -one

end subroutine grid_f12_bc
!################################################################################
!################################################################################
subroutine grid_f12_fun(xx, fun, dfun)

  use root_mod, only : softr12
  use grid_mod, only : zeta

  implicit none

  real(kind(0d0)), intent(in) :: xx
  real(kind(0d0)), intent(in) :: fun
  real(kind(0d0)), intent(out) :: dfun

  real(kind(0d0)) :: xx2
  real(kind(0d0)), parameter :: zero = 0.d+0, one = 1.d+0, two = 2.d+0

  xx2 = xx * xx

  ! screened r12
  dfun = one / sqrt(xx2 + softr12) * exp(-zeta * xx2) - fun * fun
  ! bare r12
!  dfun = one / sqrt(xx**two + one) - fun * fun

end subroutine grid_f12_fun
!################################################################################
!################################################################################
subroutine grid_f12_rk4(xx, dxx, fun, dfun)

  implicit none

  real(kind(0d0)), intent(in) :: xx, dxx
  real(kind(0d0)), intent(inout) :: fun
  real(kind(0d0)), intent(inout) :: dfun

  real(kind(0d0)) :: x1, dx1, dx2, dx3, dx6, fun0, fun1
  real(kind(0d0)), parameter :: zero = 0.d+0, one = 1.d+0, two = 2.d+0, three = 3.d+0, six = 6.d+0, half = 5.d-1

  dx1 = dxx
  dx2 = dxx / two
  dx3 = dxx / three
  dx6 = dxx / six
  fun0 = fun

  ! step 1
  fun = fun + dfun * dx6
  fun1 = fun0 + dfun * dx2

  ! step 2
  x1 = xx + dx2
  call grid_f12_fun(x1, fun1, dfun)
  fun = fun + dfun * dx3
  fun1 = fun0 + dfun * dx2

  ! step 3
  x1 = xx + dx2
  call grid_f12_fun(x1, fun1, dfun)
  fun = fun + dfun * dx3
  fun1 = fun0 + dfun * dx1

  ! step 4
  x1 = xx + dx1
  call grid_f12_fun(x1, fun1, dfun)
  fun = fun + dfun * dx6

end subroutine grid_f12_rk4
!################################################################################
!################################################################################
subroutine grid_f12_print(ixx, xx, fun, dfun)

  implicit none

  integer, intent(in) :: ixx
  real(kind(0d0)), intent(in) :: xx
  real(kind(0d0)), intent(in) :: fun, dfun

  write(6, "(i5, 3f20.10)") ixx, xx, fun, dfun

end subroutine grid_f12_print
!################################################################################
