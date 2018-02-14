program f12

  implicit none

  character(len=16) :: carg
  integer :: nx
  real(kind(0d0)) :: dx
  real(kind(0d0)) :: fun
  real(kind(0d0)) :: dfun
  integer :: type
  real(kind(0d0)) :: zeta

  integer :: ix
  real(kind(0d0)) :: x

  call getarg(1, carg); read(carg, *) nx
  call getarg(2, carg); read(carg, *) dx
  call getarg(3, carg); read(carg, *) type
  call getarg(4, carg); read(carg, *) zeta

  do ix = 0, nx

     x = dx * ix

     if (ix == 0) then
        call f12_bound(type, fun, dfun)
     else
        call f12_fun(type, zeta, x, fun, dfun)
     end if

     call f12_print(type, zeta, ix, x, fun, dfun)
     call f12_rk4(type, zeta, x, dx, fun, dfun)

  end do

end program f12
!###############
subroutine f12_bound(type, fun, dfun)

  implicit none

  integer, intent(in) :: type
  real(kind(0d0)), intent(out) :: fun
  real(kind(0d0)), intent(out) :: dfun

  real(kind(0d0)), parameter :: zero = 0.d+0, one = 1.d+0

  fun = zero
  dfun = one

!debug
!  fun = one
!  dfun = -one

end subroutine f12_bound
!###############
subroutine f12_fun(type, zeta, x, fun, dfun)

  implicit none

  integer, intent(in) :: type
  real(kind(0d0)), intent(in) :: zeta
  real(kind(0d0)), intent(in) :: x
  real(kind(0d0)), intent(in) :: fun
  real(kind(0d0)), intent(out) :: dfun

  real(kind(0d0)), parameter :: zero = 0.d+0, one = 1.d+0, two = 2.d+0

!sin
!  dfun = sqrt(one - fun * fun)
!exp
!  dfun = -fun
  if (type == 0) then
     !gauss-screened r12
     dfun = one / sqrt(x**two + one) * exp(-zeta * x**two) - fun * fun
  else if (type == 1) then
     !gauss-screened r12 with gauss-screened f12
     dfun = one / sqrt(x**two + one) - exp(-zeta * x**two) * fun * fun + two * zeta * x * fun
  else if (type == 2) then
     !slater-screened r12 with slater-screened f12
     dfun = one / sqrt(x**two + one) - exp(-zeta * x) * fun * fun + zeta * fun
  end if
!bare r12
!  dfun = one / sqrt(x**two + one) - fun * fun

end subroutine f12_fun
!###############
subroutine f12_rk4(type, zeta, x, dx, fun, dfun)

  implicit none

  integer, intent(in) :: type
  real(kind(0d0)), intent(in) :: zeta
  real(kind(0d0)), intent(in) :: x, dx
  real(kind(0d0)), intent(inout) :: fun
  real(kind(0d0)), intent(inout) :: dfun

  real(kind(0d0)) :: x1, dx1, dx2, dx3, dx6, fun0, fun1
  real(kind(0d0)), parameter :: zero = 0.d+0, one = 1.d+0, two = 2.d+0, three = 3.d+0, six = 6.d+0, half = 5.d-1

  dx1 = dx
  dx2 = dx / two
  dx3 = dx / three
  dx6 = dx / six
  fun0 = fun

  ! step 1
  fun = fun + dfun * dx6
  fun1 = fun0 + dfun * dx2

  ! step 2
  x1 = x + dx2
  call f12_fun(type, zeta, x1, fun1, dfun)
  fun = fun + dfun * dx3
  fun1 = fun0 + dfun * dx2

  ! step 3
  x1 = x + dx2
  call f12_fun(type, zeta, x1, fun1, dfun)
  fun = fun + dfun * dx3
  fun1 = fun0 + dfun * dx1

  ! step 3
  x1 = x + dx1
  call f12_fun(type, zeta, x1, fun1, dfun)
  fun = fun + dfun * dx6

end subroutine f12_rk4
!###############
subroutine f12_print(type, zeta, ix, x, fun, dfun)

  implicit none

  integer, intent(in) :: type
  real(kind(0d0)), intent(in) :: zeta
  integer, intent(in) :: ix
  real(kind(0d0)), intent(in) :: x
  real(kind(0d0)), intent(in) :: fun, dfun

  real(kind(0d0)) :: coul, scoul, ecoul, fun1, dfun1
  real(kind(0d0)), parameter :: one = 1.d+0, two = 2.d+0

  if (type == 0) then
     fun1 = fun
     dfun1 = dfun
     coul = one / sqrt(x**two + one)
     scoul = one / sqrt(x**two + one) * exp(-zeta * x ** two)
     ecoul = coul - fun1 ** two - dfun1
  else if (type == 1) then
     fun1 = exp(-zeta * x ** two) * fun
     dfun1 = exp(-zeta * x ** two) * ( -two * zeta * x * fun + dfun)
     coul = one / sqrt(x**two + one)
     scoul = one / sqrt(x**two + one) * exp(-zeta * x ** two)
     ecoul = coul - fun1 ** two - dfun1
  else if (type == 2) then
     fun1 = exp(-zeta * x) * fun
     dfun1 = exp(-zeta * x) * ( -zeta * fun + dfun)
     coul = one / sqrt(x**two + one)
     scoul = one / sqrt(x**two + one) * exp(-zeta * x)
     ecoul = coul - fun1 ** two - dfun1
  end if


  write(6, "(i5, 6f20.10)") ix, x, coul, scoul, fun1*fun1, dfun1, ecoul

end subroutine f12_print
