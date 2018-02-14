!######################################################################
subroutine util_zcopy(n, x, incx, y, incy)

  use omp_mod
  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n, incx, incy
  complex(kind(0d0)), intent(in) :: x(1:*)
  complex(kind(0d0)), intent(out) :: y(1:n)

  integer :: i

  if (incx <  0) stop 'util_zcopy: bad incx.'
  if (incy <= 0) stop 'util_zcopy: bad incy.'

  if (incx == 0) then
     do i = 1, n
        y(i) = x(1)
     end do
#ifdef USE_BLAS
  else
     call zcopy(n, x, incx, y, incy)
  end if
#else
  else
     do i = 1, n
        y(i) = x(i)
     end do
  end if
#endif

  return

end subroutine util_zcopy
!######################################################################
!######################################################################
subroutine util_zscal(n, a, x, incx)

  use omp_mod
  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n, incx
  complex(kind(0d0)), intent(in) :: a
  complex(kind(0d0)), intent(inout) :: x(1:*)

  integer :: i

  if (incx <= 0) stop 'util_zscal: bad incx.'

#ifdef USE_BLAS
  call zscal(n, a, x, incx)
#else
  do i = 1, n
     x(i) = a * x(i)
  end do
#endif

  return

end subroutine util_zscal
!######################################################################
!######################################################################
subroutine util_zaxpy(n, a, x, incx, y, incy)

  use omp_mod
  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n, incx, incy
  complex(kind(0d0)), intent(in) :: a
  complex(kind(0d0)), intent(in) :: x(1:*)
  complex(kind(0d0)), intent(inout) :: y(1:n)

  integer :: i

  if (incx <  0) stop 'util_zaxpy: bad incx.'
  if (incy <= 0) stop 'util_zaxpy: bad incy.'

  if (incx == 0) then
     do i = 1, n
        y(i) = y(i) + a * x(1)
     end do
#ifdef USE_BLAS
  else
     call zaxpy(n, a, x, incx, y, incy)
  end if
#else
  else
     do i = 1, n
        y(i) = y(i) + a * x(i)
     end do
  end if
#endif

  return

end subroutine util_zaxpy
!######################################################################
complex(kind(0d0)) function util_zdotc(n, x, incx, y, incy)

  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n, incx, incy
  complex(kind(0d0)), intent(in) :: x(1:n)
  complex(kind(0d0)), intent(in) :: y(1:n)

  integer :: i
  complex(kind(0d0)) :: tmp

  if (incx /= 1) stop 'util_zdotc: bad incx.'
  if (incy /= 1) stop 'util_zdotc: bad incy.'

  tmp = czero
  do i = 1, n
     tmp = tmp + conjg(x(i)) * y(i)
  end do

  util_zdotc = tmp
  return

end function util_zdotc
!######################################################################
!######################################################################
complex(kind(0d0)) function util_zdotcp(n, a, b)

  use omp_mod
  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: a(1:n)
  complex(kind(0d0)), intent(in) :: b(1:n)

  integer :: i
  complex(kind(0d0)) :: tmp

  tmp = czero
!$omp parallel default(shared) reduction(+:tmp)
!$omp do
  do i = 1, n
     tmp = tmp + conjg(a(i)) * b(i)
  end do
!$omp end do
!$omp end parallel

  util_zdotcp = tmp
  return

end function util_zdotcp
!######################################################################
!######################################################################
complex(kind(0d0)) function util_zdotu(n, x, incx, y, incy)

  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n, incx, incy
  complex(kind(0d0)), intent(in) :: x(1:n)
  complex(kind(0d0)), intent(in) :: y(1:n)

  integer :: i
  complex(kind(0d0)) :: tmp

  if (incx /= 1) stop 'util_zdotu: bad incx.'
  if (incy /= 1) stop 'util_zdotu: bad incy.'

  tmp = czero
  do i = 1, n
     tmp = tmp + x(i) * y(i)
  end do

  util_zdotu = tmp
  return

end function util_zdotu
!######################################################################
!######################################################################
complex(kind(0d0)) function util_zdotup(n, a, b)

  use omp_mod
  use const_mod, only : czero

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: a(1:n)
  complex(kind(0d0)), intent(in) :: b(1:n)

  integer :: i
  complex(kind(0d0)) :: tmp

  tmp = czero
!$omp parallel default(shared) reduction(+:tmp)
!$omp do
  do i = 1, n
     tmp = tmp + a(i) * b(i)
  end do
!$omp end do
!$omp end parallel

  util_zdotup = tmp
  return

end function util_zdotup
!######################################################################
