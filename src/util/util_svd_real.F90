!################################################################################
subroutine util_svd_real(m, n, a, z, u, v)

  implicit none
  integer, intent(in) :: m, n
  real(kind(0d0)), intent(in) :: a(1:m, 1:n)
  real(kind(0d0)), intent(out) :: z(1:*)
  real(kind(0d0)), intent(out) :: u(1:m, 1:m)
  real(kind(0d0)), intent(out) :: v(1:n, 1:n)

  integer :: info, i, j, len
  integer, allocatable :: iwork(:)
  real(kind(0d0)), allocatable :: rwork(:)
  real(kind(0d0)), allocatable :: tmpa(:,:)
  real(kind(0d0)), allocatable :: vt(:,:)

  len = 8*min(m, n)
!debug
!  write(6, "('util_svd_real: len = ', i10)") len
!debug
  allocate(iwork(1:len))
  allocate(rwork(1:1))
  allocate(tmpa(1:m, 1:n))
  allocate(vt(1:n, 1:n))
  call dgesdd('A', m, n, tmpa, m, z, u, m, vt, n, rwork, -1, iwork, info)
  len = int(rwork(1))
!debug
!  write(6, "('util_svd_real: len = ', i10)") len
!debug
  deallocate(rwork)
  allocate(rwork(1:len))


  vt = 0d0
  iwork = 0
  rwork = 0d0
  tmpa(1:m, 1:n) = a(1:m, 1:n)
  call dgesdd('A', m, n, tmpa, m, z, u, m, vt, n, rwork, len, iwork, info)

  if (info /= 0) then
     write(6, "('util_svd_real: info = ', i20)") info
     stop 'error in util_svd_real.'
  else
     do i = 1, n
        do j = 1, n
           v(j, i) = vt(i, j)
        end do
     end do
  end if

!DEBUG
!  write(6, "('util_svd_real. singular values:')")
!  do i = 1, min(m, n)
!     write(6, "(i10, f20.10)") i, z(i)
!  end do
!DEBUG

  deallocate(vt)
  deallocate(tmpa)
  deallocate(rwork)
  deallocate(iwork)

end subroutine util_svd_real
!################################################################################
