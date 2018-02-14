!################################################################################
subroutine util_diag_real(dsc, n, a, u)

  implicit none
  logical, intent(in) :: dsc
  integer, intent(in) :: n
  real(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  real(kind(0d0)), intent(out) :: u(1:n, 1:n)

  integer :: info, lwork, i, len
  real(kind(0d0)) :: rlwork

  real(kind(0d0)), allocatable :: tmpa(:, :)
  real(kind(0d0)), allocatable :: work(:)
  real(kind(0d0)), allocatable :: eig(:)

  len = max(n*n, 3*n-2)
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

  call dsyev("v", "u", n, tmpa, n, eig, rlwork, -1, info)
  lwork = int(rlwork)

  allocate(work(lwork))
  call dsyev("v", "u", n, tmpa, n, eig, work, lwork, info)

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag.'
  else
     a(1:n, 1:n) = 0.d+0
     u(1:n, 1:n) = 0.d+0
     if (dsc) then
        do i = 1, n
           a(i, i) = eig(n - i + 1)
           u(1:n, i) = tmpa(1:n, n - i + 1)
        end do
     else
        do i = 1, n
           a(i, i) = eig(i)
           u(1:n, i) = tmpa(1:n, i)
        end do
     end if
  end if
  deallocate(work)

  deallocate(eig)
  deallocate(tmpa)

end subroutine util_diag_real
!################################################################################
!################################################################################
subroutine util_diag_real2(dsc, n, s, a, u)

  implicit none
  logical, intent(in) :: dsc
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: s(1:n, 1:n)
  real(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  real(kind(0d0)), intent(out) :: u(1:n, 1:n)

  integer :: info, lwork, i, len
  real(kind(0d0)) :: rlwork

  real(kind(0d0)), allocatable :: tmpa(:, :)
  real(kind(0d0)), allocatable :: work(:)
  real(kind(0d0)), allocatable :: eig(:)

  len = max(n*n, 3*n-2)
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

  call dsygv(1, "v", "u", n, tmpa, n, s, n, eig, rlwork, -1, info)
  lwork = int(rlwork)

  allocate(work(lwork))
  call dsygv(1, "v", "u", n, tmpa, n, s, n, eig, work, lwork, info)

  if (info /= 0) then
     write(6, "('util_diag_real2: info = ', i20)") info
     stop 'error in util_diag_real2'
  else
     u(1:n, 1:n) = 0.d+0
     a(1:n, 1:n) = 0.d+0
     if (dsc) then
        do i = 1, n
           a(i, i) = eig(n - i + 1)
           u(1:n, i) = tmpa(1:n, n - i + 1)
        end do
     else
        do i = 1, n
           a(i, i) = eig(i)
           u(1:n, i) = tmpa(1:n, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:n)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(tmpa)

end subroutine util_diag_real2
!################################################################################
