!################################################################################
subroutine util_diag_ovlp_real(n, s, a, u)

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: s(1:n, 1:n)
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: u(1:n, 1:n)

  integer :: info, lwork, i, len
  real(kind(0d0)) :: rlwork

  real(kind(0d0)), allocatable :: tmps(:,:), tmpa(:, :)
  real(kind(0d0)), allocatable :: work(:)
  real(kind(0d0)), allocatable :: eig(:)

  len = max(n*n, 3*n-2)
  allocate(tmps(n, n))
  allocate(tmpa(n, n))
  allocate(eig(n))

  tmps(1:n, 1:n) = real(s(1:n, 1:n))
  tmpa(1:n, 1:n) = real(a(1:n, 1:n))

  call dsyev(1, "v", "u", n, tmpa, n, tmps, n, eig, rlwork, -1, info)
  lwork = int(rlwork)

  allocate(work(lwork))
  call dsyev(1, "v", "u", n, tmpa, n, tmps, n, eig, work, lwork, info)

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag.'
  else
     a(1:n, 1:n) = (0.d+0, 0.d+0)
     u(1:n, 1:n) = (0.d+0, 0.d+0)
     do i = 1, n
        a(i, i) = eig(n - i + 1)
        u(1:n, i) = tmpa(1:n, n - i + 1)
     end do
  end if
  deallocate(work)

  deallocate(eig)
  deallocate(tmpa)
  deallocate(tmps)

end subroutine util_diag_ovlp_real
!################################################################################
