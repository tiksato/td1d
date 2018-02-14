!################################################################################
subroutine util_diag_ovlp_comp(n, s, a, u)

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: s(1:n, 1:n)
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: u(1:n, 1:n)

  integer :: info, lwork, i, len
  complex(kind(0d0)) :: clwork

  complex(kind(0d0)), allocatable :: tmps(:,:), tmpa(:, :)
  complex(kind(0d0)), allocatable :: work(:)
  real(kind(0d0)), allocatable :: eig(:), rwork(:)

  len = max(n*n, 3*n-2)
  allocate(rwork(len))
  allocate(tmps(n, n))
  allocate(tmpa(n, n))
  allocate(eig(n))

  tmps(1:n, 1:n) = s(1:n, 1:n)
  tmpa(1:n, 1:n) = a(1:n, 1:n)

  call zhegv(1, "v", "u", n, tmpa, n, tmps, n, eig, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zhegv(1, "v", "u", n, tmpa, n, tmps, n, eig, work, lwork, rwork, info)

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag_comp.'
  else
     u(1:n, 1:n) = (0.d+0, 0.d+0)
     a(1:n, 1:n) = (0.d+0, 0.d+0)
     do i = 1, n
        a(i, i) = eig(n - i + 1)
        u(1:n, i) = tmpa(1:n, n - i + 1)
     end do
  end if
  deallocate(work)

  deallocate(eig)
  deallocate(tmpa)
  deallocate(tmps)
  deallocate(rwork)

end subroutine util_diag_ovlp_comp
!################################################################################
