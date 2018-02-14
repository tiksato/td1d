!################################################################################
subroutine util_svd_comp(m, n, a, z, u, v)

  implicit none
  integer, intent(in) :: m, n
  complex(kind(0d0)), intent(in) :: a(1:m, 1:n)
  real(kind(0d0)), intent(out) :: z(1:*)
  complex(kind(0d0)), intent(out) :: u(1:m, 1:m)
  complex(kind(0d0)), intent(out) :: v(1:n, 1:n)

  integer :: info, i, j, lcwork, lrwork
  real(kind(0d0)), allocatable :: rwork(:)
  complex(kind(0d0)), allocatable :: cwork(:)
  complex(kind(0d0)), allocatable :: tmpa(:,:)
  complex(kind(0d0)), allocatable :: vt(:,:)

  lrwork = 5 * min(m, n)
  lcwork = 2 * min(m, n) + max(m, n)

  allocate(rwork(1:lrwork))
  allocate(cwork(1:lcwork))
  allocate(tmpa(1:m, 1:n))
  allocate(vt(1:n, 1:n))

  call zgesvd('A', 'A', m, n, tmpa, m, z, u, m, vt, n, cwork, -1, rwork, info)
  lcwork = int(cwork(1))
  deallocate(cwork)
  allocate(cwork(1:lcwork))

  tmpa(1:m, 1:n) = a(1:m, 1:n)
  call zgesvd('A', 'A', m, n, tmpa, m, z, u, m, vt, n, cwork, lcwork, rwork, info)

  if (info /= 0) then
     write(6, "('util_svd_comp: info = ', i20)") info
     stop 'error in util_svd_comp.'
  else
     do i = 1, n
        do j = 1, n
           v(j, i) = conjg(vt(i, j))
        end do
     end do
  end if

!  write(6, "('util_svd_comp. input matrix-R:')")
!  do i = 1, m
!     do j = 1, n
!        write(6, "(f20.10)", advance = 'no') dble(a(i, j))
!     end do
!     write(6, *)
!  end do
!  write(6, "('util_svd_comp. input matrix-I:')")
!  do i = 1, m
!     do j = 1, n
!        write(6, "(f20.10)", advance = 'no') aimag(a(i, j))
!     end do
!     write(6, *)
!  end do
!  write(6, "('util_svd_comp. singular values:')")
!  do i = 1, min(m, n)
!     write(6, "(i10, f20.10)") i, z(i)
!  end do

  deallocate(vt)
  deallocate(tmpa)
  deallocate(cwork)
  deallocate(rwork)

end subroutine util_svd_comp
!################################################################################
