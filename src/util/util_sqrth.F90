!################################################################################
subroutine util_sqrth(doinv, thresh, n, hmat, sqrth)

  implicit none
  logical, intent(in) :: doinv
  real(kind(0d0)), intent(in) :: thresh
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: hmat(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: sqrth(1:n, 1:n, 1:*)

  integer :: i, j, k
  real(kind(0d0)) :: diag, do12
  complex(kind(0d0)), allocatable :: htmp(:,:)
  complex(kind(0d0)), allocatable :: uvec(:,:)
  real(kind(0d0)), parameter :: zero = 0.d+0
  real(kind(0d0)), parameter :: one = 1.d+0
  complex(kind(0d0)), parameter :: czero = (0.d+0, 0.d+0)

  allocate(htmp(1:n, 1:n))
  allocate(uvec(1:n, 1:n))

  htmp(1:n, 1:n) = hmat(1:n, 1:n)
  uvec(1:n, 1:n) = czero
  call util_diag_comp(.false., n, htmp, uvec)

  sqrth(1:n, 1:n, 1) = czero
  if (doinv) sqrth(1:n, 1:n, 2) = czero

  do k = 1, n
     diag = dble(htmp(k, k))
     if (diag < zero) stop 'Non-positive-definite H in util_sqrth.'

     do12 = sqrt(diag)
     do i = 1, n
        do j = 1, n
           sqrth(j, i, 1) = sqrth(j, i, 1) + uvec(j, k) * do12 * conjg(uvec(i, k))
        end do
     end do

     if (doinv .and. diag > thresh) then
        do12 = one / sqrt(diag)
        do i = 1, n
           do j = 1, n
              sqrth(j, i, 2) = sqrth(j, i, 2) + uvec(j, k) * do12 * conjg(uvec(i, k))
           end do
        end do
     end if
  end do

  deallocate(uvec)
  deallocate(htmp)

end subroutine util_sqrth
!################################################################################
