!################################################################################
subroutine util_expihd(n, d, hmat, exph)

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: d
  complex(kind(0d0)), intent(in) :: hmat(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: exph(1:n, 1:n)
  integer :: i, j, k
  real(kind(0d0)) :: argd
  complex(kind(0d0)) :: expd
  complex(kind(0d0)), allocatable :: htmp(:,:)
  complex(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)), parameter :: czero = (0.d+0, 0.d+0)
  complex(kind(0d0)), parameter :: runit = (1.d+0, 0.d+0)
  complex(kind(0d0)), parameter :: iunit = (0.d+0, 1.d+0)

  allocate(htmp(1:n, 1:n))
  allocate(uvec(1:n, 1:n))

  htmp(1:n, 1:n) = hmat(1:n, 1:n)
  uvec(1:n, 1:n) = czero
  call util_diag_comp(.false., n, htmp, uvec)

!debug  write(6, "('eigenvalues of hmat:')")
!debug  do j = 1, n
!debug     write(6, "(f12.5)", advance = 'no') dble(htmp(j, j))
!debug  end do
!debug  write(6, *)

  exph(1:n, 1:n) = czero
  do k = 1, n
     argd = dble(htmp(k, k)) * d
     expd = runit * cos(argd) + iunit * sin(argd)
     do i = 1, n
        do j = 1, n
           exph(j, i) = exph(j, i) + uvec(j, k) * expd * conjg(uvec(i, k))
        end do
     end do
  end do

  deallocate(uvec)
  deallocate(htmp)

end subroutine util_expihd
!################################################################################
