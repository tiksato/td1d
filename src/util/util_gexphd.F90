!################################################################################
subroutine util_gexphd(n, fac, hmat, exph)

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: fac
  complex(kind(0d0)), intent(in) :: hmat(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: exph(1:n, 1:n)
  integer :: i, j, k
  complex(kind(0d0)) :: expd
  complex(kind(0d0)), allocatable :: htmp(:,:)
  complex(kind(0d0)), allocatable :: lvec(:,:), linv(:,:)
  complex(kind(0d0)), allocatable :: rvec(:,:), rinv(:,:)
  real(kind(0d0)), parameter :: tiny = 1.D-20
  complex(kind(0d0)), parameter :: czero = (0.d+0, 0.d+0)
  complex(kind(0d0)), parameter :: runit = (1.d+0, 0.d+0)
  complex(kind(0d0)), parameter :: iunit = (0.d+0, 1.d+0)

  allocate(htmp(1:n, 1:n))
  allocate(lvec(1:n, 1:n)); allocate(linv(1:n, 1:n))
  allocate(rvec(1:n, 1:n)); allocate(rinv(1:n, 1:n))

  htmp(1:n, 1:n) = hmat(1:n, 1:n)
  call util_gdiag_comp_new(.false., n, htmp, lvec, rvec)
  call util_gmatinv(n, tiny, lvec, linv)
  call util_gmatinv(n, tiny, rvec, rinv)

  !DEBUG
  !DEBUG write(6, "('eigenvalues of hmat:')")
  !DEBUG do j = 1, n
  !DEBUG    write(6, "(i12, 2f12.8)") j, htmp(j, j)
  !DEBUG end do
  !DEBUG write(6, *)
  !DEBUG

  exph(1:n, 1:n) = czero
  do k = 1, n
     expd = exp(fac * htmp(k, k))
     !DEBUG
     !DEBUG write(6, "('util_gexphd expd: ', i10, 2f20.10)") k, expd
     !DEBUG
     do i = 1, n
        do j = 1, n
           exph(j, i) = exph(j, i) + conjg(linv(k, j)) * expd * rinv(k, i)
        end do
     end do
  end do

  !DEBUG
  !DEBUG write(6, "('gexphd:')")
  !DEBUG do i = 1, n
  !DEBUG    do j = 1, n
  !DEBUG       write(6, "(2i10, 2f12.8)") i, j, exph(i, j)
  !DEBUG    end do
  !DEBUG end do
  !DEBUG write(6, *)
  !DEBUG

  deallocate(rvec); deallocate(rinv)
  deallocate(lvec); deallocate(linv)
  deallocate(htmp)

end subroutine util_gexphd
!################################################################################
