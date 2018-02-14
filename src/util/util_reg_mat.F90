!################################################################################
subroutine util_reg_mat(n, type, thresh, a, areg)

  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n, type
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: areg(1:n, 1:n)

  if (type == 1) then
     call util_reg_mat1(n, thresh, a, areg)
  else if (type == 2) then
     call util_reg_mat2(n, thresh, a, areg)
  else
     stop 'bad type in util_reg_mat.'
  end if
end subroutine util_reg_mat
!################################################################################
!################################################################################
subroutine util_reg_mat1(n, thresh, a, areg)

  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: areg(1:n, 1:n)

  integer :: kfun
  complex(kind(0d0)) :: diag, expd

  call util_zcopy(n*n, a, 1, areg, 1)

  do kfun = 1, n
     diag = abs(areg(kfun, kfun)) / thresh
     expd = thresh * exp(-diag)
     areg(kfun, kfun) = areg(kfun, kfun) + expd
  end do

end subroutine util_reg_mat1
!################################################################################
!################################################################################
subroutine util_reg_mat2(n, thresh, a, areg)

  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: areg(1:n, 1:n)

  integer :: ifun, jfun, kfun
  complex(kind(0d0)) :: diag, expd
  complex(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)), allocatable :: work(:,:)

  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  call util_zcopy(n*n, a, 1, work, 1)
  call util_zcopy(n*n, a, 1, areg, 1)
  call util_diag_comp(.true., n, work, uvec)
!debug
!  write(6, "('diagd:')")
!  do kfun = 1, n
!     write(6, "(i5, 2f20.10)") kfun, work(kfun, kfun)
!  end do
!  write(6, "('diagd/thresh:')")
!  do kfun = 1, n
!     write(6, "(i5, 2f20.10)") kfun, work(kfun, kfun) / thresh
!  end do
!  write(6, "('expd:')")
!  do kfun = 1, n
!     write(6, "(i5, 2f20.10)") kfun, thresh * exp(-abs(work(kfun, kfun)) / thresh)
!  end do
!  stop 'util_reg_mat (1)'
!debug
  do kfun = 1, n
     diag = abs(work(kfun, kfun)) / thresh
     expd = thresh * exp(-diag)
     do jfun = 1, n
        do ifun = 1, n
           areg(ifun, jfun) = areg(ifun, jfun) + uvec(ifun, kfun) * expd &
                                       & * conjg(uvec(jfun, kfun))
        end do
     end do
  end do

  deallocate(work)
  deallocate(uvec)

end subroutine util_reg_mat2
!################################################################################
