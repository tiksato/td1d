!################################################################################
subroutine util_hlineq_reg(n, thresh, a, yvec, xvec)
!
! A * X = Y
! a * X' = Y', X' = U+ * X, Y' = U+ * Y
! X' = Y' / a'
! X = U * X'
!
  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(in) :: yvec(1:n)
  complex(kind(0d0)), intent(out) :: xvec(1:n)

  integer :: ifun, jfun, n_act, icol
  integer, parameter :: max_col = 10
  complex(kind(0d0)) :: diag, invd, inva
  complex(kind(0d0)), allocatable :: ytmp(:)
  complex(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)), allocatable :: work(:,:)

  allocate(ytmp(1:n))
  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  work(1:n, 1:n) = a(1:n, 1:n)
  uvec(1:n, 1:n) = 0.d+0
  call util_diag_comp(.true., n, work, uvec)
  
!debug
  write(6, "('Eigenvalues of a: ')")
  ifun = 1
  do
     do icol = 1, max_col
        if (ifun > n) exit
        write(6, "(f15.8)", advance = 'no') dble(work(ifun, ifun))
        ifun = ifun + 1
     end do
     write(6, *)
     if (ifun >= n) exit
  end do
!debug

!####################
  n_act = n
!####################

  ytmp(1:n) = 0.d+0
  do ifun = 1, n_act
     do jfun = 1, n
        ytmp(ifun) = ytmp(ifun) + conjg(uvec(jfun, ifun)) * yvec(jfun)
     end do
  end do

  if (thresh > 0.d+0) then
     do ifun = 1, n_act
        diag = work(ifun, ifun)
        inva = diag / (diag * diag + thresh)
        ytmp(ifun) = ytmp(ifun) * inva
     end do
  else
     do ifun = 1, n_act
        diag = work(ifun, ifun)
        inva = runit / diag
        ytmp(ifun) = ytmp(ifun) * inva
     end do
  end if

  xvec(1:n) = czero
  do ifun = 1, n
     do jfun = 1, n_act
        xvec(ifun) = xvec(ifun) + uvec(ifun, jfun) * ytmp(jfun)
     end do
  end do

  deallocate(work)
  deallocate(uvec)
  deallocate(ytmp)

end subroutine util_hlineq_reg
!################################################################################
