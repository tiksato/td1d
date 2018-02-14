!################################################################################
subroutine util_matinv_reg(n, thresh, a, inva)

  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: inva(1:n, 1:n)

  integer :: ifun, jfun, kfun
  complex(kind(0d0)) :: diag, invd
  complex(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)), allocatable :: work(:,:)

  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  call util_zcopy(n*n, a, 1, work, 1)
  call util_zcopy(n*n, czero, 0, inva, 1)
  call util_zcopy(n*n, czero, 0, uvec, 1)
  call util_diag_comp(.true., n, work, uvec)
  
  do kfun = 1, n
     diag = work(kfun, kfun)
!DEBUG
!     write(6, "('util_matinv_reg: ', i5, 2f20.10)") kfun, diag
!DEBUG
!    invd = runit / diag
     invd = diag / (diag * diag + thresh * thresh)
!    invd = runit / (diag + thresh * exp(-diag/thresh))
     do jfun = 1, n
        do ifun = 1, n
           inva(ifun, jfun) = inva(ifun, jfun) + uvec(ifun, kfun) * invd &
                                       & * conjg(uvec(jfun, kfun))
        end do
     end do
  end do

  deallocate(work)
  deallocate(uvec)

end subroutine util_matinv_reg
!################################################################################
!################################################################################
subroutine util_matinv_reg_real(n, thresh, a, inva)

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  real(kind(0d0)), intent(in) :: a(1:n, 1:n)
  real(kind(0d0)), intent(out) :: inva(1:n, 1:n)

  integer :: ifun, jfun, kfun
  real(kind(0d0)) :: diag, invd
  real(kind(0d0)), allocatable :: uvec(:,:)
  real(kind(0d0)), allocatable :: work(:,:)

  allocate(uvec(1:n, 1:n))
  allocate(work(1:n, 1:n))

  work(1:n, 1:n) = a(1:n, 1:n)
  inva(1:n, 1:n) = 0.d+0
  uvec(1:n, 1:n) = 0.d+0
  call util_diag_real(.true., n, work, uvec)
  
  do kfun = 1, n
     diag = work(kfun, kfun)
     invd = diag / (diag * diag + thresh * thresh)
     do jfun = 1, n
        do ifun = 1, n
           inva(ifun, jfun) = inva(ifun, jfun) + uvec(ifun, kfun) * invd &
                                             & * uvec(jfun, kfun)
        end do
     end do
  end do

  deallocate(work)
  deallocate(uvec)

end subroutine util_matinv_reg_real
!################################################################################
