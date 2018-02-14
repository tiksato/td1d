!################################################################################
subroutine grid_fdcoeff2

  implicit none

  call grid_fdcoeff2_center
  call grid_fdcoeff2_forward
  call grid_fdcoeff2_backward

end subroutine grid_fdcoeff2
!################################################################################
!################################################################################
subroutine grid_fdcoeff2_center
!
! centered formula for inner grid points
!
  use const_mod, only : zero, one, two
  use grid_mod, only : dgrid, fd_order, fd_ohalf, fd_coeff_in

  implicit none
  integer :: npt, k, n
  real(kind(0d0)) :: fdfac
  real(kind(0d0)), allocatable :: amat(:,:), bvec(:)
  integer, external :: util_ifact

! factor
  fdfac = - one / dgrid ** two
!debug  fdfac = two
!
  npt = fd_order                            ! NOTE: number of terms is npt + 1 = fd_order + 1
  allocate(amat(0:npt, -fd_ohalf:fd_ohalf))
  allocate(bvec(0:npt))

  bvec(0:npt) = zero
  bvec(2) = one

  do n = 0, npt
     do k = -fd_ohalf, fd_ohalf
        if (k /= 0) then
           amat(n, k) = real(k ** n)
        else
           if (n == 0) then
              amat(n, k) = one
           else
              amat(n, k) = zero
           end if
        end if
     end do
  end do

!debug
!write(6,*) 'grid_fdcoeff2: centered formula with order ', 2 * fd_ohalf
!write(6,*) 'amat:'
!write(6,"(3f10.5)") amat(0:npt, -fd_ohalf:fd_ohalf)
!write(6,*) 'bvec:'
!write(6,"(3f10.5)") bvec(0:npt)
!debug

  call util_syseq_real(npt + 1, amat, bvec)
  call dcopy(npt + 1, bvec(0), 1, fd_coeff_in(-fd_ohalf), 1)
  call dscal(npt + 1, fdfac, fd_coeff_in(-fd_ohalf), 1)

!debug
!write(6,*) 'centered:'
!write(6,"(f10.5)") fd_coeff_in(-fd_ohalf:fd_ohalf)
!debug
  deallocate(bvec)
  deallocate(amat)

end subroutine grid_fdcoeff2_center
!################################################################################
!################################################################################
subroutine grid_fdcoeff2_forward
!
! Forward formulae for grid points around the left boundary
!
  use const_mod, only : zero, one, two
  use grid_mod, only : dgrid, fd_order, fd_ohalf, fd_coeff_l

  implicit none
  integer :: icntr, npt, keff, k, n
  real(kind(0d0)) :: fdfac
  real(kind(0d0)), allocatable :: amat(:,:), bvec(:)
  integer, external :: util_ifact

! factor
  fdfac = - one / dgrid ** two
!debug  fdfac = two
!
  npt = fd_order + 1            ! NOTE: number of terms is npt + 1 = fd_order + 2
  allocate(amat(0:npt, 0:npt))
  allocate(bvec(0:npt))
!
  do icntr = 0, fd_ohalf - 1    ! NOTE: fd_ohalf points nearest the left boundary 0 -> 1 -> 2 -> ...
!
     bvec(0:npt) = zero
     bvec(2) = one

     do n = 0, npt
        do k = 0, npt
           keff = k - icntr

           if (keff /= 0) then
              amat(n, k) = real(keff ** n)
           else
              if (n == 0) then
                 amat(n, k) = one
              else
                 amat(n, k) = zero
              end if
           end if
        end do
     end do

     call util_syseq_real(npt + 1, amat, bvec)
     call dcopy(npt + 1, bvec(0), 1, fd_coeff_l(0, icntr), 1)
     call dscal(npt + 1, fdfac, fd_coeff_l(0, icntr), 1)

!debug
!write(6,*) 'forward: icntr = ', icntr
!write(6,"(f10.5)") fd_coeff_l(0:npt, icntr)
!debug
!
  end do

  deallocate(bvec)
  deallocate(amat)

end subroutine grid_fdcoeff2_forward
!################################################################################
!################################################################################
subroutine grid_fdcoeff2_backward
!
! Backward formulae for grid points around the right boundary
!
  use const_mod, only : zero, one, two
  use grid_mod, only : dgrid, fd_order, fd_ohalf, fd_coeff_r

  implicit none
  integer :: icntr, npt, keff, k, n
  real(kind(0d0)) :: fdfac
  real(kind(0d0)), allocatable :: amat(:,:), bvec(:)
  integer, external :: util_ifact

! factor
  fdfac = - one / dgrid ** two
!debug  fdfac = two
!
  npt = fd_order + 1            ! NOTE: number of terms is npt + 1 = fd_order + 2
  allocate(amat(0:npt, 0:npt))
  allocate(bvec(0:npt))
!
  do icntr = 0, fd_ohalf - 1    ! NOTE: fd_ohalf points nearest the right boundary ... <- 2 <- 1 <- 0
!
     bvec(0:npt) = zero
     bvec(2) = one

     do n = 0, npt
        do k = 0, npt
           keff = k - npt + icntr

           if (keff /= 0) then
              amat(n, k) = real(keff ** n)
           else
              if (n == 0) then
                 amat(n, k) = one
              else
                 amat(n, k) = zero
              end if
           end if
        end do
     end do

     call util_syseq_real(npt + 1, amat, bvec)
     call dcopy(npt + 1, bvec(0), 1, fd_coeff_r(0, icntr), 1)
     call dscal(npt + 1, fdfac,      fd_coeff_r(0, icntr), 1)

!debug
!write(6,*) 'backward: icntr = ', icntr
!write(6,"(f10.5)") fd_coeff_r(0:npt, icntr)
!debug
!
  end do

  deallocate(bvec)
  deallocate(amat)

end subroutine grid_fdcoeff2_backward
!################################################################################
