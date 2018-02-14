program dft

  use const_mod, only : two, half, czero

  implicit none
  integer :: n, nd, ios, id
  real(kind(0d0)), allocatable :: inp(:)
  complex(kind(0d0)), allocatable :: out(:)
  integer :: plan(8)

  include 'fftw3.f'

  ! read number of data
  nd = 0
  do 
     read(5, *, iostat=ios)
     if (ios /= 0) exit
     nd = nd + 1
  end do

  n = 2 ** int(log(dble(nd)) / log(two) + half)

  ! allocation
  allocate(inp(1:n))
  allocate(out(1:n/2))

  ! read data
  rewind(5)
  read(5, *) inp(1:nd)
  inp(nd+1:n) = czero

!mine  call dft2(nd, inp, out)

  ! make plan
  call dfftw_plan_dft_r2c_1d(plan, n, inp, out, fftw_estimate)
  ! calculation
  call dfftw_execute(plan, inp, out)
  ! destroy plan
  call dfftw_destroy_plan(plan)

  ! print
  do id = 1, nd / 2
     write(6, "(i10, 3e20.10)") id, out(id), (abs(out(id)) / nd) ** two
  end do

  ! deallocation
  deallocate(out)
  deallocate(inp)

end program dft

subroutine dft2(nd, inp, out)

  use const_mod, only : two, pi, czero, iunit

  implicit none
  integer, intent(in) :: nd
  complex(kind(0d0)), intent(in) :: inp(nd)
  complex(kind(0d0)), intent(out) :: out(nd)
  integer :: j, k
  complex(kind(0d0)) :: tmp, root, wj

  root = exp(-two * pi * iunit / nd)

  do j = 1, nd
     wj = root ** j
     tmp = czero
     do k = 1, nd
        tmp = tmp + inp(k) * wj ** k
     end do
     out(j) = tmp
  end do

end subroutine dft2
