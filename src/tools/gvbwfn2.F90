program gvbwfn2

  implicit none
  integer, parameter :: ior = 5
  integer, parameter :: iow = 6
  character(len = 256) :: fname
  character(len = 256) :: carg
  integer :: ngrid, igrid, jgrid, nfun
  real(kind(0d0)) :: r2den
  real(kind(0d0)), allocatable :: x(:)
  complex(kind(0d0)), allocatable :: wfn(:,:)
  complex(kind(0d0)), allocatable :: wfn2e(:,:)

  read(ior, "( I25)") ngrid

! two orbital gvb
  nfun = 2

  allocate(x(0:ngrid))
  allocate(wfn(0:ngrid, 1:nfun))
  allocate(wfn2e(0:ngrid, 0:ngrid))

  read(ior, "( E25.15)") x(0:ngrid)
  read(ior, "(2E25.15)") wfn(0:ngrid, 1:nfun)

  allocate(ovdef(1:nfun, 1:nfun))
  call gvb_ovdef(wfn, ovdef)
  norm0 = two * real(ovdef(1, 1) * ovdef(2, 2) + ovdef(1, 2) * ovdef(2, 1))
  deallocate(ovdef)


  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        wfn2e(jgrid, igrid) = wfn(jgrid, 1) * wfn(igrid, nspin)
     end do
  end do

  write(iow, "('# x, y, real, imag, r2den')")
  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        r2den = real(wfn2e(igrid, jgrid) * conjg(wfn2e(igrid, jgrid)))
        write(iow, "(5F20.10)") x(igrid), x(jgrid), wfn2e(igrid, jgrid), r2den
     end do
     write(iow, *)
  end do


  deallocate(wfn2e)
  deallocate(wfn)
  deallocate(x)

end program gvbwfn2
