program x2ewfn2

  implicit none
  integer, parameter :: ior = 5
  integer, parameter :: iow = 6
  character(len = 256) :: fname
  integer :: ngrid, igrid, jgrid
  real(kind(0d0)) :: r2den
  real(kind(0d0)), allocatable :: x(:)
  complex(kind(0d0)), allocatable :: wfn(:,:)

  read(ior, "( I25)") ngrid

  allocate(x(0:ngrid))
  allocate(wfn(0:ngrid, 0:ngrid))

  read(ior, "( E25.15)") x(0:ngrid)
  read(ior, "(2E25.15)") wfn(0:ngrid, 0:ngrid)

  write(iow, "('# x, y, real, imag, r2den')")
  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        r2den = real(wfn(igrid, jgrid) * conjg(wfn(igrid, jgrid)))
        write(iow, "(5F20.10)") x(igrid), x(jgrid), wfn(igrid, jgrid), r2den
     end do
     write(iow, *)
  end do


  deallocate(wfn)
  deallocate(x)

end program x2ewfn2
