!######################################################################
program cas_old2new

  implicit none
  character(len = 256) :: carg, fnamer, fnamew
  integer :: nfun, ngrid, ngrid0
  real(kind(0d0)), allocatable :: x(:)
  complex(kind(0d0)), allocatable :: wfn(:,:)
  complex(kind(0d0)), allocatable :: cic(:,:)

  integer :: igrid, ifun, jfun, ii, jj
  integer, parameter :: ior = 1, iow = 2

  call getarg(1, carg); read(carg,*) nfun
  call getarg(2, carg); read(carg,*) ngrid
  call getarg(3, fnamer)
  call getarg(4, fnamew)

  allocate(x  (0:ngrid))
  allocate(wfn(0:ngrid, 1:nfun))
  allocate(cic(0:nfun, 1:nfun))

  open(unit = ior, file = trim(fnamer), status = 'old', form = 'formatted')
  open(unit = iow, file = trim(fnamew), status = 'new', form = 'formatted')

  read(ior, "('ci coefficients:')")
  do ifun = 1, nfun
     do jfun = 1, ifun
        read(ior, "(2i5,2f20.10)") jj, ii, cic(jfun, ifun)
        cic(ifun, jfun) = cic(jfun, ifun)
     end do
  end do

  read(ior, "('mo coefficients:', i10)") ngrid0
  if (ngrid /= ngrid0) stop 'ngrid .ne. ngrid0.'

  do igrid = 0, ngrid
     read(ior, "(f20.10)", advance = 'no') x(igrid)
     do ifun = 1, nfun
        read(ior, "(2f20.10)", advance = 'no') wfn(igrid, ifun)
     end do
     read(ior, *)
  end do

  write(iow, "(I25)") ngrid
  write(iow, "(E25.15)") x(0:ngrid)
  write(iow, "(2E25.15)") wfn(0:ngrid, 1:nfun)
  write(iow, "(2E25.15)") cic(1:nfun, 1:nfun)  

  close(iow)
  close(ior)

  deallocate(cic)
  deallocate(wfn)
  deallocate(x)

end program cas_old2new
!######################################################################
