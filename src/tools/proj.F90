program proj

  use const_mod, only : zero, two, czero, half

  implicit none
  integer :: ngrid0, ngrid, nval
  real(kind(0d0)), allocatable :: x0(:), x(:)
  real(kind(0d0)), allocatable :: y0(:,:), y(:,:), numer(:)

  character(len = 256) :: fname0, fname, cnval
  integer, parameter :: io0 = 100, io = 200
  real(kind(0d0)), parameter :: thresh = 1.D-6
  integer :: ioerr, igrid, igrid0, ival, ndone, ll, ul, ll0, ul0
  real(kind(0d0)) :: dl, dr, d, diff_l, diff_r, denom

  call getarg(1, cnval)
  call getarg(2, fname0)
  call getarg(3, fname)
  read(cnval, *) nval
  open(unit=io0, file=trim(fname0), status='old', iostat=ioerr)
  if(ioerr /= 0) stop "input file error 0."
  open(unit=io, file=trim(fname), status='old', iostat=ioerr)
  if(ioerr /= 0) stop "input file error."

  allocate(numer(nval))

  ngrid0 = 0
  do
     read(io0, *, iostat=ioerr)
     if(ioerr /= 0) exit
     ngrid0 = ngrid0 + 1
  end do
  allocate(x0(0:ngrid0))
  allocate(y0(0:ngrid0, nval))
  rewind(io0)
  do igrid = 0, ngrid0
     read(io0, *) x0(igrid), y0(igrid, 1:nval)
  end do

  ngrid = 0
  do
     read(io, *, iostat=ioerr)
     if(ioerr /= 0) exit
     ngrid = ngrid + 1
  end do
  allocate(x(0:ngrid))
  allocate(y(0:ngrid, nval))
  rewind(io)
  do igrid = 0, ngrid
     read(io, *) x(igrid)
  end do

  y(0:ngrid, 1:nval) = czero

  do igrid = 0, ngrid
     do igrid0 = 0, ngrid0 - 1
        diff_l = x(igrid) - x0(igrid0)
        diff_r = x(igrid) - x0(igrid0 + 1)
        if (abs(diff_l) < thresh) then
           y(igrid, 1:nval) = y0(igrid0, 1:nval)
           exit
        else if (abs(diff_r) < thresh) then
           y(igrid, 1:nval) = y0(igrid0 + 1, 1:nval)
           exit
        else if (diff_l * diff_r < zero) then
           numer(1:nval) = y0(igrid0, 1:nval) * abs(diff_r) + y0(igrid0 + 1, 1:nval) * abs(diff_l)
           denom = abs(diff_l) + abs(diff_r)
           y(igrid, 1:nval) = numer(1:nval) / denom
           exit
        end if
     end do
  end do

!  if (ngrid /= ndone) stop 'ngrid /= ndone.'

  do igrid = 0, ngrid
     write(6,"(F20.10, F20.10)", advance = "no") x(igrid)
     do ival = 1, nval
        write(6,"(F20.10, F20.10)", advance = "no") y(igrid, ival)
     end do
     write(6, *)
  end do

  deallocate(y)
  deallocate(y0)
  deallocate(x)
  deallocate(x0)
  deallocate(numer)

end program proj

  
