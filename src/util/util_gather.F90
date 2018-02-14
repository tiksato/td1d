!######################################################################
subroutine util_gather(nvar, ndata, data)
!
! Accumulate processor-distributed data
!
  use const_mod, only : runit

  implicit none
  integer, intent(in) :: nvar, ndata
  complex(kind(0d0)), intent(inout) :: data(1:nvar, 1:ndata)

  integer :: idata

  do idata = 2, ndata
     call util_zaxpy(nvar, runit, data(1,idata), 1, data(1,1), 1)
  end do

end subroutine util_gather
!######################################################################
