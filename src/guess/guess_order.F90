!################################################################################
subroutine guess_order(wfn)

  use io_mod, only : ior
  use root_mod, only : name
  use const_mod, only : czero
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(out) :: wfn(0:ngrid, 1:nfun)
  !--------------------------------------------------------------------
  integer :: ifun, iorb
  character(len = 256) :: fname
  complex(kind(0d0)), allocatable :: twfn(:,:)

  allocate(twfn(0:ngrid, 1:nfun))

  call util_zcopy((ngrid+1)*nfun, wfn, 1, twfn, 1)
  call util_zcopy((ngrid+1)*nfun, czero, 0, wfn, 1)

  fname = trim(name)//".alter"
  open(unit=ior, file=trim(fname), status='old', form='formatted')
  do ifun = 1, nfun
     read(ior, *) iorb
     write(6, "(' # guess_order ', 2i5)") ifun, iorb
     call util_zcopy(ngrid+1, twfn(0, iorb), 1, wfn(0, ifun), 1)
  end do
  close(unit=ior)

  deallocate(twfn)

!  stop 'stop for debug in guess_order.'

end subroutine guess_order
!################################################################################
