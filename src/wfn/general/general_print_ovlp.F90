!######################################################################
subroutine general_print_ovlp(io, rmax, wfn)

  use wfn_mod, only : nfun

  implicit none
  integer, intent(in) :: io
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: ifun, jfun
  complex(kind(0d0)), allocatable :: ovlp(:,:)

  allocate(ovlp(1:nfun, 1:nfun))

  call general_ovlp(.true., rmax, wfn, wfn, ovlp)
  do ifun = 1, nfun
     do jfun = 1, ifun
        write(io, "(2E20.10)", advance='no') ovlp(jfun, ifun)
     end do
  end do
  write(io, *)

  deallocate(ovlp)

end subroutine general_print_ovlp
!######################################################################
