!######################################################################
subroutine hf_write_orb(iow, wfn)

  use grid_mod, only : ngrid, x
  use wfn_mod, only : nfun, nspin, nbiort

  implicit none
  integer, intent(in) :: iow
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)
  integer :: igrid, ifun, is

  do is = 1, nspin
     write(iow, "('# is =', i5)") is
     do igrid = 0, ngrid
        write(iow, "(F20.10)", advance = 'no') x(igrid)
        do ifun = 1, nfun
           write(iow, "(2F25.15)", advance = 'no') wfn(igrid,ifun,is)
        end do
        write(iow, *)
     end do
     write(iow, *)
     write(iow, *)
  end do

end subroutine hf_write_orb
!######################################################################
