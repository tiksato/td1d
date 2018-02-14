!######################################################################
subroutine hf_print(iunit, wfn)
  use grid_mod, only : ngrid, x
  use wfn_mod, only : nfun, nspin

  integer, intent(in) :: iunit
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  integer :: igrid, ifun, ispin

  write(iunit, "('mo coefficients:', i10)") ngrid
  do igrid = 0, ngrid
     write(iunit, "(F20.10)", advance = 'no') x(igrid)
     do ispin = 1, nspin
        do ifun = 1, nfun
           write(iunit, "(2F20.10)", advance = 'no') wfn(igrid, ifun, ispin)
        end do
     end do
     write(iunit, *)
  end do
end subroutine hf_print
!######################################################################
