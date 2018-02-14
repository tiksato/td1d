!######################################################################
subroutine hf_full_diag(fock, wfn)

  use io_mod, only : iostdo
  use root_mod, only : iprint
  use grid_mod, only : ngrid
  use wfn_mod, only : nspin, nfun

  implicit none
  complex(kind(0d0)), intent(inout) :: fock(0:ngrid, 0:ngrid, 1:nspin)
  complex(kind(0d0)), intent(out) ::    wfn(0:ngrid, 0:ngrid, 1:nspin)
  integer :: ifun, ispin

  do ispin = 1, nspin
     call util_diag_comp(.false., ngrid+1, fock(0,0,ispin), wfn(0,0,ispin))
  end do

  if (iprint > 4) then
     write(iostdo, "('# hf_full_diag: eigenvalues')")
     do ifun = 0, nfun - 1
        write(iostdo, "(i10)", advance="no") ifun
        write(iostdo, "(f20.10)", advance="no") dble(fock(ifun, ifun, 1))
        if (nspin == 2) write(iostdo, "(f20.10)", advance="no") dble(fock(ifun, ifun, 2))
        write(iostdo, *)
     end do
  end if

end subroutine hf_full_diag
!######################################################################
