!######################################################################
subroutine hf_ovdef(s)

  use const_mod, only : czero, runit
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(out) :: s(1:nfun, 1:nfun, 1:nspin)

  integer :: ifun, ispin

  s(1:nfun, 1:nfun, 1:nspin) = czero

  do ispin = 1, nspin
     do ifun = 1, nfun
        s(ifun, ifun, ispin) = runit
     end do
  end do

end subroutine hf_ovdef
!######################################################################
