!######################################################################
real(kind(0d0)) function hf_n2err(wfn1, wfn2)

  use const_mod, only : one, czero, ctwo
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, ncore, nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfn1(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: wfn2(0:ngrid, 1:nfun, 1:nspin)

  complex(kind(0d0)), allocatable :: ovlp(:,:,:)
  integer :: ifun, ispin
  complex(kind(0d0)) :: tmp, hffac

  allocate(ovlp(1:nfun, 1:nfun, 1:nspin))

  call hf_ovlp(.true., -one, wfn2, wfn2, ovlp)

  hffac = ctwo / nspin
  tmp = czero
  do ispin = 1, nspin
     do ifun = 1, ncore
        tmp = tmp + ovlp(ifun, ifun, ispin) * hffac
     end do
  end do

  hf_n2err = abs(tmp)

  deallocate(ovlp)

!old  hf_n2err = abs(hf_iprod(-one, wfn2, wfn2))
  
end function hf_n2err
!######################################################################
