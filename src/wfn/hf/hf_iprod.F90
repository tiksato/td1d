!######################################################################
complex(kind(0d0)) function hf_iprod(rmax, wfnl, wfnr)

  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin
  use mol_mod, only : ne

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid, 1:nfun, 1:nspin)

  complex(kind(0d0)), external :: util_det
  complex(kind(0d0)), allocatable :: ovlp(:,:,:)
  integer :: ispin
  complex(kind(0d0)) :: dets(2)

  allocate(ovlp(1:nfun, 1:nfun, 1:nspin))

  call hf_ovlp(.true., rmax, wfnl, wfnr, ovlp)
!2e only  hf_iprod = ovlp(1, 1, 1) * ovlp(1, 1, nspin)

  do ispin = 1, nspin
     dets(ispin) = util_det(ne(ispin), thrwfn, ovlp(1:ne(ispin), 1:ne(ispin), ispin))
  end do
  hf_iprod = dets(1) * dets(nspin)
    
  deallocate(ovlp)
  return
    
end function hf_iprod
!######################################################################
