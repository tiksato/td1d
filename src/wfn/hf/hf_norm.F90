!######################################################################
real(kind(0d0)) function hf_norm(rmax, wfn, p0, p1, p2)

  use const_mod, only : zero, one, two
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin
  use mol_mod, only : ne

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)
  real(kind(0d0)), intent(out) :: p0, p1, p2

  complex(kind(0d0)), external :: util_det
  complex(kind(0d0)), allocatable :: ovlp(:,:), ovlpi(:,:,:), ovlpo(:,:,:)
  integer :: ispin, ifun, jfun
  real(kind(0d0)) :: p0s(2), p1s(2), p2s(2)

  allocate(ovlp(1:nfun, 1:nfun))
  allocate(ovlpi(1:nfun, 1:nfun, 1:nspin))
  allocate(ovlpo(1:nfun, 1:nfun, 1:nspin))
  call hf_ovlp(.true., rmax, wfn, wfn, ovlpi)
  ovlpo(1:nfun, 1:nfun, 1:nspin) = one - ovlpi(1:nfun, 1:nfun, 1:nspin)

!2e only
!  p0 = ovlpi(1, 1, 1) * ovlpi(1, 1, nspin)
!
!  p1 = ovlpi(1, 1, 1) * (one - ovlpi(1, 1, nspin)) &
!   & + (one - ovlpi(1, 1, 1)) * ovlpi(1, 1, nspin)
!
!  p2 = (one - ovlpi(1, 1, 1)) * (one - ovlpi(1, 1, nspin))
!
!  hf_norm = p0 + p1 + p2
!2e only

  do ispin = 1, nspin
     ovlp(1:nfun, 1:nfun) = ovlpi(1:nfun, 1:nfun, ispin)
     p0s(ispin) = real(util_det(ne(ispin), thrwfn, ovlp(1:ne(ispin), 1:ne(ispin))))
  
     p1s(ispin) = zero
     do ifun = 1, ne(ispin)
        ovlp(1:nfun, 1:nfun) = ovlpi(1:nfun, 1:nfun, ispin)
        ovlp(1:nfun,   ifun) = ovlpo(1:nfun,   ifun, ispin)
        p1s(ispin) = p1s(ispin) + real(util_det(ne(ispin), thrwfn, ovlp(1:ne(ispin), 1:ne(ispin))))
     end do
  
     p2s(ispin) = zero
     do ifun = 1, ne(ispin)
        do jfun = 1, ifun - 1
           ovlp(1:nfun, 1:nfun) = ovlpi(1:nfun, 1:nfun, ispin)
           ovlp(1:nfun,   ifun) = ovlpo(1:nfun,   ifun, ispin)
           ovlp(1:nfun,   jfun) = ovlpo(1:nfun,   jfun, ispin)
           p2s(ispin) = p2s(ispin) + real(util_det(ne(ispin), thrwfn, ovlp(1:ne(ispin), 1:ne(ispin))))
        end do
     end do
  end do

  if (ne(3) == 1) then
     p0 = p0s(1)
     p1 = one - p0
     p2 = zero
     hf_norm = p0 + p1 + p2
  else
     p0 = p0s(1) * p0s(nspin)
     p1 = p1s(1) * p0s(nspin) + p0s(1) * p1s(nspin)
     p2 = p1s(1) * p1s(nspin) + p2s(1) * p0s(nspin) + p0s(1) * p2s(nspin)
     hf_norm = p0 + p1 + p2
  end if
    
  deallocate(ovlpo)
  deallocate(ovlpi)
  deallocate(ovlp)
  return
    
end function hf_norm
!######################################################################
