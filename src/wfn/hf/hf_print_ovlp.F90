!######################################################################
subroutine hf_print_ovlp(io, rmax, wfn, norm)

  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin
  use mol_mod, only : ne

  implicit none
  integer, intent(in) :: io
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)
  real(kind(0d0)), intent(out) :: norm

  complex(kind(0d0)), external :: util_det
  complex(kind(0d0)), allocatable :: ovlp(:,:,:)
  integer :: ispin, ifun, jfun
  complex(kind(0d0)) :: dets(2)

  allocate(ovlp(1:nfun, 1:nfun, 1:nspin))

  call hf_ovlp(.true., rmax, wfn, wfn, ovlp)

  do ispin = 1, nspin
     do ifun = 1, nfun
        do jfun = 1, ifun
           write(io, "(2E20.10)", advance='no') ovlp(jfun, ifun, ispin)
        end do
     end do
  end do
  write(io, *)

  do ispin = 1, nspin
     dets(ispin) = util_det(ne(ispin), thrwfn, ovlp(1:ne(ispin), 1:ne(ispin), ispin))
  end do
  norm = sqrt(dble(dets(1) * dets(nspin)))
    
  deallocate(ovlp)
    
end subroutine hf_print_ovlp
!######################################################################
