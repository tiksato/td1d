!######################################################################
real(kind(0d0)) function hf_ene2(wfn, hwfn, ng0, ng1)

  use const_mod, only : czero
  use mol_mod, only : ne
  use grid_mod, only : ngrid, dgrid
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfn (0:ngrid,  1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: hwfn(0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1
  integer :: igrid, ifun, ispin
  complex(kind(0d0)) :: tmp

  tmp = czero
  do ispin = 1, nspin
     do ifun = 1, ne(ispin)
        do igrid = ng0, ng1
           tmp = tmp + conjg( wfn(igrid, ifun, ispin)) &
                         & * hwfn(igrid, ifun, ispin)
        end do
     end do
  end do
  hf_ene2 = dble(tmp * dgrid / nspin)

  return

end function hf_ene2
!######################################################################
