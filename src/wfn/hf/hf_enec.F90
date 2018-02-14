!######################################################################
real(kind(0d0)) function hf_enec(wfn, hwfn, ng0, ng1)

  use const_mod, only : one, two, czero
  use mol_mod, only : ne
  use grid_mod, only : ngrid, dgrid
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfn (0:ngrid,  1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: hwfn(0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1
  integer :: igrid, ifun, ispin
  complex(kind(0d0)) :: tmp
  real(kind(0d0)) :: hffac

  if (nspin == 2 .or. ne(3) == 1) then
     hffac = dgrid
  else
     hffac = dgrid * two
  end if

  tmp = czero
  do ispin = 1, nspin
     do ifun = 1, ne(ispin)
        do igrid = ng0, ng1
           tmp = tmp + conjg(wfn(igrid, ifun, ispin)) &
                   & *      hwfn(igrid, ifun, ispin)
        end do
     end do
  end do
  hf_enec = dble(tmp * hffac)

  return

end function hf_enec
!######################################################################
