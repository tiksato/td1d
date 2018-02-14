!######################################################################
subroutine hf_hprod_meanfield_lda(wfnin, wfn, hwfn, ng0, ng1)

  use mol_mod, only : ne
  use const_mod, only : zero, one, two, czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul
  use wfn_mod, only : nfun, nspin, saex, crapx

  implicit none
  complex(kind(0d0)), intent(in)    :: wfnin(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)    :: wfn  (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(inout) :: hwfn (0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid, ifun, jfun, ispin
  real(kind(0d0)) :: rho, rwfn, iwfn, veff
  real(kind(0d0)) :: hffac

  if (nspin == 2) then
     hffac = dgrid
  else
     hffac = dgrid * two
  end if

  do igrid = ng0, ng1
     rho = zero
     do ispin = 1, nspin
        do jfun = 1, ne(ispin)
           rwfn = real (wfnin(igrid, jfun, ispin))
           iwfn = aimag(wfnin(igrid, jfun, ispin))
           rho = rho + rwfn * rwfn + iwfn * iwfn
        end do
     end do
     rho = rho * hffac

     if (abs(rho) > thrwfn) then
        veff = rho ** f34
     end if
  end do

  do ispin = 1, nspin
     do ifun = 1, nfun
        do igrid = ng0, ng1
           hwfn(igrid, ifun, ispin) &
       & = hwfn(igrid, ifun, ispin) + veffj(igrid) * wfn(igrid, ifun, ispin) 
        end do
     end do
  end do
  deallocate(veffj)

end subroutine hf_hprod_meanfield_lda
!######################################################################
