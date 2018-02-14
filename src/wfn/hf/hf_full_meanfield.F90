!######################################################################
subroutine hf_full_meanfield(wfn, fock, ng0, ng1)

  use mol_mod, only : ne
  use root_mod, only : icomp, nocoulomb
  use const_mod, only : zero, one, two, czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul
  use wfn_mod, only : nfun, nspin, hcore

  implicit none
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:nfun,  1:nspin)
  complex(kind(0d0)), intent(inout) :: fock(0:ngrid, 0:ngrid, 1:nspin)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid, jfun, ispin
  real(kind(0d0)) :: rhoj, rwfn, iwfn
  complex(kind(0d0)) :: rhok
  real(kind(0d0)) :: hffac

  if (icomp < 0 .or. ne(3) <= 1 .or. nocoulomb .or. hcore) return

  if (nspin == 2) then
     hffac = dgrid
  else
     hffac = dgrid * two
  end if

  !========== Coulomb ==========
  do jgrid = gll, gul
     rhoj = zero
     do ispin = 1, nspin
        do jfun = 1, ne(ispin)
           rwfn = real (wfn(jgrid, jfun, ispin))
           iwfn = aimag(wfn(jgrid, jfun, ispin))
           rhoj = rhoj + rwfn * rwfn + iwfn * iwfn
        end do
     end do
     rhoj = rhoj * hffac
!debug     rhoj = rhoj * dgrid

     if (abs(rhoj) > thrwfn) then
        do ispin = 1, nspin
           do igrid = ng0, ng1
              fock(igrid, igrid, ispin) &
          & = fock(igrid, igrid, ispin) + rhoj * v2(igrid, jgrid)
           end do
        end do
     end if
  end do

!debug  return

  !========== Exchange ==========
  do ispin = 1, nspin
     do jfun = 1, ne(ispin)
        do jgrid = gll, gul
           rhok = conjg(wfn(jgrid, jfun, ispin)) * dgrid
           if (abs(rhok) > thrwfn) then
              do igrid = ng0, ng1
                 fock(igrid, jgrid, ispin) &
             & = fock(igrid, jgrid, ispin) - rhok * wfn(igrid, jfun, ispin) * v2(igrid, jgrid)
              end do
           end if
        end do
     end do
  end do

end subroutine hf_full_meanfield
!######################################################################
