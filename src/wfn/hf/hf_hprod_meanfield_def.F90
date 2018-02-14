!######################################################################
subroutine hf_hprod_meanfield_def(wfnin, wfn, hwfn, ng0, ng1)

  use mol_mod, only : ne
  use const_mod, only : zero, one, two, czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul
  use wfn_mod, only : nfun, nspin, saex, crapx, ldacx

  implicit none
  complex(kind(0d0)), intent(in)    :: wfnin(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)    :: wfn  (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(inout) :: hwfn (0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid, ifun, jfun, ispin
  real(kind(0d0)) :: rhoj, rwfn, iwfn
  complex(kind(0d0)) :: rhok
  real(kind(0d0)), allocatable :: veffj(:)
  complex(kind(0d0)), allocatable :: veffk(:)
  real(kind(0d0)) :: hffac

  if (nspin == 2) then
     hffac = dgrid
  else
     hffac = dgrid * two
  end if

  !========== Coulomb ==========
  allocate(veffj(ng0:ng1))
  veffj(ng0:ng1) = zero
  do jgrid = gll, gul
     rhoj = zero
     do ispin = 1, nspin
        do jfun = 1, ne(ispin)
           rwfn = real (wfnin(jgrid, jfun, ispin))
           iwfn = aimag(wfnin(jgrid, jfun, ispin))
           rhoj = rhoj + rwfn * rwfn + iwfn * iwfn
        end do
     end do
     rhoj = rhoj * hffac

     if (abs(rhoj) > thrwfn) then
        do igrid = ng0, ng1
           veffj(igrid) = veffj(igrid) + rhoj * v2(igrid, jgrid)
        end do
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

  if (ldacx) return

  !========== Exchange ==========
  allocate(veffk(ng0:ng1))
  do ispin = 1, nspin
     do ifun = 1, nfun
        do jfun = 1, ne(ispin)
           veffk(ng0:ng1) = czero
           do jgrid = gll, gul
              rhok = (conjg(wfnin(jgrid, jfun, ispin)) &
                        & * wfn  (jgrid, ifun, ispin)) * dgrid
              if (abs(rhok) > thrwfn) then
                 do igrid = ng0, ng1
                    veffk(igrid) = veffk(igrid) + rhok * v2(igrid, jgrid)
                 end do
              end if
           end do

           do igrid = ng0, ng1
              hwfn(igrid, ifun, ispin) &
          & = hwfn(igrid, ifun, ispin) - veffk(igrid) * wfnin(igrid, jfun, ispin)
           end do
        end do
     end do
  end do
  deallocate(veffk)

end subroutine hf_hprod_meanfield_def
!######################################################################
