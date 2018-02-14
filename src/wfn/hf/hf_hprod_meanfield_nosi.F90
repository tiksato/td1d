!######################################################################
subroutine hf_hprod_meanfield_nosi(wfnin, wfn, hwfn, ng0, ng1)

  use mol_mod, only : ne
  use const_mod, only : zero, one, two, czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(in)    :: wfnin(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)    :: wfn  (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(inout) :: hwfn (0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid, ifun, jfun, ispin, jspin
  real(kind(0d0)) :: rhoj, rwfn, iwfn, hffac
  real(kind(0d0)), allocatable :: veffj(:,:,:), totj(:)
  complex(kind(0d0)) :: rhok
  complex(kind(0d0)), allocatable :: veffk(:)

  if (nspin == 2) then
     hffac = one
  else
     hffac = two
  end if

  !========== Coulomb ==========
  allocate(totj (0:ngrid))
  allocate(veffj(0:ngrid, 1:nfun, 1:nspin))

  totj(ng0:ng1) = zero
  veffj(ng0:ng1, 1:nfun, 1:nspin) = zero

  do jspin = 1, nspin
     do jfun = 1, ne(jspin)
        do jgrid = gll, gul
           rwfn = dble (wfnin(jgrid, jfun, jspin))
           iwfn = aimag(wfnin(jgrid, jfun, jspin))
           rhoj = (rwfn * rwfn + iwfn * iwfn) * dgrid
           if (abs(rhoj) > thrwfn) then
              do igrid = ng0, ng1
                 veffj(igrid, jfun, jspin) &
             & = veffj(igrid, jfun, jspin) + rhoj * v2(igrid, jgrid)
              end do
           end if
        end do
        totj(ng0:ng1) = totj(ng0:ng1) + veffj(ng0:ng1, jfun, jspin) * hffac
     end do
  end do

  ! remove self-interaction term
  do ispin = 1, nspin
     do ifun = 1, nfun
        veffj(ng0:ng1, ifun, ispin) = totj(ng0:ng1) - veffj(ng0:ng1, ifun, ispin)
     end do
  end do

  do ispin = 1, nspin
     do ifun = 1, nfun
        do igrid = ng0, ng1
           hwfn(igrid, ifun, ispin) &
       & = hwfn(igrid, ifun, ispin) + veffj(igrid, ifun, ispin) * wfn (igrid, ifun, ispin)
        end do
     end do
  end do
  deallocate(totj)
  deallocate(veffj)

  !========== Exchange ==========
  allocate(veffk(0:ngrid))
  do ispin = 1, nspin
     do ifun = 1, nfun
        do jfun = 1, ne(ispin)
           ! skip self-interaction term
           if (jfun == ifun) cycle

           ! exchange kernel
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

           ! applying exchange operator
           do igrid = ng0, ng1
              hwfn(igrid, ifun, ispin) &
          & = hwfn(igrid, ifun, ispin) - veffk(igrid) * wfnin(igrid, jfun, ispin)
           end do
        end do
     end do
  end do
  deallocate(veffk)

end subroutine hf_hprod_meanfield_nosi
!######################################################################
