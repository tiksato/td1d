!######################################################################
subroutine hf_ort(wfn)

  use const_mod, only : one
  use grid_mod, only : ngrid, gbasis
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  integer :: ifun, jfun, ispin
  real(kind(0d0)) :: norm
  complex(kind(0d0)) :: ovlp
  complex(kind(0d0)), external :: hf_ovlp1

  if (gbasis) then
     call hf_ort_gbasis(wfn)
     return
  end if

  do ispin = 1, nspin
     do ifun = 1, nfun
        do jfun = 1, ifun - 1
           ovlp = hf_ovlp1(-one, wfn(0, jfun, ispin), wfn(0, ifun, ispin))
           wfn(0:ngrid, ifun, ispin) = wfn(0:ngrid, ifun, ispin) - wfn(0:ngrid, jfun, ispin) * ovlp
        end do
        norm = real(hf_ovlp1(-one, wfn(0, ifun, ispin), wfn(0, ifun, ispin)))
        wfn(0:ngrid, ifun, ispin) = wfn(0:ngrid, ifun, ispin) * (one / sqrt(norm))
     end do
  end do

end subroutine hf_ort
!######################################################################
subroutine hf_ort_gbasis(wfn)

  use grid_mod, only : ngbas
  use wfn_mod, only : nfun,nspin

  implicit none
  complex(kind(0d0)), intent(inout) :: wfn(1:ngbas,1:nfun,1:nspin)

  integer :: ifun,jfun,ispin
  real(kind(0d0)) :: norm
  complex(kind(0d0)) :: ovlp
  complex(kind(0d0)), external :: hf_ovlp1_gbasis

  do ispin = 1, nspin
     do ifun = 1, nfun
        do jfun = 1, ifun - 1
           ovlp = hf_ovlp1_gbasis(wfn(1,jfun,ispin),wfn(1,ifun,ispin))
           wfn(1:ngbas,ifun,ispin) = wfn(1:ngbas,ifun,ispin) - wfn(1:ngbas,jfun,ispin)*ovlp
!DEBUG
!write(6,"('hf_ort_gbasis: ovlp ',2i5,2f20.10)") ifun,jfun,ovlp
!DEBUG
        end do
        norm = dble(hf_ovlp1_gbasis(wfn(1,ifun,ispin),wfn(1,ifun,ispin)))
        wfn(1:ngbas,ifun,ispin) = wfn(1:ngbas,ifun,ispin)/sqrt(norm)
!DEBUG
!write(6,"('hf_ort_gbasis: ovlp ',2i5,2f20.10)") ifun,ifun,norm,0d0
!DEBUG
     end do
  end do

end subroutine hf_ort_gbasis
!######################################################################
