!######################################################################
real(kind(0d0)) function x2e_hprod(dokin, calene, lfield, wfn, hwfn)

  use omp_mod
  use const_mod, only : zero
  use grid_mod, only : ngrid, gll, gul
  use fft_mod, only : dofft, fft_wfnr, fft_wfnk

  implicit none
  logical, intent(in) :: dokin
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)  :: wfn (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: hwfn(0:ngrid, 0:ngrid)

  real(kind(0d0)), external :: x2e_ene
  integer :: igrid, jgrid
  real(kind(0d0)) :: ene

  ene = zero
  call x2e_clear(hwfn)

  if (dokin .and. dofft) call x2e_hprod_kinetic_fft(lfield, wfn, hwfn, fft_wfnr, fft_wfnk)

!$omp parallel default(shared) reduction(+:ene)
  call omp_mod_thread(gll, gul)

  if (dokin .and. .not. dofft) call x2e_hprod_kinetic(lfield, wfn, hwfn, ng0, ng1)
  call x2e_hprod_potential(lfield, wfn, hwfn, ng0, ng1)

  if (calene) then
     ene = ene + x2e_ene(wfn, hwfn, ng0, ng1)
  end if
!$omp end parallel

  ! symmetrization
  do igrid = 0, ngrid
     do jgrid = igrid + 1, ngrid
        hwfn(jgrid, igrid) = hwfn(igrid, jgrid)
     end do
  end do

  x2e_hprod = ene
  return

end function x2e_hprod
!######################################################################
