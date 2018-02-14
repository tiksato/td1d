!######################################################################
real(kind(0d0)) function hf_hprod(dokin, calene, lfield, dtime, wfn0, wfn, hwfn)

  use omp_mod
  use const_mod, only : zero, czero
  use grid_mod, only : ngrid, gll, gul, gbasis
  use wfn_mod, only : nfcore, nfun, nspin, sae
  use fft_mod, only : dofft

  implicit none
  logical, intent(in) :: dokin
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in)  :: wfn0(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)  :: wfn (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: hwfn(0:ngrid, 1:nfun, 1:nspin)

  real(kind(0d0)) :: ene
  real(kind(0d0)), external :: hf_hprod_sum
  real(kind(0d0)), external :: hf_hprod_gbasis
  complex(kind(0d0)), allocatable :: wfnin(:,:,:)
  complex(kind(0d0)), allocatable :: h1wfn(:,:,:)
  complex(kind(0d0)), allocatable :: h2wfn(:,:,:)

  if (gbasis) then
     !hf_hprod = hf_hprod_gbasis_old(dokin, calene, lfield, dtime, wfn0, wfn, hwfn)
     hf_hprod = hf_hprod_gbasis(dokin, calene, lfield, dtime, wfn0, wfn, hwfn)
     return
  end if

  allocate(wfnin(0:ngrid, 1:nfun, 1:nspin))
  allocate(h1wfn(0:ngrid, 1:nfun, 1:nspin))
  allocate(h2wfn(0:ngrid, 1:nfun, 1:nspin))

  ene = zero
  call hf_clear(hwfn)
  call hf_clear(h1wfn)
  call hf_clear(h2wfn)
  call hf_setfroz(wfn0, wfn, wfnin)

  if (dokin .and. dofft) call general_hprod_kinetic_fft(nfun, nfcore, nspin, lfield, wfn, h1wfn)

!$omp parallel default(shared) reduction(+:ene)

  call omp_mod_thread(gll, gul)

  if (dokin .and. .not. dofft) then
     call general_hprod_kinetic(lfield, wfn(0,1,1), h1wfn(0,1,1), ng0, ng1)
     if (nspin == 2) call general_hprod_kinetic(lfield, wfn(0,1,2), h1wfn(0,1,2), ng0, ng1)
  end if

  call general_hprod_potential(lfield, wfn(0,1,1), h1wfn(0,1,1), ng0, ng1)
  if (nspin == 2) call general_hprod_potential(lfield, wfn(0,1,2), h1wfn(0,1,2), ng0, ng1)

  call hf_hprod_meanfield(wfnin, wfn, h2wfn, ng0, ng1)
  ene = ene + hf_hprod_sum(calene, wfn, h1wfn, h2wfn, hwfn, ng0, ng1)

!$omp end parallel

  call hf_hprod_proj(dtime, wfn, hwfn)
  if (sae) hwfn(0:ngrid, 1:nfun, 1) = czero

  deallocate(h2wfn)
  deallocate(h1wfn)
  deallocate(wfnin)

  hf_hprod = ene
  return

end function hf_hprod
!######################################################################
