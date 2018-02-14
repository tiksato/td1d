!######################################################################
subroutine hf_fulldiag(lfield, wfn, fullwfn)

  use omp_mod
  use const_mod, only : czero, one
  use grid_mod, only : ngrid, dgrid, gll, gul, gbasis
  use wfn_mod, only : nfun, nspin

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: fullwfn(0:ngrid, 0:ngrid, 1:nspin)

  complex(kind(0d0)), allocatable :: fock(:,:,:)
  integer :: ispin

  if (gbasis) then
     call hf_fulldiag_gbasis(lfield, wfn, fullwfn)
     return
  end if

  allocate(fock(0:ngrid, 0:ngrid, 1:nspin))
  call util_zcopy((ngrid+1)*(ngrid+1)*nspin, czero, 0, fock, 1)

!$omp parallel default(shared)

  call omp_mod_thread(gll, gul)

  call hf_full_kinetic(fock, ng0, ng1)
  call hf_full_potential(lfield, fock, ng0, ng1)
  call hf_full_meanfield(wfn, fock, ng0, ng1)

!$omp end parallel

  call hf_full_diag(fock, fullwfn)
  do ispin = 1, nspin
     call zscal((ngrid+1)*nfun, one/sqrt(dgrid), fullwfn(0,0,ispin), 1)
     call util_zcopy((ngrid+1)*nfun, fullwfn(0,0,ispin), 1, wfn(0,1,ispin), 1)
  end do

  deallocate(fock)

end subroutine hf_fulldiag
!######################################################################
subroutine hf_fulldiag_gbasis(lfield, wfn, fullwfn)

  use mol_mod, only : ne
  use const_mod, only : zero, one, czero, iunit, runit
  use grid_mod, only : ngbas, gcore, geris, gmap2, gmap4
  use wfn_mod, only : nfcore, nfun, nspin

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(inout) :: wfn(1:ngbas, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: fullwfn(1:ngbas, 1:ngbas, 1:nspin)

  real(kind(0d0)) :: ene
  complex(kind(0d0)) :: tfac
  complex(kind(0d0)), allocatable :: fock(:,:,:)  ! AO Fock matrix
  integer :: mu,nu,mu2,nu2,ifun,jfun,ispin

  allocate(fock(1:ngbas, 1:ngbas, 1:nspin))
  call hf_mkfock_gbasis(lfield,wfn,fock)

  call hf_full_diag(fock, fullwfn)
  do ispin = 1, nspin
     call util_zcopy(ngbas*nfun,fullwfn(1,1,ispin),1,wfn(1,1,ispin),1)
  end do

  deallocate(fock)

end subroutine hf_fulldiag_gbasis
!######################################################################
subroutine hf_mkfock_gbasis(lfield, wfn, fock)

  use mol_mod, only : ne
  use const_mod, only : zero, one, czero, iunit, runit
  use grid_mod, only : ngbas, gcore, geris, gmap2, gmap4
  use wfn_mod, only : nfcore, nfun, nspin

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(inout) :: wfn(1:ngbas, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: fock(1:ngbas, 1:ngbas, 1:nspin)

  real(kind(0d0)) :: ene
  complex(kind(0d0)) :: tfac
  complex(kind(0d0)), allocatable :: pmat(:,:,:)  ! AO density matrix: P = CC*
  integer :: mu,nu,mu2,nu2,ifun,jfun,ispin

  allocate(pmat(1:ngbas, 1:ngbas, 1:nspin))
  ! AO spin density
  pmat = czero
  !$omp parallel default(shared)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do mu = 1, ngbas
        do nu = 1, ngbas
           do jfun = 1, ne(ispin)
              pmat(mu,nu,ispin)=pmat(mu,nu,ispin)+wfn(mu,jfun,ispin)*conjg(wfn(nu,jfun,ispin))
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! Fock matrix
  fock(:,:,1) = gcore
  if (nspin==2) fock(:,:,2) = gcore
  !$omp parallel default(shared)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do mu = 1, ngbas
        do nu = 1, ngbas
           do mu2 = 1, ngbas
              do nu2 = 1, ngbas
                 fock(mu,nu,ispin)=fock(mu,nu,ispin)+(2d+0*geris(gmap4(gmap2(mu,nu ),gmap2(mu2,nu2))) &
                                                          -geris(gmap4(gmap2(mu,nu2),gmap2(mu2,nu ))))*pmat(nu2,mu2,ispin)
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(pmat)

end subroutine hf_mkfock_gbasis
!######################################################################
