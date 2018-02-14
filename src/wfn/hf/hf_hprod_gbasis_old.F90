!######################################################################
real(kind(0d0)) function hf_hprod_gbasis_old(dokin, calene, lfield, dtime, wfn0, wfn, hwfn)

  use const_mod, only : zero, czero
  use grid_mod, only : ngbas, goinv, gcore, geris, gmap2, gmap4
  use wfn_mod, only : nfcore, nfun, nspin, hf_doproj, sae

  implicit none
  logical, intent(in) :: dokin
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in)  :: wfn0(1:ngbas, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)  :: wfn (1:ngbas, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: hwfn(1:ngbas, 1:nfun, 1:nspin)

  real(kind(0d0)) :: ene
  complex(kind(0d0)), allocatable :: fock(:,:,:)  ! AO Fock matrix
  complex(kind(0d0)), allocatable :: pmat(:,:,:)  ! AO density matrix: P = CC*
  complex(kind(0d0)), allocatable :: qmat(:,:,:)  ! Q-space projector: Q = S^{-1}(1-SP) = S^{-1} - P
  complex(kind(0d0)), allocatable :: h1wfn(:,:,:)
  integer :: mu,nu,mu2,nu2,ifun,ispin

  if (sae) stop "sae not supported for hf_hprod_gbasis"
  if (.not.dokin) stop "split-operator not supported for hf_hprod_gbasis"

  allocate(fock(1:ngbas, 1:ngbas, 1:nspin))
  allocate(pmat(1:ngbas, 1:ngbas, 1:nspin))
  allocate(qmat(1:ngbas, 1:ngbas, 1:nspin))
  allocate(h1wfn(1:ngbas, 1:nfun, 1:nspin))

  ! Q-space projector
  pmat = czero
  qmat = czero
  !$omp parallel default(shared)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do mu = 1, ngbas
        do nu = 1, ngbas
           do ifun = 1, nfun
              pmat(mu,nu,ispin)=pmat(mu,nu,ispin)+wfn(mu,ifun,ispin)*conjg(wfn(nu,ifun,ispin))
           end do
           qmat(mu,nu,ispin)=goinv(mu,nu)-pmat(mu,nu,ispin)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! Fock matrix
  ene = 0d0
  fock(:,:,1) = gcore
  if (nspin==2) fock(:,:,2) = gcore
  !$omp parallel default(shared) reduction(+:ene)
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
           ene=ene+(gcore(mu,nu)+fock(mu,nu,ispin))*pmat(nu,mu,ispin)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
!BUG  if (nspin==1) ene=ene*2d0
  hf_hprod_gbasis = ene

  h1wfn = czero
  !$omp parallel default(shared)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do ifun = 1, nfun
        do mu = 1, ngbas
           do nu = 1, ngbas
              h1wfn(mu,ifun,ispin)=h1wfn(mu,ifun,ispin)+fock(mu,nu,ispin)*wfn(nu,ifun,ispin)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  if (.not. hf_doproj) then
     hwfn = h1wfn
  else
     hwfn = czero
     !$omp parallel default(shared)
     !$omp do collapse(3)
     do ispin = 1, nspin
        do ifun = 1, nfun
           do mu = 1, ngbas
              do nu = 1, ngbas
                 hwfn(mu,ifun,ispin)=hwfn(mu,ifun,ispin)+qmat(mu,nu,ispin)*h1wfn(nu,ifun,ispin)
              end do
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
  end if

  write(6,"('MO coeff:')")
  do mu=1, ngbas
     write(6,"(i5)",advance='no') mu
     do ifun = 1, nfun
        write(6,"(E20.12)",advance='no') dble(wfn(mu,ifun,1))
     end do
     write(6,*)
  end do
  write(6,"('PMAT:')")
  do mu=1, ngbas
     do nu=1, ngbas
        write(6,"(2i5,E20.12)") mu,nu,dble(pmat(mu,nu,1))
     end do
  end do
  write(6,"('QMAT:')")
  do mu=1, ngbas
     do nu=1, ngbas
        write(6,"(2i5,E20.12)") mu,nu,dble(qmat(mu,nu,1))
     end do
  end do
  write(6,"('FOCK:')")
  do mu=1, ngbas
     do nu=1, ngbas
        write(6,"(2i5,E20.12)") mu,nu,dble(fock(mu,nu,1))
     end do
  end do
  !STOP "for debug @ hf_hprod_gbasis 1."

  deallocate(h1wfn)
  deallocate(qmat)
  deallocate(pmat)
  deallocate(fock)
  return

end function hf_hprod_gbasis_old
!######################################################################
