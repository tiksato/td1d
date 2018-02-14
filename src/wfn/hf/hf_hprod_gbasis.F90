!######################################################################
real(kind(0d0)) function hf_hprod_gbasis(dokin, calene, lfield, dtime, wfn0, wfn, hwfn)

  use root_mod, only : icomp
  use mol_mod, only : ne
  use const_mod, only : zero, one, czero, iunit, runit
  use grid_mod, only : ngbas, gcore, geris, gmap2, gmap4
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
  complex(kind(0d0)) :: tfac
  complex(kind(0d0)), allocatable :: fock(:,:,:)  ! AO Fock matrix
  complex(kind(0d0)), allocatable :: pmat(:,:,:)  ! AO density matrix: P = CC*
  complex(kind(0d0)), allocatable :: whw(:,:,:)
  integer :: mu,nu,mu2,nu2,ifun,jfun,ispin

  if (sae) stop "sae not supported for hf_hprod_gbasis"
  if (.not.dokin) stop "split-operator not supported for hf_hprod_gbasis"

  allocate(fock(1:ngbas, 1:ngbas, 1:nspin))
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
                 fock(mu,nu,ispin)=fock(mu,nu,ispin) &
                      +(2d+0*geris(gmap4(gmap2(mu,nu ),gmap2(mu2,nu2))) &
                            -geris(gmap4(gmap2(mu,nu2),gmap2(mu2,nu ))))*pmat(nu2,mu2,ispin)
              end do
           end do
           ene=ene+(gcore(mu,nu)+fock(mu,nu,ispin))*pmat(nu,mu,ispin)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  if (nspin==2) ene=ene*0.5d0
  hf_hprod_gbasis = ene
!DEBUG
!write(6,"('ene=',f20.10)") ene
!DEBUG

  hwfn = czero
  !$omp parallel default(shared)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do ifun = 1, nfun
        do mu = 1, ngbas
           do nu = 1, ngbas
              hwfn(mu,ifun,ispin)=hwfn(mu,ifun,ispin)+fock(mu,nu,ispin)*wfn(nu,ifun,ispin)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
!DEBUG
!write(6,"('hf_hprod_gbasis: hwfn before Q application:')")
!do mu=1, ngbas
!   write(6,"(i5)",advance='no') mu
!   do ifun = 1, nfun
!      write(6,"(E20.12)",advance='no') dble(hwfn(mu,ifun,1))
!   end do
!   write(6,*)
!end do
!DEBUG

  ! Q-space projection
  if (hf_doproj) then
     !WHY DOESN'T THIS WORK?
     !call hf_hprod_proj(dtime, wfn, hwfn)
     !WHY DOESN'T THIS WORK?
     allocate(whw(1:nfun,1:nfun,1:nspin))
     call hf_ovlp_gbasis(wfn,hwfn,whw)
     do ispin = 1, nspin
        ! occupied orbitals
        do ifun = nfcore + 1, ne(ispin)
           do jfun = 1, ne(ispin)
              tfac = - whw(jfun,ifun,ispin)
              call zaxpy(ngbas,tfac,wfn(1,jfun,ispin),1,hwfn(1,ifun,ispin),1)
           end do
        end do
     
        ! virtual orbitals
        do ifun = ne(ispin) + 1, nfun
           do jfun = 1, ne(ispin)
              tfac = conjg(whw(ifun,jfun,ispin))
              call util_zcopy(ngbas,czero,0,hwfn(1,ifun,ispin),1)
              call zaxpy(ngbas,tfac,wfn(1,jfun,ispin),1,hwfn(1,ifun,ispin),1)
           end do
        end do
     end do
     deallocate(whw)
  end if

  ! multipy time step
  if (icomp == 1) then
     tfac = -iunit*dtime
  else
     tfac = -runit*dtime
  end if
  !$omp parallel default(shared)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do ifun = 1, nfun
        do mu = 1, ngbas
           hwfn(mu,ifun,ispin) = tfac*hwfn(mu,ifun,ispin)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  !write(6,"('MO coeff:')")
  !do mu=1, ngbas
  !   write(6,"(i5)",advance='no') mu
  !   do ifun = 1, nfun
  !      write(6,"(E20.12)",advance='no') dble(wfn(mu,ifun,1))
  !   end do
  !   write(6,*)
  !end do
  !write(6,"('MO coeff derivative:')")
  !do mu=1, ngbas
  !   write(6,"(i5)",advance='no') mu
  !   do ifun = 1, nfun
  !      write(6,"(E20.12)",advance='no') dble(hwfn(mu,ifun,1))
  !   end do
  !   write(6,*)
  !end do
  !stop
  !write(6,"('PMAT:')")
  !do mu=1, ngbas
  !   do nu=1, ngbas
  !      write(6,"(2i5,E20.12)") mu,nu,dble(pmat(mu,nu,1))
  !   end do
  !end do
  !write(6,"('FOCK:')")
  !do mu=1, ngbas
  !   do nu=1, ngbas
  !      write(6,"(2i5,E20.12)") mu,nu,dble(fock(mu,nu,1))
  !   end do
  !end do
  !STOP "for debug @ hf_hprod_gbasis 1."

  deallocate(pmat)
  deallocate(fock)
  return

end function hf_hprod_gbasis
!######################################################################
