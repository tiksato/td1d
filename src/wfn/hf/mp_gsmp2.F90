!######################################################################
subroutine mp_gsmp2(lfield, wfn0, wfn1)

  use omp_mod
  use io_mod, only : iostdo
  use const_mod, only : czero, runit, zero, one, two
  use grid_mod, only : ngrid, dgrid, gll, gul
  use wfn_mod, only : nfun

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)  :: wfn0(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(out) :: wfn1(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid, nvir, ifun, jfun
  real(kind(0d0)), external :: x2e_hprod
  complex(kind(0d0)), external :: general_ovlp1
  complex(kind(0d0)), allocatable :: fock(:,:)
  complex(kind(0d0)), allocatable :: full(:,:)
  complex(kind(0d0)), allocatable :: wfnvir(:,:)
  complex(kind(0d0)), allocatable :: hwfn(:,:)
  complex(kind(0d0)), allocatable :: twfn(:,:)
  complex(kind(0d0)), allocatable :: tab(:,:)
  real(kind(0d0)) :: ene0, ene01, ene2, ene3
  complex(kind(0d0)) :: abij, ijab, denom, tmp, wfac, norm
!debug
!  complex(kind(0d0)), allocatable :: focktmp(:,:)
!  complex(kind(0d0)), allocatable :: eigvir(:)
!debug

  nvir = ngrid + 1 - nfun

  allocate(fock  (0:ngrid, 0:ngrid))
  allocate(full  (0:ngrid, 0:ngrid))
  allocate(wfnvir(0:ngrid, 1:nvir))
  allocate(hwfn  (0:ngrid, 1:nvir))
  allocate(twfn  (0:ngrid, 1:nvir))
  allocate(tab   (1:nvir, 1:nvir))
  call util_zcopy((ngrid+1)*(ngrid+1), czero, 0, wfn1, 1)
  call util_zcopy((ngrid+1)*(ngrid+1), czero, 0, fock, 1)
  call util_zcopy((ngrid+1)*(ngrid+1), czero, 0, full, 1)
  call util_zcopy((ngrid+1)*nvir,      czero, 0, hwfn, 1)
  call util_zcopy((ngrid+1)*nvir,      czero, 0, twfn, 1)
  call util_zcopy(nvir*nvir,           czero, 0, tab,  1)

! ===== (1) construct fock matrix =====

!$omp parallel default(shared)

  call omp_mod_thread(gll, gul)
  call hf_full_kinetic(fock, ng0, ng1)
  call hf_full_potential(lfield, fock, ng0, ng1)
  call hf_full_meanfield(wfn0, fock, ng0, ng1)

!$omp end parallel

! ===== (2) diagonalize fock matrix =====

!debug
!  allocate(focktmp(0:ngrid,0:ngrid))
!  focktmp(0:ngrid,0:ngrid) = fock(0:ngrid,0:ngrid)
!  call hf_full_diag(focktmp, full)
!  ene0 = two * focktmp(0,0)
!  deallocate(focktmp)
!  allocate(eigvir(1:nvir))
!  call mrmp_fockdiag(nvir, wfn0, fock, wfnvir, eigvir)
!  deallocate(eigvir)  
!debug
  call hf_full_diag(fock, full)
  wfac = runit / sqrt(dgrid)
  call zscal((ngrid+1)*(ngrid+1), wfac, full(0,0), 1)
  call util_zcopy((ngrid+1)*nvir, full(0,nfun), 1, wfnvir(0,1), 1)
  ene0 = two * dble(fock(0,0))
!debug
!write(6,"('eigenvalues:')")
!do igrid = 0, 20
!   write(6,"(i20,2f20.10)") igrid, fock(igrid,igrid)
!end do
!debug

! ===== (3) compute two-electron integrals =====

!$omp parallel default(shared)

  call omp_mod_thread(gll, gul)
  call mp_moeri(nvir, wfn0, wfnvir, hwfn, ng0, ng1)

!$omp end parallel

! ===== (3) MP1 amplitudes and MP2 energy =====

  ene2 = zero
  norm = czero
  do ifun = 1, nvir
     do jfun = 1, nvir
        abij = general_ovlp1(-one, wfnvir(0,ifun), hwfn(0,jfun))
        ijab = conjg(abij)
        denom = one / (ene0 - fock(ifun,ifun) - fock(jfun,jfun))
        tab(ifun, jfun) = abij * denom
        ene2 = ene2 + dble(ijab * tab(ifun, jfun))
        norm = norm + conjg(tab(ifun,jfun)) * tab(ifun,jfun)
     end do
  end do
!debug
  write(iostdo, "('# weight of wfnpt:                     ', f20.10)") norm
!debug

! ===== (3) normalized (0-1)th order wavefunction =====
  
  do ifun = 1, nvir
     do jfun = 1, nvir
        do igrid = gll, gul
           twfn(igrid, ifun) = twfn  (igrid, ifun) &
                           & + wfnvir(igrid, jfun) * tab(jfun, ifun)
        end do
     end do
  end do

  call util_zcopy((ngrid+1)*(ngrid+1), czero, 0, full, 1)
  do igrid = gll, gul
     do jgrid = gll, gul
        tmp = wfn0(igrid, 1) * wfn0(jgrid, 1)
        full(igrid, jgrid) = full(igrid, jgrid) + tmp
     end do
  end do
  call util_zcopy((ngrid+1)*(ngrid+1), czero, 0, fock, 1)
  ene01 = x2e_hprod(.true., .true., zero, full, fock)

  do ifun = 1, nvir
     do igrid = gll, gul
        do jgrid = gll, gul
           tmp = wfnvir(igrid, ifun) * twfn(jgrid, ifun)
           full(igrid, jgrid) = full(igrid, jgrid) + tmp
           wfn1(igrid, jgrid) = wfn1(igrid, jgrid) + tmp
        end do
     end do
  end do

  norm = one / sqrt(one + norm)
  call zscal((ngrid+1)*(ngrid+1), norm, full, 1)
  call util_zcopy((ngrid+1)*(ngrid+1), czero, 0, fock, 1)
  ene3 = x2e_hprod(.true., .true., zero, full, fock)

  write(iostdo, "('# second-order energy:                 ', f20.10)") ene2
  write(iostdo, "('# 0th order energy:                    ', f20.10)") ene0
  write(iostdo, "('# scf energy:                          ', f20.10)") ene01
  write(iostdo, "('# mp2 energy:                          ', f20.10)") ene01 + ene2
  write(iostdo, "('# energy of wfn + wfnpt:               ', f20.10)") ene3
!debug
!  call x2e_den2e(full, fock)
!  call x2e_print(iostdo, fock)
!debug

  deallocate(tab)
  deallocate(twfn)
  deallocate(hwfn)
  deallocate(wfnvir)
  deallocate(full)
  deallocate(fock)

end subroutine mp_gsmp2
!######################################################################
