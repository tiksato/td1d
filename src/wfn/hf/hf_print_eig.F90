!######################################################################
subroutine hf_print_eig(wfn)

  use omp_mod
  use io_mod, only : iostdo
  use const_mod, only : czero, zero
  use grid_mod, only : ngrid, dgrid, gll, gul, gbasis
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  complex(kind(0d0)), allocatable :: fock(:,:,:)
  complex(kind(0d0)), allocatable :: full(:,:,:)
  integer :: ispin, ifun

  allocate(fock(0:ngrid, 0:ngrid, 1:nspin))
  allocate(full(0:ngrid, 0:ngrid, 1:nspin))
  call util_zcopy((ngrid+1)*(ngrid+1)*nspin, czero, 0, fock, 1)
  call util_zcopy((ngrid+1)*(ngrid+1)*nspin, czero, 0, full, 1)

  if (.not. gbasis) then
     !$omp parallel default(shared)
     call omp_mod_thread(gll, gul)
     call hf_full_kinetic(fock, ng0, ng1)
     call hf_full_potential(zero, fock, ng0, ng1)
     call hf_full_meanfield(wfn, fock, ng0, ng1)
     !$omp end parallel
  else
     call hf_mkfock_gbasis(zero, wfn, fock)
  end if

  call hf_full_diag(fock(0,0,1), full(0,0,1))
  write(iostdo, "('# orbital energies:')")
  do ifun = 0, nfun - 1
     write(iostdo, "(21x,i10)", advance="no") ifun + 1
     write(iostdo, "(8x,f20.10)", advance="no") dble(fock(ifun, ifun, 1))
     if (nspin == 2) write(iostdo, "(f20.10)", advance="no") dble(fock(ifun, ifun, 2))
     write(iostdo, *)
  end do

  deallocate(full)
  deallocate(fock)

end subroutine hf_print_eig
!######################################################################
