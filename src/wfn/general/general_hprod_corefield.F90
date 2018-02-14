subroutine general_hprod_corefield(wfn, fwfn, ene, ng0, ng1)
!
! Complete operator-vector product fwfn, where f is the core Fock.
! Array fwfn, on entry, holds core Hamiltonian contributions, and
! Coulomb and exchange contributions are added here.
!
  use const_mod, only : czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul
  use wfn_mod, only : nfcore, ncore, nfun

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun)     ! orbitals
  complex(kind(0d0)), intent(inout) :: fwfn(0:ngrid, 1:nfun) ! core Fock times orbitals
  real(kind(0d0)), intent(inout) :: ene                      ! core energy
  integer, intent(in) :: ng0, ng1                    ! OMP grid boundaries

  integer :: npt, igrid, jgrid, ifun, jfun
  real(kind(0d0)) :: ene0, ene1
  complex(kind(0d0)) :: rho, tmp
  complex(kind(0d0)), external :: util_zdotu
  complex(kind(0d0)), allocatable :: veff(:), cwfn(:,:)

  if (ncore == 0) return

  npt = ng1 - ng0 + 1
  allocate(veff(ng0:ng1))
  allocate(cwfn(0:ngrid, 1:ncore))

  do jfun = 1, ncore
     do jgrid = gll, gul
        cwfn(jgrid, jfun) = conjg(wfn(jgrid, jfun)) * dgrid
     end do
  end do

  ! bare Hamiltonian expectation value (half of which)
  tmp = czero
  do ifun = 1, ncore
     tmp = tmp + util_zdotu(npt, cwfn(ng0,ifun), 1, fwfn(ng0,ifun), 1)
  end do
  ene0 = dble(tmp)
!debug
!print *, 'general_hprod_corefield: ene0 =', ene0
!debug
  ! Coulomb
  call util_zcopy(npt, czero, 0, veff(ng0), 1)
  do jgrid = gll, gul
     rho = czero
     do jfun = 1, ncore
        rho = rho + cwfn(jgrid, jfun) * wfn(jgrid, jfun)
     end do
     rho = rho + rho

     if (abs(rho) > thrwfn) then
        do igrid = ng0, ng1
           veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
        end do
     end if
  end do

  do ifun = nfcore + 1, nfun
     do igrid = ng0, ng1
        fwfn(igrid, ifun) = fwfn(igrid, ifun) + veff(igrid) * wfn(igrid, ifun)
     end do
  end do

  ! exchange
  do ifun = nfcore + 1, nfun
     do jfun = 1, ncore
        call util_zcopy(npt, czero, 0, veff(ng0), 1)
        do jgrid = gll, gul
           rho = cwfn(jgrid, jfun) * wfn(jgrid, ifun)
           if (abs(rho) > thrwfn) then
              do igrid = ng0, ng1
                 veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
              end do
           end if
        end do

        do igrid = ng0, ng1
           fwfn(igrid, ifun) = fwfn(igrid, ifun) - veff(igrid) * wfn(igrid, jfun)
        end do
     end do
  end do

  ! core Fock expectation value (half of which)
  tmp = czero
  do ifun = 1, ncore
     tmp = tmp + util_zdotu(npt, cwfn(ng0,ifun), 1, fwfn(ng0,ifun), 1)
  end do
  ene1 = dble(tmp)
  ene = ene + ene0 + ene1
!debug
!print *, 'general_hprod_corefield: ene1 =', ene1
!debug

  deallocate(cwfn)
  deallocate(veff)

end subroutine general_hprod_corefield
!######################################################################
