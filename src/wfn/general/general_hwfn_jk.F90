!######################################################################
subroutine general_hwfn_jk(nfin, nfout, wfn1, wfn2, wfnket, scale, &
     & hwfn, ng0, ng1)

  use const_mod, only : czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul

  implicit none
  integer, intent(in) :: nfin, nfout
  complex(kind(0d0)), intent(in) :: wfn1(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: wfn2(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: wfnket(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: scale
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:*)
  integer, intent(in) :: ng0, ng1

  integer :: npt, igrid, jgrid, ifun, jfun
  complex(kind(0d0)) :: rho, fac
  complex(kind(0d0)), allocatable :: veff(:)

  fac = dgrid * scale
  npt = ng1 - ng0 + 1
  allocate(veff(ng0:ng1))

  !========== Coulomb ==========
  call util_zcopy(npt, czero, 0, veff(ng0), 1)
  do jgrid = gll, gul
     rho = czero
     do jfun = 1, nfin
        rho = rho + conjg(wfn1(jgrid, jfun)) &
                      & * wfn2(jgrid, jfun) * fac
     end do
     rho = rho + rho

     if (abs(rho) > thrwfn) then
        do igrid = ng0, ng1
           veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
        end do
     end if
  end do

  do ifun = 1, nfout
     do igrid = ng0, ng1
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + veff(igrid) * wfnket(igrid, ifun)
     end do
  end do

  !========== Exchange ==========
  do ifun = 1, nfout
     do jfun = 1, nfin
        call util_zcopy(npt, czero, 0, veff(ng0), 1)
        do jgrid = gll, gul
           rho = conjg(wfn1(jgrid, jfun)) &
                   * wfnket(jgrid, ifun) * fac
           if (abs(rho) > thrwfn) then
              do igrid = ng0, ng1
                 veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
              end do
           end if
        end do

        do igrid = ng0, ng1
           hwfn(igrid, ifun) = hwfn(igrid, ifun) - veff(igrid) * wfn2(igrid, jfun)
        end do
     end do
  end do

  deallocate(veff)

end subroutine general_hwfn_jk
!######################################################################
