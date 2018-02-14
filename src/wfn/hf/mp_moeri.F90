!!!!! Closed-shell only !!!!!
!######################################################################
subroutine mp_moeri(nvir, wfnocc, wfnvir, hwfn, ng0, ng1)

  use mol_mod, only : ne
  use const_mod, only : one, two, czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul
  use wfn_mod, only : nfun

  implicit none
  integer, intent(in) :: nvir
  complex(kind(0d0)), intent(in)    :: wfnocc(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(in)    :: wfnvir(0:ngrid, 1:nvir)
  complex(kind(0d0)), intent(inout) :: hwfn  (0:ngrid, 1:nvir)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid, ifun, jfun
  complex(kind(0d0)) :: rho, test
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(ng0:ng1))
    
  do ifun = 1, nvir
     do jfun = 1, nfun

        ! kernel
        veff(ng0:ng1) = czero
        do jgrid = gll, gul
           rho = (conjg(wfnvir(jgrid, ifun)) &
                    & * wfnocc(jgrid, jfun)) * dgrid
           if (abs(rho) > thrwfn) then
              do igrid = ng0, ng1
                 veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
              end do
           end if
        end do

        ! operator
        do igrid = ng0, ng1
           hwfn(igrid, ifun) &
       & = hwfn(igrid, ifun) + veff(igrid) * wfnocc(igrid, jfun)
        end do
     end do
  end do

  deallocate(veff)

end subroutine mp_moeri
!######################################################################
