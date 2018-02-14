!######################################################################
real(kind(0d0)) function general_ecore(norb, wfn, ng0, ng1)

  use const_mod, only : zero, two, czero
  use grid_mod, only : ngrid, dgrid, v2
  use thresh_mod, only : thrwfn

  implicit none
  integer, intent(in) :: norb, ng0, ng1
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:*)

  integer :: igrid, jgrid, ifun, jfun
  real(kind(0d0)) :: enej, enek
  complex(kind(0d0)) :: rho
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(ng0:ng1))

  !========== Coulomb ==========
  veff(ng0:ng1) = czero
  do jgrid = 0, ngrid
     rho = czero
     do jfun = 1, norb
        rho = rho + conjg(wfn(jgrid, jfun)) * wfn(jgrid, jfun)
     end do
     rho = rho * dgrid * two

     if (abs(rho) > thrwfn) then
        do igrid = ng0, ng1
           veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
        end do
     end if
  end do

  enej = zero
  do igrid = ng0, ng1
     rho = czero
     do ifun = 1, norb
        rho = rho + conjg(wfn(igrid, ifun)) * wfn(igrid, ifun)
     end do
     rho = rho * dgrid
     enej = enej + dble(rho * veff(igrid))
  end do

  !========== Exchange ==========
  enek = zero
  do ifun = 1, norb
     do jfun = 1, norb
        veff(ng0:ng1) = czero
        do jgrid = 0, ngrid
           rho = (conjg(wfn(jgrid, jfun)) &
                    & * wfn(jgrid, ifun)) * dgrid
           if (abs(rho) > thrwfn) then
              do igrid = ng0, ng1
                 veff(igrid) = veff(igrid) + rho * v2(igrid, jgrid)
              end do
           end if
        end do

        do igrid = ng0, ng1
           rho = conjg(wfn(igrid, ifun)) * wfn(igrid, jfun) * dgrid
           enek = enek - dble(rho * veff(igrid))
        end do
     end do
  end do

  general_ecore = enej + enek

  deallocate(veff)

end function general_ecore
!######################################################################
