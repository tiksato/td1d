!######################################################################
subroutine hf_print_vlocal(iout, istep, time, lfield, wfn)

  use wfn_mod, only : ehf, nfun, nspin
  use root_mod, only : iprint
  use field_mod, only : period
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, x, v1
  use const_mod, only : czero, one, half, runit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: iout, istep
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:*)
  !--------------------------------------------------------------------
  integer :: igrid, ifun
  real(kind(0d0)) :: optcyc
  complex(kind(0d0)), allocatable :: vcoul(:,:)

  allocate(vcoul(0:ngrid, 1:nfun))

  vcoul(0:ngrid, 1:nfun) = czero
  call hf_print_vlocal_coul(wfn, vcoul)

  ! print
  optcyc = dble(time) / period
  write(iout, "('# vlocal:     ', i15, E15.6)") istep, optcyc
  do igrid = 0, ngrid
     write(iout, "(E15.6)", advance = 'no') x(igrid)
     write(iout, "(E15.6)", advance = 'no') v1(igrid)
     write(iout, "(E15.6)", advance = 'no') -lfield * x(igrid)
     do ifun = 1, nfun
        write(iout, "(2E15.6)", advance = 'no') vcoul(igrid, ifun)
     end do
     write(iout, *)
  end do
  write(iout, *)

  deallocate(vcoul)
  
end subroutine hf_print_vlocal
!######################################################################
!######################################################################
subroutine hf_print_vlocal_coul(wfnin, vcoul)

  use omp_mod
  use root_mod, only : nocoulomb
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2
  use wfn_mod, only : nfun, crapola, nfroz
  use const_mod, only : czero, runit, ctwo

  implicit none
  complex(kind(0d0)), intent(in)    :: wfnin(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: vcoul(0:ngrid, 1:nfun)

  integer :: igrid, jgrid, ifun, jfun
  complex(kind(0d0)) :: rho, fac

  !$omp parallel default(shared) private(ifun, jfun, igrid, jgrid, rho)
  call omp_mod_thread(0, ngrid)
  do ifun = 1, nfun
     do jfun = 1, nfun
        if (ifun == jfun) then
           fac = runit
        else
           fac = ctwo
        end if
        do jgrid = 0, ngrid
           rho = (conjg(wfnin(jgrid, jfun)) &
             & *        wfnin(jgrid, jfun)) * dgrid
           if (abs(rho) > thrwfn) then
              rho = rho * fac
              do igrid = ng0, ng1
                 vcoul(igrid, ifun) = vcoul(igrid, ifun) + rho * v2(igrid, jgrid)
              end do
           end if
        end do
     end do
  end do
  !$omp end parallel

end subroutine hf_print_vlocal_coul
!######################################################################
