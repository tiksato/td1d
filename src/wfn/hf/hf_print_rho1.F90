!######################################################################
subroutine hf_print_rho1(iout, istep, time, wfn)

  use root_mod, only : iprint
  use field_mod, only : period
  use thresh_mod, only : thrwfn
  use const_mod, only : zero, ctwo
  use grid_mod, only : ngrid, dgrid, x
  use wfn_mod, only : nfun, nocc
  use mol_mod, only : ne

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: iout, istep
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:*)
  !--------------------------------------------------------------------
  integer :: igrid, ifun
  real(kind(0d0)) :: optcyc, nel
  real(kind(0d0)), allocatable :: rho1(:)

  allocate(rho1(1:nfun))

  optcyc = dble(time) / period
  write(iout, "('# rho1:   ', i10, E15.6)") istep, optcyc

  nel = zero
  do igrid = 0, ngrid

     rho1(1:nfun) = zero
     do ifun = 1, nocc
        rho1(ifun) = dble(wfn(igrid, ifun) * conjg(wfn(igrid, ifun)) * ctwo)
        if (ifun <= ne(1)) then
           nel = nel + rho1(ifun) * dgrid
        end if
     end do

     write(iout, "(E15.6)", advance = 'no') x(igrid)
     do ifun = 1, nfun
        if (abs(rho1(ifun)) < thrwfn) rho1(ifun) = zero
        write(iout, "(E15.6)", advance = 'no') rho1(ifun)
     end do
     write(iout, *)
  end do

  if (iprint > 1) then
     write(iout, "('# nel:    ', E15.6)") nel
  end if

  deallocate(rho1)
  
end subroutine hf_print_rho1
