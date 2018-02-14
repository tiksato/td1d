!######################################################################
subroutine hf_print_gvb(iout, istep, time, wfn)

  use wfn_mod, only : ehf, nfun, nspin
  use root_mod, only : iprint
  use field_mod, only : period
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, x
  use const_mod, only : zero, one, half, runit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: iout, istep
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:*)
  !--------------------------------------------------------------------
  integer :: igrid, ifun
  real(kind(0d0)) :: optcyc

  optcyc = dble(time) / period
  write(iout, "('# gvb:        ', i15, E15.6)") istep, optcyc

  do igrid = 0, ngrid
     write(iout, "(E15.6)", advance = 'no') x(igrid)
     do ifun = 1, nfun * nspin
        write(iout, "(2E15.6)", advance = 'no') wfn(igrid, ifun)
     end do
     write(iout, *)
  end do
  
  write(iout, *)
  
end subroutine hf_print_gvb
