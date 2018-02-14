!################################################################################
subroutine tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)

  use root_mod, only : name, nprint_ene, nprint_op1e, nprint_op1x, nprint_ionp, nprint_rho1
  use const_mod, only : zero, one, four
  use grid_mod, only : rmax, x0, xmask
  use field_mod, only : period
  use io_mod, only : iostdo, iow, io_rho1

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: istep, imethod
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(inout) :: lfield, ene
  complex(kind(0d0)), intent(in) :: wfn0(*)
  complex(kind(0d0)), intent(in) :: wfn(*)
  !--------------------------------------------------------------------
  logical :: do_ene, do_op1e, do_op1x, do_ionp, do_rho1
  real(kind(0d0)) :: optcyc
  character(len = 256) :: fname

  do_ene  = nprint_ene > 0 .and. mod(istep, nprint_ene) == 0
  do_op1e = nprint_op1e > 0 .and. mod(istep, nprint_op1e) == 0
  do_op1x = nprint_op1x > 0 .and. mod(istep, nprint_op1x) == 0
  do_ionp = nprint_ionp > 0 .and. mod(istep, nprint_ionp) == 0
  do_rho1 = nprint_rho1 > 0 .and. mod(istep, nprint_rho1) == 0

  optcyc = dble(time) / period  
  if (do_ene) write(iostdo, "(' ene:    ', i10, 3E20.10)") istep, optcyc, lfield, ene
  if (do_op1e) call wfn_print_op1e(iostdo, istep, time, x0, wfn, imethod)
  if (do_op1x) call wfn_print_op1x(iostdo, istep, time, wfn, imethod)
  if (do_ionp) call wfn_print_ionp(iostdo, istep, time, rmax, wfn0, wfn, imethod)
!disabled  if (do_rho1) then
!disabled     fname = trim(name)//".rho1"
!disabled     open(io_rho1, file = trim(fname), status = 'old', access = 'append')
!disabled     call wfn_print_rho1(io_rho1, istep, time, wfn, imethod)
!disabled     close(io_rho1)
!disabled
!disabled     fname = trim(name)//".gvb"
!disabled     open(iow, file = trim(fname), status = 'unknown', access = 'append')
!disabled     call wfn_print_gvb(iow, istep, time, wfn, imethod)
!disabled     close(iow)
!disabled
!disabled     fname = trim(name)//".vlocal"
!disabled     open(iow, file = trim(fname), status = 'unknown', access = 'append')
!disabled     call wfn_print_vlocal(iow, istep, time, lfield, wfn, imethod)
!disabled     close(iow)
!disabled
!disabled!     fname = trim(name)//".gvb"
!disabled!     open(iow, file = trim(fname), status = 'unknown', access = 'append')
!disabled!     call wfn_print_gvb(iow, istep, time, wfn, imethod)
!disabled!     close(iow)
!disabled  end if

end subroutine tdse_print
!################################################################################
