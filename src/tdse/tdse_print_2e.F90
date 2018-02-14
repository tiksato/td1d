!################################################################################
subroutine tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)

  use root_mod, only : name, nprint_rho1
  use const_mod, only : one, four
  use grid_mod, only : rmax, x0, xmask
  use field_mod, only : period
  use io_mod, only : iostdo, iow

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: istep, imethod
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(in) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)), intent(in) :: wfn0(*)
  complex(kind(0d0)), intent(in) :: wfn(*)
  !--------------------------------------------------------------------
  character(len = 256) :: fname
  logical :: do_rho1
  real(kind(0d0)), external :: wfn_op1e
  real(kind(0d0)), external :: wfn_norm
  complex(kind(0d0)), external :: wfn_iprod
  real(kind(0d0)) :: optcyc, p0, p1, p2, ptot, dip, vel, acc
  complex(kind(0d0)) :: c0

  do_rho1 = nprint_rho1 > 0 .and. mod(istep, nprint_rho1) == 0

  optcyc = real(time) / period
  c0 = wfn_iprod(x0, wfn0, wfn, imethod)
  dip = wfn_op1e(0, -one, wfn, imethod)
  vel = wfn_op1e(2, -one, wfn, imethod)
  acc = wfn_op1e(1, -one, wfn, imethod)
  ptot = wfn_norm(rmax, wfn, p0, p1, p2, imethod)

  write(6,"(' step  ', i10, 15E15.6)") istep, real(time), optcyc, lfield, &
       &  c0, p0, p1, p2, q0, q1, q2, dip, vel, acc, ene

!disabled  if (do_rho1) then
!disabled     fname = trim(name)//".gvb"
!disabled     open(iow, file = trim(fname), status = 'unknown', access = 'append')
!disabled     call wfn_print_gvb(iow, istep, time, wfn, imethod)
!disabled     close(iow)
!disabled
!disabled     fname = trim(name)//".vlocal"
!disabled     open(iow, file = trim(fname), status = 'unknown', access = 'append')
!disabled     call wfn_print_vlocal(iow, istep, time, lfield, wfn, imethod)
!disabled     close(iow)
!disabled  end if

end subroutine tdse_print_2e
!################################################################################
