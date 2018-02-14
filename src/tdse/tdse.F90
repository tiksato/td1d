!################################################################################
subroutine tdse

  use root_mod, only : name, icomp, iprint
  use wfn_mod, only : wfn, imethod, docip2

  implicit none

  real(kind(0d0)), external :: util_clock
  character(len = 256) :: fname
  integer :: itime
  real(kind(0d0)) :: rtime

  rtime = util_clock(.true., itime)

  ! hack for split propagator
  docip2 = .false.

  ! real-time simulation
  icomp = 1
  call tdse_cycle(wfn, imethod)

  ! print
  call tdse_summary(wfn, imethod)

  ! write
  fname = trim(name)//".wfn1"
  call wfn_write(fname, wfn, imethod)

  rtime = util_clock(.false., itime)
  write(6, "('# wall clock time for tdse:   ', E12.5)") rtime

end subroutine tdse
!################################################################################
