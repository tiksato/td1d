!################################################################################
subroutine init

  use root_mod, only : name, icomp, iprint
  use wfn_mod, only : wfn, imethod
  use const_mod, only : czero
  use io_mod, only : iow

  implicit none

  real(kind(0d0)), external :: util_clock
  character(len = 256) :: fname
  integer :: itime
  real(kind(0d0)) :: rtime

  rtime = util_clock(.true., itime)

  ! stationary state optimization
  icomp = 0
  call init_cycle(wfn, imethod)
! call init_cycle_old(wfn, imethod)

  ! print
  call init_summary(wfn, imethod)

  ! write
  fname = trim(name)//".wfn0"
  call wfn_write(fname, wfn, imethod)
  fname = trim(name)//".cic"
  call wfn_write_cic(fname, wfn, imethod)
  fname = trim(name)//".orb"
  call wfn_write_orb(fname, wfn, imethod)
!  fname = trim(name)//".den1"
!  call wfn_write_den1(fname, wfn, imethod)
!  fname = trim(name)//".den2"
!  call wfn_write_den2(fname, wfn, imethod)

!disabled  fname = trim(name)//".gvb0"
!disabled  open(iow, file = trim(fname), status = 'unknown', access = 'append')
!disabled  call wfn_print_gvb(iow, 0, czero, wfn, imethod)
!disabled  close(iow)

  rtime = util_clock(.false., itime)
  write(6, "('# wall clock time for init:   ', E12.5)") rtime

end subroutine init
!################################################################################
