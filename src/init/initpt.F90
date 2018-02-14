!################################################################################
subroutine initpt(wfn, imethod)

  use const_mod, only : czero
  use root_mod, only : name, icomp, iprint
  use init_mod, only : guess_type, prop_type_init
  use wfn_mod, only : dopt

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)

  real(kind(0d0)), external :: util_clock
  character(len = 256) :: fname
  integer :: itime
  real(kind(0d0)) :: rtime

  if (.not. dopt) return

  rtime = util_clock(.true., itime)

  if (trim(guess_type) == 'READ') then
     fname = trim(name)//".wfnptin"
     call wfn_readpt(fname, wfn, imethod)
  else
     call wfn_gsmp2(wfn, imethod)
  end if

  fname = trim(name)//".wfnpt0"
  call wfn_writept(fname, wfn, imethod)
  
  rtime = util_clock(.false., itime)
  write(6, "('# wall clock time for initpt: ', E12.5)") rtime

end subroutine initpt
!################################################################################
