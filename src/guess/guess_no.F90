!################################################################################
subroutine guess_no(wfn, imethod)

  use root_mod, only : name

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(*)
  !--------------------------------------------------------------------
  character(len = 256) :: fname

  ! read from disk and transform to natural orbital representation
  fname = trim(name)//".wfnin"
  call wfn_read(fname, wfn, imethod)
  call wfn_transno(wfn, imethod)

end subroutine guess_no
!################################################################################
