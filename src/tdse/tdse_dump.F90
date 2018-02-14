!################################################################################
subroutine tdse_dump(istep, time, lfield, next, wfn, imethod)

  use io_mod, only : iow
  use root_mod, only : name
  use prop_mod, only : cyctot
  use field_mod, only : period

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: istep, imethod
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(in) :: lfield
  integer, intent(inout) :: next
  complex(kind(0d0)), intent(in) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: optcyc
  character(len = 16) :: tag
  character(len = 256) :: fname

  optcyc = dble(time) / period

!debug
!  write(6, "('WRITE wavefunction at each iteration!')")
!  fname = trim(name)//".wfn"
!  call wfn_write(fname, wfn, imethod)
!  open(unit = iow, file = trim(fname), status = 'old', access = 'append')
!  write(iow,"(' #  info: ', i10, 3E25.15)") istep, dble(time), optcyc, lfield
!  close(unit = iow)
!debug

  if (optcyc > dble(next)) then

     if (imethod /= -1) then
        call util_itochar(next, cyctot, tag)
        fname = trim(name)//"."//trim(tag)//".wfn"
        call wfn_write(fname, wfn, imethod)
        
        write(6,"(' #     ', i10, 3E20.10, ': wfn is written to ', a)") &
             & istep, dble(time), optcyc, lfield, trim(fname)
     else
        write(6,"(' #     ', i10, 3E20.10, ': wfn is NOT written!')")
     end if

     next = next + 1

  end if

end subroutine tdse_dump
!################################################################################
