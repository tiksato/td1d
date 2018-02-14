!################################################################################
subroutine prop_euler_new(len, wfn, hwfn)

  use const_mod, only : iunit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len)
  !--------------------------------------------------------------------

  wfn(1:len) = wfn(1:len) + hwfn(1:len)

end subroutine prop_euler_new
!################################################################################
subroutine prop_euler(dt, len, wfn, hwfn)

  use const_mod, only : iunit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len
  complex(kind(0d0)), intent(in) :: dt
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len)
  !--------------------------------------------------------------------
  complex(kind(0d0)) :: dt1

  dt1 = - iunit * dt
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt1

end subroutine prop_euler
!################################################################################
