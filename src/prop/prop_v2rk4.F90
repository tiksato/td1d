!################################################################################
subroutine prop_v2rk4(dt, len, wfn, hwfn, imethod)

  use const_mod, only : zero, one, two, three, six, half, iunit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: dt
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: ene, dtime, dinit
  complex(kind(0d0)) :: dt1, dt2, dt3, dt6
  complex(kind(0d0)), allocatable :: wfnini(:)
  complex(kind(0d0)), allocatable :: wfntmp(:)

  dtime = abs(dt)
  dt1 = one
  dt2 = one / two
  dt3 = one / three
  dt6 = one / six
  allocate(wfnini(len))
  allocate(wfntmp(len))
  wfnini(1:len) = wfn(1:len)

! first step
  call wfn_v2prod(dtime, wfn, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt6
  wfntmp(1:len) = wfnini(1:len) + hwfn(1:len) * dt2

! second step
  call wfn_v2prod(dtime, wfntmp, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt3
  wfntmp(1:len) = wfnini(1:len) + hwfn(1:len) * dt2

! third step
  call wfn_v2prod(dtime, wfntmp, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt3
  wfntmp(1:len) = wfnini(1:len) + hwfn(1:len) * dt1

! fourth step
  call wfn_v2prod(dtime, wfntmp, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt6

  deallocate(wfntmp)
  deallocate(wfnini)

end subroutine prop_v2rk4
!################################################################################
