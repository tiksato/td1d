!################################################################################
subroutine prop_exp(time, dt, len, wfn0, wfn, hwfn, imethod)

  use const_mod, only : zero, one, two, three, six, half, iunit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(in) :: wfn0(1:len)
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: lfield, ene, dtime, dinit
  complex(kind(0d0)) :: ttime, dt1, dt2, dt3, dt6
  real(kind(0d0)), external :: field
  real(kind(0d0)), external :: wfn_hprod, wfn_hprod2
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
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt6
  wfntmp(1:len) = wfnini(1:len) + hwfn(1:len) * dt2

! second step
!debug
!  write(6, "('prop_rk4-2: wfntmp')")
!  call tchf_print(6, wfntmp)
!debug
  ttime = time + dt * half
  lfield = field(ttime)
  ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt3
  wfntmp(1:len) = wfnini(1:len) + hwfn(1:len) * dt2
!debug  call wfn_print_cic(6, wfn, imethod)
!debug
!  write(6, "('prop_rk4-2: hwfn')")
!  call tchf_print(6, hwfn)
!debug

! third step
!debug
!  write(6, "('prop_rk4-3: wfntmp')")
!  call tchf_print(6, wfntmp)
!debug
  ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt3
  wfntmp(1:len) = wfnini(1:len) + hwfn(1:len) * dt1
!debug  call wfn_print_cic(6, wfn, imethod)
!debug
!  write(6, "('prop_rk4-3: hwfn')")
!  call tchf_print(6, hwfn)
!debug

! fourth step
!debug
!  write(6, "('prop_rk4-4: wfntmp')")
!  call tchf_print(6, wfntmp)
!debug
  ttime = time + dt
  lfield = field(ttime)
  ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
  wfn(1:len) = wfn(1:len) + hwfn(1:len) * dt6
!debug  call wfn_print_cic(6, wfn, imethod)
!debug
!  write(6, "('prop_rk4-4: hwfn')")
!  call tchf_print(6, hwfn)
!debug

  deallocate(wfntmp)
  deallocate(wfnini)

end subroutine prop_exp
!################################################################################
