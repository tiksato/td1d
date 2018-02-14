!################################################################################
subroutine prop_vrk5_new(time, dt, len, wfn0, wfn, hwfn, imethod)

  use const_mod, only : zero, one, two, three, six, half, iunit, czero
  use prop_mod, only : rk5_a, rk5_b, rk5_crk5, rk5_crk4, &
       & prop_tol, prop_safety, prop_ulscal, dstep_ll, dstep_ul

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(inout) :: dt
  complex(kind(0d0)), intent(in) :: wfn0(1:len)
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len)
  !--------------------------------------------------------------------

  logical :: acc_enough
  integer :: itry
  real(kind(0d0)) :: lfield, ene, scale, dtime, dinit
  real(kind(0d0)) :: n2err, msdev, sqerr
  complex(kind(0d0)) :: ttime, rk5_cerr(1:6)
  real(kind(0d0)), external :: field
  real(kind(0d0)), external :: wfn_hprod
  real(kind(0d0)), external :: wfn_n2err
  complex(kind(0d0)), external :: util_zdotc
  complex(kind(0d0)), allocatable :: wfn_rk5(:), wfn_err(:)
  complex(kind(0d0)), allocatable :: wfnini(:), hwfnini(:), wfntmp(:)
  complex(kind(0d0)), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)

  allocate(wfnini(len))
  allocate(hwfnini(len))
  allocate(k1(len))
  allocate(k2(len))
  allocate(k3(len))
  allocate(k4(len))
  allocate(k5(len))
  allocate(k6(len))
  allocate(wfntmp(len))
  allocate(wfn_rk5(len))
  allocate(wfn_err(len))

  itry = 0
  dinit = abs(dt)
  acc_enough = .false.
  wfnini(1:len) = wfn(1:len)
  hwfnini(1:len) = hwfn(1:len)
  rk5_cerr(1:6) = rk5_crk5(1:6) - rk5_crk4(1:6)

  do while (.not. acc_enough)

     itry = itry + 1
     dtime = abs(dt)
     scale = dtime / dinit

     wfn_rk5(1:len) = wfnini(1:len)
     wfn_err(1:len) = czero
     k1(1:len) = hwfnini(1:len) * scale

     ! step 1
     wfn_rk5(1:len) = wfn_rk5(1:len) + k1(1:len) * rk5_crk5(1)
     wfn_err(1:len) = wfn_err(1:len) + k1(1:len) * rk5_cerr(1)

     ! step 2
     ttime = time + dt * rk5_a(2)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 2)
     ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
     k2(1:len) = hwfn(1:len)
     wfn_rk5(1:len) = wfn_rk5(1:len) + k2(1:len) * rk5_crk5(2)
     wfn_err(1:len) = wfn_err(1:len) + k2(1:len) * rk5_cerr(2)

     ! step 3
     ttime = time + dt * rk5_a(3)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 3) &
                                 & + k2(1:len) * rk5_b(2, 3)
     ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
     k3(1:len) = hwfn(1:len)
     wfn_rk5(1:len) = wfn_rk5(1:len) + k3(1:len) * rk5_crk5(3)
     wfn_err(1:len) = wfn_err(1:len) + k3(1:len) * rk5_cerr(3)

     ! step 4
     ttime = time + dt * rk5_a(4)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 4) &
                                 & + k2(1:len) * rk5_b(2, 4) &
                                 & + k3(1:len) * rk5_b(3, 4)
     ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
     k4(1:len) = hwfn(1:len)
     wfn_rk5(1:len) = wfn_rk5(1:len) + k4(1:len) * rk5_crk5(4)
     wfn_err(1:len) = wfn_err(1:len) + k4(1:len) * rk5_cerr(4)

     ! step 5
     ttime = time + dt * rk5_a(5)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 5) &
                                 & + k2(1:len) * rk5_b(2, 5) &
                                 & + k3(1:len) * rk5_b(3, 5) &
                                 & + k4(1:len) * rk5_b(4, 5)
     ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
     k5(1:len) = hwfn(1:len)
     wfn_rk5(1:len) = wfn_rk5(1:len) + k5(1:len) * rk5_crk5(5)
     wfn_err(1:len) = wfn_err(1:len) + k5(1:len) * rk5_cerr(5)

     ! step 6
     ttime = time + dt * rk5_a(6)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 6) &
                                 & + k2(1:len) * rk5_b(2, 6) &
                                 & + k3(1:len) * rk5_b(3, 6) &
                                 & + k4(1:len) * rk5_b(4, 6) &
                                 & + k5(1:len) * rk5_b(5, 6)
     ene = wfn_hprod(.false., lfield, dtime, wfn0, wfntmp, hwfn, imethod)
     k6(1:len) = hwfn(1:len)
     wfn_rk5(1:len) = wfn_rk5(1:len) + k6(1:len) * rk5_crk5(6)
     wfn_err(1:len) = wfn_err(1:len) + k6(1:len) * rk5_cerr(6)

! try and error
     n2err = wfn_n2err(wfn_rk5, wfn_err, imethod)
     msdev = dble(util_zdotc(len, wfn_err, 1, wfn_err, 1)) / len
     sqerr = n2err
!     sqerr = msdev
! try and error

     if (abs(dt) < dstep_ll) then
        write(6, "('WARNING: exit VRK5 for too small dt.')")
        acc_enough = .true.
        if (abs(sqerr) < prop_tol) then
           scale = prop_safety * (prop_tol / abs(sqerr)) ** (one / 10.D+0)
           scale = min(prop_ulscal, max(one, scale))
        else
           scale = one
        end if
     else
        if (abs(sqerr) < prop_tol) then
           acc_enough = .true.
           scale = prop_safety * (prop_tol / abs(sqerr)) ** (one / 10.D+0)
           scale = min(prop_ulscal, max(one, scale))
        else
           acc_enough = .false.
           scale = prop_safety * (prop_tol / abs(sqerr)) ** (one / 8.D+0)
        end if
     end if
     !debug
     write(6, "('VRK5: ',I10, 4E20.5, 9X, l1)") itry, prop_tol, &
          & abs(n2err), abs(msdev), real(dt), acc_enough
     !debug

     dt = dt * scale
     if (abs(dt) > dstep_ul) dt = dt / abs(dt) * dstep_ul
  end do

  wfn(1:len) = wfn_rk5(1:len)  

  deallocate(wfn_err)
  deallocate(wfn_rk5)
  deallocate(wfntmp)
  deallocate(k6)
  deallocate(k5)
  deallocate(k4)
  deallocate(k3)
  deallocate(k2)
  deallocate(k1)
  deallocate(hwfnini)
  deallocate(wfnini)

!debug
!write(6,"('rk5_a2 = ',f12.5)") rk5_a(2)
!write(6,"('rk5_a3 = ',f12.5)") rk5_a(3)
!write(6,"('rk5_a4 = ',f12.5)") rk5_a(4)
!write(6,"('rk5_a5 = ',f12.5)") rk5_a(5)
!write(6,"('rk5_a6 = ',f12.5)") rk5_a(6)
!write(6,"('debug: rk5_b(:,2) = ',5f12.5)") rk5_b(1:5, 2)
!write(6,"('debug: rk5_b(:,3) = ',5f12.5)") rk5_b(1:5, 3)
!write(6,"('debug: rk5_b(:,4) = ',5f12.5)") rk5_b(1:5, 4)
!write(6,"('debug: rk5_b(:,5) = ',5f12.5)") rk5_b(1:5, 5)
!write(6,"('debug: rk5_b(:,6) = ',5f12.5)") rk5_b(1:5, 6)
!write(6,"('rk5_c1 = ',f12.5)") rk5_crk5(1)
!write(6,"('rk5_c2 = ',f12.5)") rk5_crk5(2)
!write(6,"('rk5_c3 = ',f12.5)") rk5_crk5(3)
!write(6,"('rk5_c4 = ',f12.5)") rk5_crk5(4)
!write(6,"('rk5_c5 = ',f12.5)") rk5_crk5(5)
!write(6,"('rk5_c6 = ',f12.5)") rk5_crk5(6)
!debug
end subroutine prop_vrk5_new
!################################################################################
subroutine prop_vrk5(time, dt, len, wfn0, wfn, hwfn, imethod)

  use const_mod, only : zero, one, two, three, six, half, iunit, czero
  use prop_mod, only : rk5_a, rk5_b, rk5_crk5, rk5_crk4, &
       & prop_tol, prop_safety, prop_ulscal, dstep_ll, dstep_ul

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(inout) :: dt
  complex(kind(0d0)), intent(in) :: wfn0(1:len)
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len)
  !--------------------------------------------------------------------

  logical :: acc_enough
  integer :: itry
  real(kind(0d0)) :: lfield, ene, scale
  real(kind(0d0)) :: n2err, msdev, sqerr
  complex(kind(0d0)) :: ttime, idt, rk5_cerr(1:6)
  real(kind(0d0)), external :: field
  real(kind(0d0)), external :: wfn_hprod
  real(kind(0d0)), external :: wfn_n2err
  complex(kind(0d0)), external :: util_zdotc
  complex(kind(0d0)), allocatable :: wfn_rk5(:), wfn_err(:)
  complex(kind(0d0)), allocatable :: wfnini(:), hwfnini(:), wfntmp(:)
  complex(kind(0d0)), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)

  allocate(wfnini(len))
  allocate(hwfnini(len))
  allocate(k1(len))
  allocate(k2(len))
  allocate(k3(len))
  allocate(k4(len))
  allocate(k5(len))
  allocate(k6(len))
  allocate(wfntmp(len))
  allocate(wfn_rk5(len))
  allocate(wfn_err(len))

  itry = 0
  acc_enough = .false.
  wfnini(1:len) = wfn(1:len)
  hwfnini(1:len) = hwfn(1:len)
  rk5_cerr(1:6) = rk5_crk5(1:6) - rk5_crk4(1:6)

  do while (.not. acc_enough)

     itry = itry + 1
     idt = - iunit * dt

     wfn_rk5(1:len) = wfnini(1:len)
     wfn_err(1:len) = czero
     k1(1:len) = hwfnini(1:len) * idt

     ! step 1
     wfn_rk5(1:len) = wfn_rk5(1:len) + k1(1:len) * rk5_crk5(1)
     wfn_err(1:len) = wfn_err(1:len) + k1(1:len) * rk5_cerr(1)

     ! step 2
     ttime = time + dt * rk5_a(2)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 2)
     ene = wfn_hprod(.false., lfield, dt, wfn0, wfntmp, hwfn, imethod)
     k2(1:len) = hwfn(1:len) * idt
     wfn_rk5(1:len) = wfn_rk5(1:len) + k2(1:len) * rk5_crk5(2)
     wfn_err(1:len) = wfn_err(1:len) + k2(1:len) * rk5_cerr(2)

     ! step 3
     ttime = time + dt * rk5_a(3)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 3) &
                                 & + k2(1:len) * rk5_b(2, 3)
     ene = wfn_hprod(.false., lfield, dt, wfn0, wfntmp, hwfn, imethod)
     k3(1:len) = hwfn(1:len) * idt
     wfn_rk5(1:len) = wfn_rk5(1:len) + k3(1:len) * rk5_crk5(3)
     wfn_err(1:len) = wfn_err(1:len) + k3(1:len) * rk5_cerr(3)

     ! step 4
     ttime = time + dt * rk5_a(4)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 4) &
                                 & + k2(1:len) * rk5_b(2, 4) &
                                 & + k3(1:len) * rk5_b(3, 4)
     ene = wfn_hprod(.false., lfield, dt, wfn0, wfntmp, hwfn, imethod)
     k4(1:len) = hwfn(1:len) * idt
     wfn_rk5(1:len) = wfn_rk5(1:len) + k4(1:len) * rk5_crk5(4)
     wfn_err(1:len) = wfn_err(1:len) + k4(1:len) * rk5_cerr(4)

     ! step 5
     ttime = time + dt * rk5_a(5)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 5) &
                                 & + k2(1:len) * rk5_b(2, 5) &
                                 & + k3(1:len) * rk5_b(3, 5) &
                                 & + k4(1:len) * rk5_b(4, 5)
     ene = wfn_hprod(.false., lfield, dt, wfn0, wfntmp, hwfn, imethod)
     k5(1:len) = hwfn(1:len) * idt
     wfn_rk5(1:len) = wfn_rk5(1:len) + k5(1:len) * rk5_crk5(5)
     wfn_err(1:len) = wfn_err(1:len) + k5(1:len) * rk5_cerr(5)

     ! step 6
     ttime = time + dt * rk5_a(6)
     lfield = field(ttime)
     wfntmp(1:len) = wfnini(1:len) + k1(1:len) * rk5_b(1, 6) &
                                 & + k2(1:len) * rk5_b(2, 6) &
                                 & + k3(1:len) * rk5_b(3, 6) &
                                 & + k4(1:len) * rk5_b(4, 6) &
                                 & + k5(1:len) * rk5_b(5, 6)
     ene = wfn_hprod(.false., lfield, dt, wfn0, wfntmp, hwfn, imethod)
     k6(1:len) = hwfn(1:len) * idt
     wfn_rk5(1:len) = wfn_rk5(1:len) + k6(1:len) * rk5_crk5(6)
     wfn_err(1:len) = wfn_err(1:len) + k6(1:len) * rk5_cerr(6)

! try and error
     n2err = wfn_n2err(wfn_rk5, wfn_err, imethod)
     msdev = dble(util_zdotc(len, wfn_err, 1, wfn_err, 1)) / len
     sqerr = n2err
!     sqerr = msdev
! try and error

     if (abs(dt) < dstep_ll) then
        write(6, "('WARNING: exit VRK5 for too small dt.')")
        acc_enough = .true.
        if (abs(sqerr) < prop_tol) then
           scale = prop_safety * (prop_tol / abs(sqerr)) ** (one / 10.D+0)
           scale = min(prop_ulscal, max(one, scale))
        else
           scale = one
        end if
     else
        if (abs(sqerr) < prop_tol) then
           acc_enough = .true.
           scale = prop_safety * (prop_tol / abs(sqerr)) ** (one / 10.D+0)
           scale = min(prop_ulscal, max(one, scale))
        else
           acc_enough = .false.
           scale = prop_safety * (prop_tol / abs(sqerr)) ** (one / 8.D+0)
        end if
     end if
     !debug
     write(6, "('VRK5: ',I10, 4E20.5, 9X, l1)") itry, prop_tol, &
          & abs(n2err), abs(msdev), real(dt), acc_enough
     !debug

     dt = dt * scale
     if (abs(dt) > dstep_ul) dt = dt / abs(dt) * dstep_ul
  end do

  wfn(1:len) = wfn_rk5(1:len)  

  deallocate(wfn_err)
  deallocate(wfn_rk5)
  deallocate(wfntmp)
  deallocate(k6)
  deallocate(k5)
  deallocate(k4)
  deallocate(k3)
  deallocate(k2)
  deallocate(k1)
  deallocate(hwfnini)
  deallocate(wfnini)

!debug
!write(6,"('rk5_a2 = ',f12.5)") rk5_a(2)
!write(6,"('rk5_a3 = ',f12.5)") rk5_a(3)
!write(6,"('rk5_a4 = ',f12.5)") rk5_a(4)
!write(6,"('rk5_a5 = ',f12.5)") rk5_a(5)
!write(6,"('rk5_a6 = ',f12.5)") rk5_a(6)
!write(6,"('debug: rk5_b(:,2) = ',5f12.5)") rk5_b(1:5, 2)
!write(6,"('debug: rk5_b(:,3) = ',5f12.5)") rk5_b(1:5, 3)
!write(6,"('debug: rk5_b(:,4) = ',5f12.5)") rk5_b(1:5, 4)
!write(6,"('debug: rk5_b(:,5) = ',5f12.5)") rk5_b(1:5, 5)
!write(6,"('debug: rk5_b(:,6) = ',5f12.5)") rk5_b(1:5, 6)
!write(6,"('rk5_c1 = ',f12.5)") rk5_crk5(1)
!write(6,"('rk5_c2 = ',f12.5)") rk5_crk5(2)
!write(6,"('rk5_c3 = ',f12.5)") rk5_crk5(3)
!write(6,"('rk5_c4 = ',f12.5)") rk5_crk5(4)
!write(6,"('rk5_c5 = ',f12.5)") rk5_crk5(5)
!write(6,"('rk5_c6 = ',f12.5)") rk5_crk5(6)
!debug
end subroutine prop_vrk5
!################################################################################
