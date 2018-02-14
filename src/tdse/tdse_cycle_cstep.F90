!################################################################################
subroutine tdse_cycle_cstep(wfn, imethod)

  use root_mod, only : name
  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use mol_mod, only : ne
  use field_mod, only : period
  use prop_mod, only : totstep, dstep, initime, prop_type
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & noci, nomo, nolag1, nolag2, nolag3, normci, sep_fc

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  complex(kind(0d0)), allocatable :: wfn0(:), hwfn(:), wfnp(:)

  integer :: len, istep, next
  real(kind(0d0)) :: dtime, lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, time
  character(len = 256) :: fname

  if (prop_type == 'ABM' .or. prop_type == 'VABM') then
     call tdse_cycle_cstep_abm(wfn, imethod)
     return
  else if (prop_type == 'SPLIT') then
!    call tdse_cycle_cstep_split1(wfn, imethod)
     call tdse_cycle_cstep_split2(wfn, imethod)
!    call tdse_cycle_cstep_split_debug(wfn, imethod)
     return
  end if

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(1:len))
  allocate(wfnp(1:len))
  wfn0(1:len) = wfn(1:len)
  hwfn(1:len) = czero
  wfnp(1:len) = czero

  ! overwrite frozen-core orbitals 
  if (sep_fc .ne. 0) then
     fname = trim(name)//".ofc"
     call wfn_read_fcore(fname, wfn, imethod)
  end if

  doci   = .not. noci
  domo   = .not. nomo
  domop1 = .not. nolag1
  domop2 = .not. nolag2
  domoq  = .not. nolag3
  domop  = domop1 .or. domop2

  dt = dstep * runit
  q0 = zero
  q1 = zero
  q2 = zero

  time = initime
  next = int(dble(initime/period)) + 1

  do istep = 0, totstep

!DEBUG
!if (istep > 1000) stop "FOR DEBUG @ tdse_cycle_step."
!DEBUG

     dtime = abs(dt)
     lfield = field(time)

     ! this is not necessary for split operator propagators...
     ene = wfn_hprod(.true., lfield, dtime, wfn0, wfn, hwfn, imethod)

     ! output
     if (imethod == -1 .or. imethod == 4 .or. imethod == -4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        !if (ne(3) == 2) then
        !   call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        !end if
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     ! propagation from t to t + dt
     call prop_new(time, dt, len, wfn0, wfn, hwfn, imethod, prop_type)

     if (normci) then
        call wfn_ort_ci(wfn, imethod)
     end if

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

     time = time + dt
!DEBUG
!     stop 'for debug @ tdse_cycle_cstep'
!DEBUG
  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_cstep
!################################################################################
!################################################################################
subroutine tdse_cycle_cstep_abm(wfn, imethod)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : totstep, dstep, initime, prop_type
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  complex(kind(0d0)), allocatable :: wfn0(:), hwfn(:,:), twfn(:), wfnp(:)

  integer :: len, istep, next
  real(kind(0d0)) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(twfn(1:len))
  allocate(hwfn(1:len, 0:3))
  allocate(wfnp(1:len))
  call util_zcopy(len, wfn,   1, wfn0, 1)
  call util_zcopy(len, czero, 0, hwfn, 1)
  call util_zcopy(len, czero, 0, wfnp, 1)

  doci   = .not. noci
  domo   = .not. nomo
  domop1 = .not. nolag1
  domop2 = .not. nolag2
  domoq  = .not. nolag3
  domop  = domop1 .or. domop2

  dt = dstep * runit
  q0 = zero
  q1 = zero
  q2 = zero

  time = initime
  next = int(dble(initime/period)) + 1

  ! generate first three points by RK4 method
  do istep = 0, 2

     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, twfn, imethod)
     hwfn(1:len, 3-istep) = twfn(1:len)

     ! output
     if (imethod <= 1 .or. imethod == 4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     ! propagation from t to t + dt
     call prop(time, dt, len, wfn0, wfn, twfn, imethod, 'RK4')

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

     time = time + dt

  end do


  do istep = 3, totstep

     lfield = field(time)

     ! this is not necessary for split operator propagators...
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, twfn, imethod)
     hwfn(1:len, 0) = twfn(1:len)

     ! output
     if (imethod <= 1 .or. imethod == 4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     ! propagation from t to t + dt
     call prop_abm4(prop_type, time, dt, len, wfn0, wfn, hwfn, imethod)

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

     time = time + dt

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(twfn)
  deallocate(wfn0)

end subroutine tdse_cycle_cstep_abm
!################################################################################
!################################################################################
subroutine tdse_cycle_cstep_split1(wfn, imethod)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : totstep, dstep, initime, prop_type
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod, wfn_prop_p2
  complex(kind(0d0)), allocatable :: wfn0(:), hwfn(:), wfnp(:)

  integer :: len, istep, next
  real(kind(0d0)) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, dt2, time

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(1:len))
  allocate(wfnp(1:len))
  call util_zcopy(len, wfn,   1, wfn0, 1)
  call util_zcopy(len, czero, 0, hwfn, 1)
  call util_zcopy(len, czero, 0, wfnp, 1)

  dt = dstep * runit
  dt2 = dt / 2.d+0
  q0 = zero
  q1 = zero
  q2 = zero

  next = int(dble(initime/period)) + 1

  do istep = 0, totstep

     ! propagation from t to t + dt / 2
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .not. nolag1
     domop2 = .false.
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2

     time = initime + istep * dt
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)

     ! output
     if (imethod <= 1 .or. imethod == 4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     call prop(time, dt2, len, wfn0, wfn, hwfn, imethod, prop_type)

     ! propagation from t to t + dt 
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .false.
     domop2 = .not. nolag2
     domoq  = .false.
     domop  = domop1 .or. domop2

     time = initime + istep * dt + dt2
     lfield = field(time)
     ene = wfn_prop_p2(lfield, dt, wfn, hwfn, imethod)

     ! propagation from t + dt /2 to t + dt
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .not. nolag1
     domop2 = .false.
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2

     time = initime + istep * dt + dt2
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)
     call prop(time, dt2, len, wfn0, wfn, hwfn, imethod, prop_type)

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_cstep_split1
!################################################################################
!################################################################################
subroutine tdse_cycle_cstep_split2(wfn, imethod)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : totstep, dstep, initime, prop_type
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod, wfn_prop_p2
  complex(kind(0d0)), allocatable :: wfn0(:), hwfn(:), wfnp(:)

  integer :: len, istep, next
  real(kind(0d0)) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, dt2, time

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(1:len))
  allocate(wfnp(1:len))
  call util_zcopy(len, wfn,   1, wfn0, 1)
  call util_zcopy(len, czero, 0, hwfn, 1)
  call util_zcopy(len, czero, 0, wfnp, 1)

  dt = dstep * runit
  dt2 = dt / 2.d+0
  q0 = zero
  q1 = zero
  q2 = zero

  next = int(dble(initime/period)) + 1

  do istep = 0, totstep

     ! propagation from t to t + dt / 2
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2

     time = initime + istep * dt
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)

     ! output
     if (imethod <= 1 .or. imethod == 4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     call prop(time, dt2, len, wfn0, wfn, hwfn, imethod, prop_type)

     ! propagation from t to t + dt 
     doci   = .false.
     domo   = .not. nomo
     domop1 = .false.
     domop2 = .not. nolag2
     domoq  = .false.
     domop  = domop1 .or. domop2

     time = initime + istep * dt + dt2
     lfield = field(time)
     ene = wfn_prop_p2(lfield, dt, wfn, hwfn, imethod)

     ! propagation from t + dt /2 to t + dt
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2

     time = initime + istep * dt + dt2
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)
     call prop(time, dt2, len, wfn0, wfn, hwfn, imethod, prop_type)

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_cstep_split2
!################################################################################
!################################################################################
subroutine tdse_cycle_cstep_split_debug(wfn, imethod)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : totstep, dstep, initime, prop_type
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, docip2, &
       & noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod, wfn_prop_p2
  complex(kind(0d0)), allocatable :: wfn0(:), hwfn(:), wfnp(:)

  integer :: len, istep, next
  real(kind(0d0)) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, dt2, time

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(1:len))
  allocate(wfnp(1:len))
  call util_zcopy(len, wfn,   1, wfn0, 1)
  call util_zcopy(len, czero, 0, hwfn, 1)
  call util_zcopy(len, czero, 0, wfnp, 1)

  dt = dstep * runit
  dt2 = dt / 2.d+0
  q0 = zero
  q1 = zero
  q2 = zero

  next = int(dble(initime/period)) + 1

  do istep = 0, totstep

     ! propagation from t to t + dt / 2
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .not. nolag1
     domop2 = .false.
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2
     docip2 = .false.

     time = initime + istep * dt
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)

     ! output
     if (imethod <= 1 .or. imethod == 4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     call prop(time, dt2, len, wfn0, wfn, hwfn, imethod, prop_type)

     ! propagation from t to t + dt 
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .false.
     domop2 = .not. nolag2
     domoq  = .false.
     domop  = domop1 .or. domop2
     docip2 = .true.

     time = initime + istep * dt
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn0, wfn, hwfn, imethod, prop_type)

     ! propagation from t + dt /2 to t + dt
     doci   = .not. noci
     domo   = .not. nomo
     domop1 = .not. nolag1
     domop2 = .false.
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2
     docip2 = .false.

     time = initime + istep * dt + dt2
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)
     call prop(time, dt2, len, wfn0, wfn, hwfn, imethod, prop_type)

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_cstep_split_debug
!################################################################################
