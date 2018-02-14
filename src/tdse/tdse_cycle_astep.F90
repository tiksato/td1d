!################################################################################
subroutine tdse_cycle_astep(wfn, imethod)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : totstep, cyctot, dstep, initime, prop_type
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
  complex(kind(0d0)), allocatable :: wfn0(:), hwfn(:), wfnp(:)

  integer :: len, istep, next
  real(kind(0d0)) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt0, dt, time

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(1:len))
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

  dt0 = dstep * runit
  dt = dt0
  q0 = zero
  q1 = zero
  q2 = zero

  istep = 0
  time = initime
  next = int(dble(initime/period)) + 1

  do

     lfield = field(time)

     ene = wfn_hprod(.true., lfield, dt, wfn0, wfn, hwfn, imethod)
     call wfn_adjdt(dt, dt0, wfn, hwfn, imethod)

     ! output
     if (imethod <= 1 .or. imethod == 4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     ! propagation from t to t + dt
     call prop(time, dt, len, wfn0, wfn, hwfn, imethod, prop_type)

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethod)
     end if

     if (dble(time) > cyctot*period) exit
     istep = istep + 1
     time = time + dt

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_astep
!################################################################################
