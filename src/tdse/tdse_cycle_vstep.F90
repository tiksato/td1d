!################################################################################
subroutine tdse_cycle_vstep(wfn, imethod)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : totstep, cyctot, dstep, initime, prop_type
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & noci, nomo, nolag1, nolag2, nolag3, normci, dpsi_dyreg0, dpsi_smin

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
  real(kind(0d0)) :: optcyc, dtime, lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(len))
  allocate(wfnp(len))
  wfn0(1:len) = wfn(1:len)
  hwfn(1:len) = czero
  wfnp(1:len) = czero

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

  istep = 0
  time = initime
  next = int(dble(initime/period)) + 1

  do

     dtime = abs(dt)
     lfield = field(time)
     domop2 = .not. nolag2
     ene = wfn_hprod(.true., lfield, dtime, wfn0, wfn, hwfn, imethod)

!disabled     ! minimum singular value of aa-rotation matrix
!disabled     if (dpsi_smin < dpsi_dyreg0) then
!disabled        domop2 = .false.
!disabled        optcyc = dble(time) / period  
!disabled        write(6, "(' DoMOP2 = F: ', i10, 3E15.6)") istep, optcyc, dpsi_smin, dpsi_dyreg0
!disabled     end if

     ! output
     if (imethod == -1 .or. imethod == 4 .or. imethod == -4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
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

     if (dble(time) > cyctot*period) exit
     istep = istep + 1
     time = time + dt

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_vstep
!################################################################################
