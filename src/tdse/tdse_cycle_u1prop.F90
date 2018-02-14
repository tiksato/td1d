!################################################################################
subroutine tdse_cycle_u1prop(wfn, imethod)

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

  len = wfn_size(imethod)
  allocate(wfn0(1:len))
  allocate(hwfn(1:len))
  allocate(wfnp(1:len))
  call util_zcopy(len, wfn,   1, wfn0, 1)
  call util_zcopy(len, czero, 0, hwfn, 1)
  call util_zcopy(len, czero, 0, wfnp, 1)

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

  q0 = zero
  q1 = zero
  q2 = zero
  time = initime
  dt = dstep * runit
  next = int(dble(initime/period)) + 1

  do istep = 0, totstep

     dtime = abs(dt)
     lfield = field(time)

     ! this is not necessary for split operator propagators...
     ene = wfn_hprod(.true., lfield, dtime, wfn0, wfn, hwfn, imethod)

     ! output
     if (imethod <= 1 .or. imethod == 4 .or. imethod == -4) then
        call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     else
        if (ne(3) == 2) call tdse_print_2e(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, imethod)
        call tdse_print(istep, time, lfield, ene, wfn0, wfn, imethod)
        call tdse_dump(istep, time, lfield, next, wfn, imethod)
     end if

     ! propagation from t to t + dt
     call wfn_u1prop(time, dt, wfn, hwfn, imethod)

     if (normci) then
        call wfn_ort_ci(wfn, imethod)
     end if

     if (domask) then
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)
     else if (docap) then
        stop 'tdse_cycle_u1prop: wfn_capped nyi.'
     end if

     time = time + dt

  end do

  deallocate(wfnp)
  deallocate(hwfn)
  deallocate(wfn0)

end subroutine tdse_cycle_u1prop
!################################################################################
