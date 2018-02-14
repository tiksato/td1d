!################################################################################
subroutine input_prop(ioin)

! input propagation settings

  use const_mod, only : zero, one, two, pi, czero
  use prop_mod, only : prop_type, cyctot, nstep, totstep, dstep, initime, inistep, &
       & prop_tol, prop_safety, prop_ulscal, dstep_ll, dstep_ul
  use field_mod, only : period
  use fft_mod, only : dofft

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr
  namelist /prop/ prop_type, cyctot, nstep, dstep, dofft, totstep, initime, inistep, &
       & prop_tol, prop_safety, prop_ulscal, dstep_ll, dstep_ul

  prop_type = 'RK4'     ! ABM, RK4, VRK5, EULER, SPLIT_RK4, SPLIT_EXP, ARNOLDI
  cyctot = 0            ! total simulation time in unit of period
  nstep = -1            ! number of steps in a period
  totstep = -1          ! total number of steps (for debug)
  dofft = .false.       ! use FFT for kinetic energy
  prop_tol = 1.D-15     ! error tolerance for variable step-size propagator
  prop_safety = 0.50d+0 ! safety factor for vrk propagator
  prop_ulscal = 1.20d+0 ! maximum allowed enhancement factor in vrk
  initime = czero       ! initial time
  inistep = 0           ! initial step
  dstep_ll = -one
  dstep_ul = -one

  rewind(ioin)
  read(unit=ioin, nml=prop, iostat=ioerr)
  if(ioerr /= 0) stop "error in namelist prop."

  call util_transchar(prop_type)

  if (nstep < 0) stop "error in prop.nstep."
  if (totstep < 0) totstep = nstep * cyctot
  dstep = period / nstep;
  if (dstep_ul < zero) dstep_ul = dstep

  call input_prop_setrk

  write(6, nml=prop)

end subroutine input_prop
!################################################################################
!################################################################################
subroutine input_prop_setrk

  use prop_mod, only :rk5_a, rk5_b, rk5_crk5, rk5_crk4

  implicit none

  rk5_a(2) = 1.d0/5.d0
  rk5_a(3) = 3.d0/10.d0
  rk5_a(4) = 3.d0/5.d0
  rk5_a(5) = 1.d0
  rk5_a(6) = 7.d0/8.d0

  rk5_b(1:5, 2) = (/   1.d0/5.d0,       0.d0,          0.d0,              0.d0,             0.d0        /)
  rk5_b(1:5, 3) = (/   3.d0/40.d0,      9.d0/40.d0,    0.d0,              0.d0,             0.d0        /)
  rk5_b(1:5, 4) = (/   3.d0/10.d0,     -9.d0/10.d0,    6.d0/5.d0,         0.d0,             0.d0        /)
  rk5_b(1:5, 5) = (/ -11.d0/54.d0,      5.d0/2.d0,   -70.d0/27.d0,       35.d0/27,          0.d0        /)
  rk5_b(1:5, 6) = (/1631.d0/55296.d0, 175.d0/512.d0, 575.d0/13824.d0, 44275.d0/110592.d0, 253.d0/4096.d0/)

  rk5_crk5(1) =  37.d0/378.d0
  rk5_crk5(2) =   0.d0
  rk5_crk5(3) = 250.d0/621.d0
  rk5_crk5(4) = 125.d0/594.d0
  rk5_crk5(5) =   0.d0
  rk5_crk5(6) = 512.d0/1771.d0

  rk5_crk4(1) =  2825.d0/27648.d0
  rk5_crk4(2) =     0.d0
  rk5_crk4(3) = 18575.d0/48384.d0
  rk5_crk4(4) = 13525.d0/55296.d0
  rk5_crk4(5) =   277.d0/14336.d0
  rk5_crk4(6) =     1.d0/4.d0

end subroutine input_prop_setrk
!################################################################################
