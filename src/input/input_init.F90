!################################################################################
subroutine input_init(ioin)

  ! input for initial wavefunction calculation

  use const_mod, only : zero, one
  use init_mod, only : maxcyc, maxdav, distep, distepci, distepmo, distep1, distep2, &
       & distep3, distepx, init_type, diag_type, guess_type, prop_type_init, cic_11, &
       & broken_symmetry, guess_interpolate, guess_alter, init_ah_trust, init_ah_alpha, &
       & init_ah_xnew, init_ah_maxvec, init_ah_maxvec2, init_ah_couple, init_ah_norm

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr
  namelist /init/ maxcyc, maxdav, distep, distepci, distepmo, distep1, distep2, distep3, &
       & distepx, init_type, diag_type, guess_type, prop_type_init, cic_11, broken_symmetry, &
       & guess_interpolate, guess_alter, init_ah_trust, init_ah_alpha, init_ah_xnew, init_ah_maxvec, &
       & init_ah_maxvec2, init_ah_couple, init_ah_norm

  maxcyc = 100000           ! Maximum number of initial optimization cycles
  maxdav = 100              ! Maximum number of Davidson micro cycles
  distep = 0.1d+0           ! Imaginary propagation step size
  distepci = -one           ! Imaginary propagation step size for ci
  distepmo = -one           ! Imaginary propagation step size for mo
  distep1 = -one            ! Imaginary propagation step size for mop1
  distep2 = -one            ! Imaginary propagation step size for mop2
  distep3 = -one            ! Imaginary propagation step size for moq
  distepx = -one            ! Imaginary propagation step size for unitary occ-occ rotation
  init_type = 'PROP1'       ! DIAG_DIAG, DIAG_PROP, DIAG_PROP_GBT, PROP1, PROP2, PROP3
  diag_type = 'LANCZOS'        ! CI diagonalization method (FULL, DAVIDSON, OLSEN)
  guess_type = 'LCAO'       ! LCAO, READ, HF, CAS, APSG, X2E
  prop_type_init = 'RK4'    ! RK4, SPLIT_RK4, SPLIT_EXP, DIAG
  cic_11(1:2) = 1           ! CIC(i,j) is set 1. i=cic_11(1), j=cic_11(2)
  broken_symmetry = .false. ! Force broken symmetry UHF
  guess_interpolate = .false. ! Grid interpolation
  guess_alter = .false.     ! To alternate order of guess orbitals
  init_ah_xnew = 0          ! algorighm for generating new trial vectors
  init_ah_maxvec = 10       ! algorighm for generating new trial vectors
  init_ah_maxvec2(1:3) = 0  ! algorighm for generating new trial vectors
  init_ah_trust = 0.6d+0    ! Initial trust radius augmented Hessian method
  init_ah_alpha = 1.0d+0    ! Initial step size of augmented Hessian method
  init_ah_couple(1:3) = 1
  init_ah_norm = .true.

  rewind(ioin)
  read(unit=ioin, nml=init, iostat=ioerr)
  if(ioerr /= 0) stop "error in namelist init."

  call util_transchar(guess_type)
  call util_transchar(init_type)
  call util_transchar(diag_type)
  call util_transchar(prop_type_init)

  distep1 = max(distep1, distep)
  distep2 = max(distep2, distep)
  distep3 = max(distep3, distep)
  distepci = max(distepci, distep)
  distepmo = min(min(distep1, distep2), distep3)
  distep = min(distepci, distepmo)
  if (trim(init_type) /= 'AUGHESS' .and. distep < zero) stop 'error in init.distep.'
  if (trim(init_type) /= 'AUGHESS' .and. distepci < zero) stop 'error in init.distepci.'
  if (trim(init_type) /= 'AUGHESS' .and. distepmo < zero) stop 'error in init.distepmo.'
  if (trim(init_type) /= 'AUGHESS' .and. distep1 < zero) stop 'error in init.distep1.'
  if (trim(init_type) /= 'AUGHESS' .and. distep2 < zero) stop 'error in init.distep2.'
  if (trim(init_type) /= 'AUGHESS' .and. distep3 < zero) stop 'error in init.distep3.'

  write(6, nml=init)
!debug
!  stop 'for debug in input_init.'
!debug

end subroutine input_init
!################################################################################
