!################################################################################
subroutine input_root(ioin)

  use omp_mod
  use const_mod, only : zero, one
  use root_mod, only : name, ical, iprint, isoftnuc, isoftr12, softnuc, softr12, &
       & nocoulomb, nprint_ene, nprint_op1e, nprint_op1x, nprint_ionp, nprint_rho1, root_debug, use_blas
  use thresh_mod, only : threne, thrwfn, throcc, thrdet, threqn, thrgbt, thrgrd, thrdav

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr
  namelist /root/ name, ical, iprint, isoftnuc, isoftr12, softnuc, softr12, &
                & nocoulomb, threne, thrwfn, throcc, thrdet, threqn, thrgbt, thrgrd, thrdav, &
                & nproc, nprint_ene, nprint_op1e, nprint_op1x, nprint_ionp, nprint_rho1, root_debug, use_blas

  name = 'test'       ! basename of scratch files
  ical = 0            ! -1: init, 0: init+tdse, 1: tdse
  iprint = 0          ! amount of debug printing
  isoftnuc = 0        ! 0: 1 / sqrt(x**2 + d), 1: 1 / (|x| + d)
  isoftr12 = 0        ! 0: 1 / sqrt(x**2 + d), 1: 1 / (|x| + d)
  softnuc = one       ! parameter d of soft n-e potential
  softr12 = one       ! parameter d of soft e-e potential
  nocoulomb = .false. ! omit e-e (for debug)
  threne = 1.D-13     ! energy threshold
  thrwfn = 1.D-20     ! wavefunction threshold
  throcc = 1.D-5      ! occupation number denominator threshold (cas and apsg)
  thrdet = 1.D-12     ! orthonormality threshold used in det_detsx
  threqn = 1.D-15     ! iterative linear equation solver threshold
  thrgbt = 1.D-5      ! generalized brillouin condition threshold
  thrgrd = 1.D-5      ! energy gradient threshold
  thrdav = 1.D-6      ! norm of residual vector for davidson diagonalization
  nprint_ene = 1
  nprint_op1e = -1    ! printing op1e at every nprint_op1e step
  nprint_op1x = -1    ! printing domain-devided op1e at every nprint_op1x step
  nprint_ionp = -1    ! printing ionp at every nprint_ionp step
  nprint_rho1 = -1    ! printing rho1 at every nprint_rho1 step
  root_debug = .false.
  use_blas = .true.
  
  rewind(ioin)
  read(unit=ioin, nml=root, iostat=ioerr)
  if (ioerr /= 0) stop "error in namelist root."
!debug  if (ioerr /= 0) then
!debug     write(6, "('input_root: ioerr = ', i5)") ioerr
!debug  end if

!$omp parallel default(shared)
  if (iproc == 0) nproc = omp_get_num_threads()
!$omp end parallel

  if (isoftnuc /=0 .and. isoftnuc /= 1 .and. isoftnuc /= -1) stop "error in root.isoftnuc."
  if (isoftr12 /=0 .and. isoftr12 /= 1 .and. isoftr12 /= -1) stop "error in root.isoftr12."
  if (threne < zero) stop "error in thresh.threne."
  if (thrwfn < zero) stop "error in thresh.thrwfn."
  if (throcc < zero) stop "error in thresh.throcc."
  if (thrdav < zero) stop "error in thresh.thrdav."
  if (nproc > maxproc) stop "error in nproc = $OMP_NUM_THREADS."

  write(6, nml=root)
!debug
!  stop 'stop for debug in input_root.'
!debug

end subroutine input_root
!################################################################################
