!######################################################################
subroutine util_krylov(time, dtime, dim, nvec, maxcyc, thresh, ntot, norm1, alph, beta, vec, hprod1, ovlp1)

  implicit none
  complex(kind(0d0)), intent(in) :: time, dtime
  real(kind(0d0)), intent(in) :: thresh
  integer, intent(in) :: dim, nvec, maxcyc
  integer, intent(out) :: ntot
  real(kind(0d0)), intent(out) :: norm1(1:nvec)
  real(kind(0d0)), intent(out) :: alph(1:maxcyc)
  real(kind(0d0)), intent(out) :: beta(1:maxcyc)
  complex(kind(0d0)), intent(inout) :: vec(1:dim, 1:maxcyc)
  external hprod1
  complex(kind(0d0)), external :: ovlp1
 !complex(kind(0d0)), external :: util_zdotc

  integer :: ivec, icyc, jcyc, isub, nvecx
  real(kind(0d0)) :: oores  ! one over residue (1/res)
  real(kind(0d0)) :: error  ! error?
  real(kind(0d0)), parameter :: one = 1d0
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable :: hvec(:)

  ! ##### initialization #####
  alph(1:maxcyc) = 0d0
  beta(1:maxcyc) = 0d0

  ! ##### normalize first vectors #####
  nvecx = 0
  do ivec = 1, nvec
!     tmp = util_zdotc(dim, vec(1,ivec), 1, vec(1,ivec), 1)
!     tmp = util_zdotc(dim, vec(1,ivec), 1, vec(1,ivec), 1) * dgrid
     tmp = ovlp1(-one, vec(1,ivec), vec(1,ivec))
     norm1(ivec) = sqrt(abs(tmp))
     if (abs(norm1(ivec)) > thresh) then
        !write(6, "('krylov: norm1:', i5, e20.10)") ivec, norm1(ivec)
        nvecx = nvecx + 1
        tmp = 1d0/norm1(ivec)
        vec(1:dim,ivec) = vec(1:dim,ivec)*tmp
     end if
  end do
  if (nvecx == 0) return
  if (nvecx .ne. nvec) stop 'util_krylov: vanishing and non-vanishing first vectors...'

  ntot = nvecx
  allocate(hvec(1:dim))
  do icyc = 1, maxcyc

     ! new sigma vectors
     hvec = 0d0
     call hprod1(time, vec(1,icyc), hvec)

     ! new trial functions
!     tmp = util_zdotc(dim, vec(1,icyc), 1, hvec, 1)
!     tmp = util_zdotc(dim, vec(1,icyc), 1, hvec, 1) * dgrid
     tmp = ovlp1(-one, vec(1,icyc), hvec)
     alph(icyc) = dble(tmp)

     if (icyc == 1) then
        hvec(1:dim) = &
        hvec(1:dim) - vec(1:dim, icyc) * alph(icyc)
     else
        hvec(1:dim) = &
        hvec(1:dim) - vec(1:dim, icyc) * alph(icyc) &
                    - vec(1:dim, icyc - 1) * beta(icyc - 1)
        ! ##### EXACT orthogonalization ######################
        do jcyc = icyc - 2, 1, -1
!           tmp = util_zdotc(dim, vec(1,jcyc), 1, hvec, 1)
!           tmp = util_zdotc(dim, vec(1,jcyc), 1, hvec, 1) * dgrid
           tmp = ovlp1(-one, vec(1,jcyc), hvec)
           hvec(1:dim) = hvec(1:dim) - vec(1:dim, jcyc) * tmp
        end do
        ! ####################################################
     end if

!     tmp = util_zdotc(dim, hvec, 1, hvec, 1)
!     tmp = util_zdotc(dim, hvec, 1, hvec, 1) * dgrid
     tmp = ovlp1(-one, hvec, hvec)
     beta(icyc) = sqrt(dble(tmp))

     error = (icyc - 1) * log(abs(dtime))
     do isub = 1, icyc - 1
        error = error - log(dble(isub))
        error = error + log(beta(isub))
     end do
     error = error + log(beta(icyc))
     error = exp(error)
     error = error * error
     !DEBUG
     write(6, "('krylov: ', i10, 2f20.10, e20.10)") icyc, alph(icyc), beta(icyc), error
     !DEBUG

     if (error < thresh .and. icyc >= nvec) then
        !DEBUG
        write(6, "('krylov: converged', i10, e20.10)") icyc, error
        !DEBUG
        ntot = icyc
        exit
     else if (error < thresh) then
        !DEBUG
        write(6, "('krylov: converged', i10, e20.10)") icyc, error
        !DEBUG
     else if (icyc == maxcyc) then
        write(6, "('krylov: not converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else if (icyc > nvecx) then
        oores = 1d0 / beta(icyc)
        vec(1:dim,icyc+1) = hvec(1:dim) * oores
     end if
  end do

  deallocate(hvec)

end subroutine util_krylov
!######################################################################
subroutine util_krylov2(time, dtime, dim, maxcyc, thresh, ntot, norm1, alph, beta, vec, hprod1, ovlp1)
                        
  implicit none
  complex(kind(0d0)), intent(in) :: time, dtime
  real(kind(0d0)), intent(in) :: thresh
  integer(8), intent(in) :: dim, maxcyc
  integer(8), intent(out) :: ntot
  real(kind(0d0)), intent(out) :: norm1
  real(kind(0d0)), intent(out) :: alph(1:maxcyc)
  real(kind(0d0)), intent(out) :: beta(1:maxcyc)
  complex(kind(0d0)), intent(inout) :: vec(1:dim, 1:*)
  external hprod1
  complex(kind(0d0)), external :: ovlp1

  integer(8) :: icyc, jcyc, isub
  real(kind(0d0)) :: oores  ! one over residue (1/res)
  real(kind(0d0)) :: error  ! error?
  real(kind(0d0)), parameter :: one = 1d0
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable :: hvec(:)

  ! ##### initialization #####
  alph(1:maxcyc) = 0d0
  beta(1:maxcyc) = 0d0

  ! ##### normalize first vector #####
  !call zdotc_omp(dim, vec(1,1), vec(1,1), tmp)
  !tmp = util_zdot(dim, vec(1,1), vec(1,1))
  tmp = ovlp1(-one, vec(1,1), vec(1,1))
  norm1 = sqrt(abs(tmp))
  if (abs(norm1) < thresh) then
     !write(6, "('krylov: norm1 = ', e20.10)") norm1
     ntot = 0
     return
  end if
  tmp = 1d0 / norm1
  vec(1:dim, 1) = vec(1:dim, 1) * tmp

  allocate(hvec(1:dim))

  ntot = 1

  do icyc = 1, maxcyc

     ! new sigma vectors
     hvec(1:dim) = 0d0
     call hprod1(time, vec(1,icyc), hvec)

     ! new trial functions
     !call zdotc_omp(dim, vec(1, icyc), hvec, tmp)
     !tmp = util_zdot(dim, vec(1,icyc), hvec)
     tmp = ovlp1(-one, vec(1,icyc), hvec)
     alph(icyc) = dble(tmp)

     if (icyc == 1) then
        hvec(1:dim) = &
        hvec(1:dim) - vec(1:dim, icyc) * alph(icyc)
     else
        hvec(1:dim) = &
        hvec(1:dim) - vec(1:dim, icyc) * alph(icyc) &
                    - vec(1:dim, icyc - 1) * beta(icyc - 1)
        ! ##### EXACT orthogonalization ######################
        do jcyc = icyc - 2, 1, -1
           !call zdotc_omp(dim, vec(1, jcyc), hvec, tmp)
           !tmp = util_zdot(dim, vec(1,jcyc), hvec)
           tmp = ovlp1(-one, vec(1,jcyc), hvec)
           hvec(1:dim) = hvec(1:dim) - vec(1:dim, jcyc) * tmp
        end do
        ! ####################################################
     end if

     !call zdotc_omp(dim, hvec, hvec, tmp)
     !tmp = util_zdot(dim, hvec, hvec)
     tmp = ovlp1(-one, hvec, hvec)
     beta(icyc) = sqrt(dble(tmp))

     !error = krylov_error(icyc, alph, beta, dtime)
     !error = (icyc - 1) * log(dtime)
     error = (icyc - 1) * log(abs(dtime))
     do isub = 1, icyc - 1
        error = error - log(dble(isub))
        error = error + log(beta(isub))
     end do
     error = error + log(beta(icyc))
     error = exp(error)
     error = error * error
     !DEBUG
     !write(6, "('krylov: ', i10, 2f20.10, e20.10)") icyc, alph(icyc), beta(icyc), error
     !DEBUG

!    if (icyc == maxcyc .or. beta(icyc) < thresh) then
     if (error < thresh) then
        !DEBUG
        write(6, "('krylov: converged', i10, e20.10)") icyc, error
        !DEBUG
        ntot = icyc
        exit
     else if (icyc == maxcyc) then
        write(6, "('krylov: not converged', i10, e20.10)") icyc, error
        ntot = icyc
        exit
     else
        oores = 1d0 / beta(icyc)
        vec(1:dim, icyc + 1) = hvec(1:dim) * oores
!       vec(1:dim, icyc + 1) = hvec(1:dim)
     end if
  end do

  deallocate(hvec)

end subroutine util_krylov2
!######################################################################
