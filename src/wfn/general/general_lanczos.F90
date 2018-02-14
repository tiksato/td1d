!######################################################################
subroutine general_lanczos(time, dt, wfn)

  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nfcore

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 1:nfun)
  !--------------------------------------------------------------------
  integer :: ifun

!STOP "general_lanczos still involves bug, NYI for nfun > 1."

  do ifun = nfcore+1, nfun
     call general_lanczos1(time, dt, wfn(0,ifun))
    !call general_lanczos2(time, dt, wfn(0,ifun))
  end do

end subroutine general_lanczos
!######################################################################
subroutine general_lanczos1(time, dt, wfn)

  use grid_mod, only : ngrid

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid)
  !--------------------------------------------------------------------
  integer,parameter :: krylov_nmic = 10
  integer,parameter :: krylov_maxcyc = 18
  real(kind(0d0)),parameter :: krylov_thresh = 1D-15
  real(kind(0d0)) :: norm1
  real(kind(0d0)) :: alph(1:krylov_maxcyc)
  real(kind(0d0)) :: beta(1:krylov_maxcyc)
  real(kind(0d0)), allocatable :: hmat(:,:)
  real(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)) :: dt1, timeh, val, tmp
  complex(kind(0d0)), allocatable :: func(:)
  complex(kind(0d0)), allocatable :: cvec(:,:)
  integer :: imic, igrid, nkry, isub, jsub, ii, nbas
  logical, parameter :: debug = .false.
  external :: general_hprod1
  complex(kind(0d0)), external :: general_ovlp1

  nbas = ngrid+1
  dt1 = dt/krylov_nmic
  allocate(cvec(0:ngrid, 1:krylov_maxcyc))

  do imic = 1, krylov_nmic
     ! generate Krylov subspace
     timeh = time + dt1*(imic-0.5d0)
     cvec(0:ngrid,1) = wfn(0:ngrid)
     !BUG call util_krylov(timeh, dt1, nbas, 1, krylov_maxcyc, krylov_thresh, &
     !BUG                  nkry, norm1, alph, beta, cvec, general_hprod1, general_ovlp1)
     call util_krylov2(timeh, dt1, nbas, krylov_maxcyc, krylov_thresh, &
                      nkry, norm1, alph, beta, cvec, general_hprod1, general_ovlp1)

     ! subspace exponentialization
     if (nkry > 0) then
        allocate(hmat(1:nkry, 1:nkry))
        allocate(uvec(1:nkry, 1:nkry))
        allocate(func(1:nkry))
        hmat(1:nkry, 1:nkry) = 0d0
        uvec(1:nkry, 1:nkry) = 0d0
        do isub = 1, nkry
           hmat(isub, isub) = alph(isub)
           if (isub < nkry) hmat(isub, isub + 1) = beta(isub)
           if (isub >    1) hmat(isub, isub - 1) = beta(isub-1)
        end do
        if (debug) then
           !write(6, "('Lanczos vectors')")
           !do isub = 1, nkry
           !   write(6, "(i5, 2f15.8)") isub, alph(isub), beta(isub)
           !end do
           write(6, "('CI Hamiltonian in ', i5, ' dimensional Krylov space:')") nkry
           do isub = 1, nkry
              do jsub = 1, nkry
                 write(6, "(f15.8)", advance = 'no') hmat(jsub, isub)
              end do
              write(6, *)
           end do
        end if

        call util_diag_real(.false., nkry, hmat, uvec)
        if (debug) write(6, "('Eigenvalues and operator diagonals in the Krylov space:')")

        func(1:nkry) = 0d0
        do isub = 1, nkry
          !val = exp(-(0d0,1d0)*hmat(isub,isub)*dt)
           val = exp(-(0d0,1d0)*hmat(isub,isub)*dt1)
           if (debug) write(6, "(i5, 3f15.8)") isub, hmat(isub, isub), val
           do jsub = 1, nkry
              func(jsub) = func(jsub) + uvec(jsub, isub) * val * uvec(1, isub) * norm1
           end do
        end do

        if (debug) then
           write(6, "('Weights of Krylov vectors in the solution:')")
           do isub = 1, nkry
              write(6, "(i5, 2f15.8)") isub, func(isub)
           end do
        end if

        wfn(:) = 0d0
        do isub = 1, nkry
           !$omp parallel default(shared)
           !$omp do
           do igrid = 0, ngrid
              wfn(igrid) = wfn(igrid) + cvec(igrid, isub) * func(isub)
           end do
           !$omp end do
           !$omp end parallel
        end do

        deallocate(func)
        deallocate(uvec)
        deallocate(hmat)
     end if
  end do

  deallocate(cvec)

end subroutine general_lanczos1
!######################################################################
subroutine general_lanczos2(time, dt, wfn)

  use grid_mod, only : ngrid

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn(1:(ngrid+1))
  external :: general_hprod1
  complex(kind(0d0)), external :: general_ovlp1
  !--------------------------------------------------------------------
  integer,parameter :: krylov_nmic = 10
  integer,parameter :: krylov_maxcyc = 18
  real(kind(0d0)),parameter :: krylov_thresh = 1D-10
  real(kind(0d0)) :: norm1
  real(kind(0d0)) :: alph(1:krylov_maxcyc)
  real(kind(0d0)) :: beta(1:krylov_maxcyc)
  real(kind(0d0)), allocatable :: hmat(:,:)
  real(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)) :: dt1, timeh, val, tmp
  complex(kind(0d0)), allocatable :: func(:)
  complex(kind(0d0)), allocatable :: cvec(:,:)
  integer :: imic, idim, nkry, isub, jsub, ii
  integer :: lwfn
  logical, parameter :: debug = .false.

  lwfn = ngrid + 1
  dt1 = dt/krylov_nmic
  allocate(cvec(1:lwfn, 1:krylov_maxcyc))

  do imic = 1, krylov_nmic
     ! generate Krylov subspace
     timeh = time + dt1*(imic-0.5)
     cvec(1:lwfn,1) = wfn(1:lwfn)
     ! call krylov(timeh, dt1, lwfn, krylov_maxcyc, krylov_thresh, &
     !             nkry, norm1, alph, beta, cvec, hprod)
     call util_krylov(timeh, dt1, lwfn, 1, krylov_maxcyc, krylov_thresh, &
                      nkry, norm1, alph, beta, cvec, general_hprod1, general_ovlp1)

     ! subspace exponentialization
     if (nkry > 0) then
        allocate(hmat(1:nkry, 1:nkry))
        allocate(uvec(1:nkry, 1:nkry))
        allocate(func(1:nkry))
        hmat(1:nkry, 1:nkry) = 0d0
        uvec(1:nkry, 1:nkry) = 0d0
        do isub = 1, nkry
           hmat(isub, isub) = alph(isub)
           if (isub < nkry) hmat(isub, isub + 1) = beta(isub)
           if (isub >    1) hmat(isub, isub - 1) = beta(isub-1)
        end do
        if (debug) then
           write(6, "('Lanczos vectors')")
           do isub = 1, nkry
              write(6, "(i5, 2f15.8)") isub, alph(isub), beta(isub)
           end do
           write(6, "('CI Hamiltonian in ', i5, ' dimensional Krylov space:')") nkry
           do isub = 1, nkry
              do jsub = 1, nkry
                 write(6, "(f15.8)", advance = 'no') hmat(jsub, isub)
              end do
              write(6, *)
           end do
        end if

!        call util_dsyev(nkry, hmat, uvec)
        call util_diag_real(.false., nkry, hmat, uvec)
        if (debug) write(6, "('Eigenvalues and operator diagonals in the Krylov space:')")

        func(1:nkry) = 0d0
        do isub = 1, nkry
           val = exp(-(0d0,1d0)*hmat(isub,isub)*dt1)
           if (debug) write(6, "(i5, 3f15.8)") isub, hmat(isub, isub), val
           do jsub = 1, nkry
              func(jsub) = func(jsub) + uvec(jsub, isub) * val * uvec(1, isub) * norm1
           end do
        end do

        wfn = 0d0
        do isub = 1, nkry
           !$omp parallel default(shared)
           !$omp do
           do ii = 1, lwfn
              wfn(ii) = wfn(ii) + cvec(ii, isub) * func(isub)
           end do
           !$omp end do
           !$omp end parallel
        end do

        deallocate(func)
        deallocate(uvec)
        deallocate(hmat)
     end if
  end do

  deallocate(cvec)

end subroutine general_lanczos2
!######################################################################
