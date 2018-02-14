!################################################################################
subroutine guess_core(wfn, imethod)

  use grid_mod, only : ngbas,govlp,gcore
  use wfn_mod, only : nfun,nspin

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(1:ngbas, 1:nfun, 1:nspin)
  !--------------------------------------------------------------------
  integer :: mu,ifun,ispin
  real(kind(0d0)) :: xn, x2, alpha, gauss
  real(kind(0d0)), allocatable :: ctmp(:,:), umat(:,:)

  allocate(ctmp(1:ngbas,1:ngbas))
  allocate(umat(1:ngbas,1:ngbas))

  ctmp = gcore
  !call util_diag_real2(.false.,ngbas,govlp,ctmp,umat)
  call util_diag_real(.false.,ngbas,ctmp,umat)

  do ispin = 1, nspin
     do ifun = 1, nfun
        wfn(1:ngbas,ifun,ispin)=umat(1:ngbas,ifun)
     end do
  end do

  !DEBUG
  !write(6,"('Eigenvalues of HCORE:')")
  !do mu = 1, ngbas
  !   write(6,"(i5,f20.10)") mu,ctmp(mu,mu)
  !end do
  !stop "for debug @ guess_core"
  !DEBUG

  call hf_ort_gbasis(wfn)

  deallocate(umat)
  deallocate(ctmp)

end subroutine guess_core
!######################################################################
