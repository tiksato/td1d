!######################################################################
subroutine hf_ovlp(inner, rmax, wfnl, wfnr, slr)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid, gll, gul
  use wfn_mod, only : nfun, nspin

  implicit none
  logical, intent(in) :: inner
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: slr(1:nfun, 1:nfun, 1:nspin)

  complex(kind(0d0)) :: tmp
  integer :: igrid, ifun, jfun, ispin, llgrid, ulgrid

  slr(1:nfun, 1:nfun, 1:nspin) = czero
  call get_irmax(rmax, llgrid, ulgrid, gll, gul)

  if (inner) then
     do ispin = 1, nspin
        do ifun = 1, nfun
           do jfun = 1, nfun
              tmp = czero
              do igrid = llgrid, ulgrid
                 tmp = tmp + conjg(wfnl(igrid, ifun, ispin)) &
                         & *       wfnr(igrid, jfun, ispin)
              end do
              slr(ifun, jfun, ispin) = tmp * dgrid
           end do
        end do
     end do
  else
     do ispin = 1, nspin
        do ifun = 1, nfun
           do jfun = 1, nfun
              tmp = czero
              do igrid = 1, llgrid - 1
                 tmp = tmp + conjg(wfnl(igrid, ifun, ispin)) &
                               & * wfnr(igrid, jfun, ispin)
              end do
              do igrid = ulgrid + 1, gul
                 tmp = tmp + conjg(wfnl(igrid, ifun, ispin)) &
                               & * wfnr(igrid, jfun, ispin)
              end do
              slr(ifun, jfun, ispin) = tmp * dgrid
           end do
        end do
     end do
  end if

end subroutine hf_ovlp
!######################################################################
subroutine hf_ovlp_gbasis(wfnl, wfnr, ovlp)

  use const_mod, only : czero
  use grid_mod, only : ngbas
  use wfn_mod, only : nfun,nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfnl(1:ngbas,1:nfun,1:nspin)
  complex(kind(0d0)), intent(in) :: wfnr(1:ngbas,1:nfun,1:nspin)
  complex(kind(0d0)), intent(out) :: ovlp(1:nfun,1:nfun,1:nspin)

  complex(kind(0d0)) :: tmp
  integer :: mu,ifun,jfun,ispin

  !$omp parallel default(shared) private(tmp)
  !$omp do collapse(3)
  do ispin = 1, nspin
     do jfun = 1, nfun
        do ifun = 1, nfun
           tmp = czero
           do mu = 1, ngbas
              tmp = tmp + conjg(wfnl(mu,ifun,ispin)) &
                      & *       wfnr(mu,jfun,ispin)
           end do
           ovlp(ifun,jfun,ispin) = tmp
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine hf_ovlp_gbasis
!######################################################################
