!######################################################################
subroutine general_diff(wfn, dwfn, ohalf, coeff_in, coeff_l, coeff_r, ng0, ng1)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid
  use wfn_mod, only : nfun, nfroz

  implicit none
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: dwfn(0:ngrid, 1:nfun)
  integer, intent(in) :: ohalf
  real(kind(0d0)), intent(in) :: coeff_in(-ohalf:ohalf)
  real(kind(0d0)), intent(in) :: coeff_l(0:(2*ohalf+1), 0:ohalf)
  real(kind(0d0)), intent(in) :: coeff_r(0:(2*ohalf+1), 0:ohalf)
  integer, intent(in) :: ng0, ng1

  complex(kind(0d0)) :: tmp
  integer :: k, kgrid, igrid, ifun, icoeff, kmax
!
  kmax = 2 * ohalf + 1
!
  do ifun = nfroz(1) + 1, nfun

!debug     do igrid = ng0, ng1
!debug        tmp = wfn(igrid, ifun) * coeff_in(0)
!debug        do k = 1, ohalf
!debug           if(igrid - k >= 0)     tmp = tmp + wfn(igrid - k, ifun) * coeff_in(-k)
!debug           if(igrid + k <= ngrid) tmp = tmp + wfn(igrid + k, ifun) * coeff_in(+k)
!debug        end do
!debug        dwfn(igrid, ifun) = dwfn(igrid, ifun) + tmp
!debug     end do

     ! inside ===============================================
     !         q <= igrid <= ngrid - q (ngrid - 2q + 1 terms)
     ! igrid - q <= kgrid <= igrid + q (        2q + 1 terms)
     ! ======================================================
     do igrid = max(ng0, ohalf), min(ng1, ngrid - ohalf)
        tmp = wfn(igrid, ifun) * coeff_in(0)
        do k = 1, ohalf
           tmp = tmp + wfn(igrid - k, ifun) * coeff_in(-k) &
                   & + wfn(igrid + k, ifun) * coeff_in(+k)
        end do
        dwfn(igrid, ifun) = dwfn(igrid, ifun) + tmp
     end do

!no boundary     ! left boundary =====================
!no boundary     ! 0 <= igrid <=  q - 1 (     q terms)
!no boundary     ! 0 <= kgrid <= 2q + 1 (2q + 2 terms)
!no boundary     ! ===================================
!no boundary     do igrid = ng0, min(ng1, ohalf - 1)
!no boundary        icoeff = igrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = k
!no boundary           tmp = tmp + wfn(kgrid, ifun) * coeff_l(k, icoeff)
!no boundary        end do
!no boundary        dwfn(igrid, ifun) = dwfn(igrid, ifun) + tmp
!no boundary     end do
!no boundary
!no boundary     ! right boundary ================================
!no boundary     ! ngrid -  q + 1 <= igrid <= ngrid (     q terms)
!no boundary     ! ngrid - 2q - 1 <= kgrid <= ngrid (2q + 2 terms)
!no boundary     ! ===============================================
!no boundary     do igrid = max(ng0, ngrid - ohalf + 1), ng1
!no boundary        icoeff = ngrid - igrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = ngrid - kmax + k
!no boundary           tmp = tmp + wfn(kgrid, ifun) * coeff_r(k, icoeff)
!no boundary        end do
!no boundary        dwfn(igrid, ifun) = dwfn(igrid, ifun) + tmp
!no boundary     end do

  end do

end subroutine general_diff
!######################################################################
