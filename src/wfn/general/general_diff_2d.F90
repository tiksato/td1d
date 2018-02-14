!######################################################################
subroutine general_diff_2d(wfn, dwfn, ohalf, coeff_in, coeff_l, coeff_r, ng0, ng1)

!     0,1,..,j,..,N-1,N
!   0 -----------------
!   1 -----------------
! ... -------c---------
! ng0 -------c---------
! ... -------c---------
!   i ----rrr*rrr------
! ... -------c---------
! ng1 -------c---------
! ... -------c---------
! N-1 -----------------
!   N -----------------

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: dwfn(0:ngrid, 0:ngrid)
  integer, intent(in) :: ohalf
  real(kind(0d0)), intent(in) :: coeff_in(-ohalf:ohalf)
  real(kind(0d0)), intent(in) :: coeff_l(0:(2*ohalf+1), 0:ohalf)
  real(kind(0d0)), intent(in) :: coeff_r(0:(2*ohalf+1), 0:ohalf)
  integer, intent(in) :: ng0, ng1

  complex(kind(0d0)) :: tmp
  integer :: k, kgrid, igrid, jgrid, icoeff, kmax

  kmax = 2 * ohalf + 1

! second derivative w.r.t. igrid

  ! i: inside
  do igrid = max(ng0, ohalf), min(ng1, ngrid - ohalf)
     do jgrid = 0, igrid
        tmp = wfn(jgrid, igrid) * coeff_in(0)
        do k = 1, ohalf
           tmp = tmp + wfn(jgrid, igrid - k) * coeff_in(-k) &
                   & + wfn(jgrid, igrid + k) * coeff_in(+k)
        end do
        dwfn(jgrid, igrid) = dwfn(jgrid, igrid) + tmp
     end do
  end do
!
!no boundary  ! i: left
!no boundary  do igrid = ng0, min(ng1, ohalf - 1)
!no boundary     icoeff = igrid
!no boundary     do jgrid = 0, igrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = k
!no boundary!2012.10.06           tmp = tmp + wfn(kgrid, jgrid) * coeff_l(k, icoeff)
!no boundary           tmp = tmp + wfn(jgrid, kgrid) * coeff_l(k, icoeff)
!no boundary        end do
!no boundary        dwfn(jgrid, igrid) = dwfn(jgrid, igrid) + tmp
!no boundary     end do
!no boundary  end do
!no boundary!
!no boundary  ! i: right
!no boundary  do igrid = max(ng0, ngrid - ohalf + 1), ng1
!no boundary     icoeff = ngrid - igrid     
!no boundary     do jgrid = 0, igrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = ngrid - kmax + k
!no boundary!2012.10.06           tmp = tmp + wfn(kgrid, jgrid) * coeff_r(k, icoeff)
!no boundary           tmp = tmp + wfn(jgrid, kgrid) * coeff_r(k, icoeff)
!no boundary        end do
!no boundary        dwfn(jgrid, igrid) = dwfn(jgrid, igrid) + tmp
!no boundary     end do
!no boundary  end do
!
! second derivative w.r.t. jgrid
!
  do igrid = ng0, ng1

     ! j: inside
     do jgrid = ohalf, min(igrid, ngrid - ohalf)
        tmp = wfn(jgrid, igrid) * coeff_in(0)
        do k = 1, ohalf
           tmp = tmp + wfn(jgrid - k, igrid) * coeff_in(-k) &
                   & + wfn(jgrid + k, igrid) * coeff_in(+k)
        end do
        dwfn(jgrid, igrid) = dwfn(jgrid, igrid) + tmp
     end do

!no boundary     ! j: left
!no boundary     do jgrid = 0, min(igrid, ohalf - 1)
!no boundary        icoeff = jgrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = k
!no boundary           tmp = tmp + wfn(kgrid, igrid) * coeff_l(k, icoeff)
!no boundary        end do
!no boundary        dwfn(jgrid, igrid) = dwfn(jgrid, igrid) + tmp
!no boundary     end do
!no boundary
!no boundary     ! j: right
!no boundary     do jgrid = ngrid - ohalf + 1, igrid
!no boundary        icoeff = ngrid - jgrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = ngrid - kmax + k
!no boundary           tmp = tmp + wfn(kgrid, igrid) * coeff_r(k, icoeff)
!no boundary        end do
!no boundary
!no boundary        dwfn(jgrid, igrid) = dwfn(jgrid, igrid) + tmp
!no boundary     end do

  end do

end subroutine general_diff_2d
!######################################################################
