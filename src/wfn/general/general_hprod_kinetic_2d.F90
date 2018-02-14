!######################################################################
subroutine general_hprod_kinetic_2d(lfield, wfn, hwfn, ng0, ng1)

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

  use root_mod, only : icomp
  use const_mod, only : czero, iunit
  use field_mod, only : gauge
  use grid_mod, only : ngrid, dgrid
  use grid_mod, only : fd_order, fd_ohalf, fd_coeff_in, fd_coeff_l, fd_coeff_r, &
                     & fd_coeff1_in, fd_coeff1_l, fd_coeff1_r

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 0:ngrid)
  integer, intent(in) :: ng0, ng1
  complex(kind(0d0)) :: tmp
  integer :: k, kgrid, igrid, jgrid, icoeff, kmax
  complex(kind(0d0)) :: coeff_in(-fd_ohalf:fd_ohalf)
  complex(kind(0d0)) :: coeff_l(0:(fd_order + 1), 0:fd_ohalf)
  complex(kind(0d0)) :: coeff_r(0:(fd_order + 1), 0:fd_ohalf)

  kmax = 2 * fd_ohalf + 1
  if (icomp /= 1 .or. trim(gauge) == 'L') then
     ! kinetic energy
     coeff_in(-fd_ohalf:fd_ohalf) = fd_coeff_in(-fd_ohalf:fd_ohalf)
  else
     ! kinetic energy + velocity gauge field
     coeff_in(-fd_ohalf:fd_ohalf) = fd_coeff_in (-fd_ohalf:fd_ohalf) &
               & - iunit * lfield * fd_coeff1_in(-fd_ohalf:fd_ohalf)
  end if

! second derivative w.r.t. igrid

  ! i: inside
  do igrid = max(ng0, fd_ohalf), min(ng1, ngrid - fd_ohalf)
     do jgrid = 0, igrid
        tmp = wfn(jgrid, igrid) * coeff_in(0)
        do k = 1, fd_ohalf
           tmp = tmp + wfn(jgrid, igrid - k) * coeff_in(-k) &
                   & + wfn(jgrid, igrid + k) * coeff_in(+k)
        end do
        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
     end do
  end do
!
!no boundary  ! i: left
!no boundary  do igrid = ng0, min(ng1, fd_ohalf - 1)
!no boundary     icoeff = igrid
!no boundary     do jgrid = 0, igrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = k
!no boundary!2012.10.06           tmp = tmp + wfn(kgrid, jgrid) * fd_coeff_l(k, icoeff)
!no boundary           tmp = tmp + wfn(jgrid, kgrid) * fd_coeff_l(k, icoeff)
!no boundary        end do
!no boundary        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
!no boundary     end do
!no boundary  end do
!no boundary!
!no boundary  ! i: right
!no boundary  do igrid = max(ng0, ngrid - fd_ohalf + 1), ng1
!no boundary     icoeff = ngrid - igrid     
!no boundary     do jgrid = 0, igrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = ngrid - kmax + k
!no boundary!2012.10.06           tmp = tmp + wfn(kgrid, jgrid) * fd_coeff_r(k, icoeff)
!no boundary           tmp = tmp + wfn(jgrid, kgrid) * fd_coeff_r(k, icoeff)
!no boundary        end do
!no boundary        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
!no boundary     end do
!no boundary  end do
!
! second derivative w.r.t. jgrid
!
  do igrid = ng0, ng1

     ! j: inside
     do jgrid = fd_ohalf, min(igrid, ngrid - fd_ohalf)
        tmp = wfn(jgrid, igrid) * coeff_in(0)
        do k = 1, fd_ohalf
           tmp = tmp + wfn(jgrid - k, igrid) * coeff_in(-k) &
                   & + wfn(jgrid + k, igrid) * coeff_in(+k)
        end do
        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
     end do

!no boundary     ! j: left
!no boundary     do jgrid = 0, min(igrid, fd_ohalf - 1)
!no boundary        icoeff = jgrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = k
!no boundary           tmp = tmp + wfn(kgrid, igrid) * fd_coeff_l(k, icoeff)
!no boundary        end do
!no boundary        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
!no boundary     end do
!no boundary
!no boundary     ! j: right
!no boundary     do jgrid = ngrid - fd_ohalf + 1, igrid
!no boundary        icoeff = ngrid - jgrid
!no boundary        tmp = czero
!no boundary        do k = 0, kmax
!no boundary           kgrid = ngrid - kmax + k
!no boundary           tmp = tmp + wfn(kgrid, igrid) * fd_coeff_r(k, icoeff)
!no boundary        end do
!no boundary
!no boundary        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
!no boundary     end do

  end do

end subroutine general_hprod_kinetic_2d
!######################################################################
