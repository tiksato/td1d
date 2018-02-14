!################################################################################
subroutine grid_v2
  use const_mod, only : zero, one, two
  use root_mod, only : isoftr12, softr12
  use grid_mod, only : ngrid, x, v2, yukawa

  implicit none
  integer :: igrid, jgrid

  ! santra-2006
  if(isoftr12 == 1) then
     do igrid = 0, ngrid
        do jgrid = 0, ngrid
           v2(jgrid, igrid) = one / (softr12 + abs(x(igrid) - x(jgrid)))
        end do
     end do

  ! yukawa
  else if (isoftr12 == -1) then
     do igrid = 0, ngrid
        do jgrid = 0, ngrid
           v2(jgrid, igrid) = one / sqrt((softr12 + (x(igrid) - x(jgrid))**two)) &
                & * exp(-yukawa * abs(x(igrid) - x(jgrid)))
        end do
     end do

  ! default
  else 
     do igrid = 0, ngrid
        do jgrid = 0, ngrid
           v2(jgrid, igrid) = one / sqrt((softr12 + (x(igrid) - x(jgrid))**two))
        end do
     end do
  end if

!  write(6, "('WARNING: v2 is set ZERO!')")
!  do igrid = 0, ngrid
!     do jgrid = 0, ngrid
!        v2(jgrid, igrid) = zero
!     end do
!  end do

!test
!  call grid_v2_diag
!  stop 'debug for v2_diag'
!test

end subroutine grid_v2
!################################################################################
!################################################################################
subroutine grid_v2_diag

  use const_mod, only : one, two
  use grid_mod, only : ngrid, x, v2

  implicit none
  integer :: igrid, ng0, ngng0
  real(kind(0d0)), allocatable :: v2e(:,:), v2u(:,:)

  allocate(v2e(0:ngrid,0:ngrid))
  allocate(v2u(0:ngrid,0:ngrid))

  ng0 = ngrid + 1
  ngng0 = ng0 * ng0

  call dcopy(ngng0, v2, 1, v2e, 1)
  call util_diag_real(.false., ng0, v2e, v2u)

  write(6,"(10E17.10)") (v2e(igrid, igrid), igrid = 0, ngrid)

  deallocate(v2u)
  deallocate(v2e)  

end subroutine grid_v2_diag
!################################################################################
