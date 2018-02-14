!################################################################################
subroutine grid_mkinner

  use mol_mod, only: natom, c
  use grid_mod, only : ngrid, x0, x, rmax_atom, grid_inner, llrmax, ulrmax, rmax, dgrid

  implicit none
  integer :: igrid, iatom
  real(kind(0d0)) :: dmin, dist

  llrmax = 0
  ulrmax = ngrid
  do igrid = 0, ngrid
     if (x(igrid) > -rmax-dgrid*0.1d0) then
        llrmax = igrid
        exit
     end if
  end do
  do igrid = ngrid, 0, -1
     if (x(igrid) < rmax+dgrid*0.1d0) then
        ulrmax = igrid
        exit
     end if
  end do

  ! as a real mask function
  grid_inner(0:ngrid) = 1
  do igrid = 0, ngrid
     dmin = x0
     do iatom = 1, natom
        dist = abs(x(igrid) - c(iatom))
        if (dist < dmin) dmin = dist
     end do
     if (dmin < rmax_atom) grid_inner(igrid) = 0
     !DEBUG
     !DEBUG write(6, "('grid_mkinner: ', i10, 2f20.5, i10)") igrid, x(igrid), dmin, grid_inner(igrid)
     !DEBUG
  end do

end subroutine grid_mkinner
!################################################################################
