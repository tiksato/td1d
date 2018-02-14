!################################################################################
subroutine grid_v1

  use const_mod, only : zero, one, two, three
  use root_mod, only : isoftnuc, softnuc
  use grid_mod, only : ngrid, x, v1, gv1, yukawa
  use mol_mod, only: natom, c, z

  implicit none
  integer :: iatom, igrid
  real(kind(0d0)) :: xn, denom

  v1(0:ngrid) = zero
  gv1(0:ngrid) = zero

  ! santra-2006
  if(isoftnuc == 1) then
    do iatom = 1, natom
       do igrid = 0, ngrid
          xn = x(igrid) - c(iatom)
          denom = softnuc + abs(xn)
          v1(igrid) = v1(igrid) - z(iatom) / denom
          gv1(igrid) = gv1(igrid) + z(iatom) / denom ** two
       end do
    end do

  ! yukawa
  else if (isoftnuc == -1) then
    do iatom = 1, natom
       do igrid = 0, ngrid
          xn = x(igrid) - c(iatom)
          denom = sqrt((softnuc + xn * xn))
          v1(igrid) = v1(igrid) - z(iatom) / denom * exp(-yukawa * abs(xn))
          gv1(igrid) = gv1(igrid) + z(iatom) / denom ** three * (xn + yukawa * denom ** two)
       end do
    end do

  ! default
  else 
    do iatom = 1, natom
       do igrid = 0, ngrid
          xn = x(igrid) - c(iatom)
          denom = sqrt((softnuc + xn * xn))
          v1(igrid) = v1(igrid) - z(iatom) / denom
          gv1(igrid) = gv1(igrid) + z(iatom) * xn / denom ** three
       end do
    end do
  end if

end subroutine grid_v1
!################################################################################
