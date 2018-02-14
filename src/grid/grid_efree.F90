!#######################################################################
subroutine grid_efree

! p(ikpt)         = +dp * ikpt (0 <= ikpt <  ngrid/2)
! p(ngrid - ikpt) = -dp * ikpt (1 <= ikpt <= ngrid/2)
! dp = 2 * pi / lbox

  use const_mod, only : zero, two, half, half, pi
  use grid_mod, only : ngrid, x0, p, p2, efree, efree2

  implicit none
  integer :: kgrid, lgrid, ikpt
  real(kind(0d0)) :: lbox, fac, mom, ene

  lbox = two * x0
  fac = (two * pi) / lbox

  p(0:ngrid-1) = zero
  efree(0:ngrid-1) = zero
  efree2(0:ngrid-1, 0:ngrid-1) = zero

  if (mod(ngrid, 2) == 0) then
     do ikpt = 1, ngrid / 2 - 1
        mom = fac * dble(ikpt)
        ene = mom ** two * half

        p(ikpt) = mom
        p(ngrid - ikpt) = -mom

        efree(ikpt) = ene
        efree(ngrid - ikpt) = ene
     end do
     mom = fac * dble(ngrid / 2)
     ene = mom ** two * half

!old     p(ngrid / 2) = mom
     p(ngrid / 2) = -mom
     efree(ngrid / 2) = ene
  else
     do ikpt = 1, ngrid / 2
        mom = fac * dble(ikpt)
        ene = mom ** two * half

        p(ikpt) = mom
        p(ngrid - ikpt) = -mom

        efree(ikpt) = ene
        efree(ngrid - ikpt) = ene
     end do
  end if

  do kgrid = 0, ngrid - 1
     do lgrid = 0, kgrid - 1
        p2(lgrid, kgrid) = p(lgrid) + p(kgrid)
        p2(kgrid, lgrid) = p2(lgrid, kgrid)
        efree2(lgrid, kgrid) = efree(lgrid) + efree(kgrid)
        efree2(kgrid, lgrid) = efree2(lgrid, kgrid)
     end do
     p2(kgrid, kgrid) = p(kgrid) + p(kgrid)
     efree2(kgrid, kgrid) = efree(kgrid) + efree(kgrid)
  end do

end subroutine grid_efree
!#######################################################################
