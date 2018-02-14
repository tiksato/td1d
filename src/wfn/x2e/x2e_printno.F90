!######################################################################
subroutine x2e_printno(iow, wfn)

  use const_mod, only : two, czero
  use grid_mod, only : ngrid, x, dgrid

  implicit none
  integer, intent(in) :: iow
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid, kgrid
  complex(kind(0d0)) :: fac
  complex(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)), allocatable :: den(:,:)

  allocate(uvec(0:ngrid, 0:ngrid))
  allocate(den(0:ngrid, 0:ngrid))

  call util_zcopy(ngrid*ngrid, czero, 0, uvec, 1)
  call util_zcopy(ngrid*ngrid, czero, 0, den, 1)

  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        do kgrid = 0, ngrid
           den(jgrid, igrid) = den(jgrid, igrid) + &
                & wfn(jgrid, kgrid) * conjg(wfn(igrid, kgrid))
        end do
     end do
  end do

  fac = two * dgrid * dgrid
  do igrid = 0, ngrid
     do jgrid = 0, ngrid
        den(jgrid, igrid) = den(jgrid, igrid) * fac
     end do
  end do

  call util_diag_comp(.true., ngrid+1, den, uvec)

  do igrid = 0, ngrid
     write(iow, "(i10, F20.10)") igrid+1, dble(den(igrid, igrid))
  end do

  deallocate(uvec)
  deallocate(den)

end subroutine x2e_printno
!######################################################################
