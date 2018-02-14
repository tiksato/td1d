!######################################################################
subroutine x2e_read(iunit, wfn)

  use init_mod, only : guess_type
  use grid_mod, only : ngrid, x

  implicit none
  integer, intent(in) :: iunit
  complex(kind(0d0)), intent(out) :: wfn(0:ngrid, 0:ngrid)

  integer :: ngrid0, nval, nshift
  integer :: igrid, jgrid, igrid0, jgrid0
  real(kind(0d0)), allocatable :: xchoi(:)
  complex(kind(0d0)), allocatable :: wfnchoi(:,:)

  read(iunit, "(I25)") ngrid0
  allocate(xchoi(0:ngrid0))
  allocate(wfnchoi(0:ngrid0, 0:ngrid0))

  read(iunit, "( E25.15)") xchoi(0:ngrid0)
  read(iunit, "(2E25.15)") wfnchoi(0:ngrid0, 0:ngrid0)

  if (trim(guess_type) == 'READOPT') then
!nyi
     write(6, "('WARNING: x2e_readproj not yet implemented.')")
!nyi
     nshift = (ngrid - ngrid0) / 2
     do igrid0 = 0, ngrid0
        igrid = igrid0 + nshift
        do jgrid0 = 0, ngrid0
           jgrid = jgrid0 + nshift
           wfn(igrid, jgrid) = wfnchoi(igrid0, jgrid0)
        end do
     end do
  else
!old     wfn(0:ngrid, 0:ngrid) = wfnchoi(0:ngrid, 0:ngrid)
     nshift = (ngrid - ngrid0) / 2
     do igrid0 = 0, ngrid0
        igrid = igrid0 + nshift
        do jgrid0 = 0, ngrid0
           jgrid = jgrid0 + nshift
           wfn(igrid, jgrid) = wfnchoi(igrid0, jgrid0)
        end do
     end do
  end if

  deallocate(wfnchoi)
  deallocate(xchoi)

end subroutine x2e_read
!######################################################################
