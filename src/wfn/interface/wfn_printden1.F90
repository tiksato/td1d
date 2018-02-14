!######################################################################
subroutine wfn_printden1(iunit, wfn, imethod)

  use const_mod, only : zero, two
  use grid_mod, only : ngrid, dgrid, x

  implicit none
  integer, intent(in) :: iunit, imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: igrid, jgrid
  real(kind(0d0)) :: denx
  real(kind(0d0)), parameter :: tiny = 1.E-20
  real(kind(0d0)), allocatable  :: den(:)
  complex(kind(0d0)), allocatable  :: wfn2e(:,:)
  integer :: poscic
  integer, external :: wfn_poscic

end subroutine wfn_printden1
!######################################################################
