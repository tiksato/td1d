!################################################################################
subroutine get_irmax(rmax, llgrid, ulgrid, ng0, ng1)

  use const_mod, only : zero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, x

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: rmax
  integer, intent(out) :: llgrid, ulgrid
  integer, intent(in) :: ng0, ng1
  !--------------------------------------------------------------------
  integer :: igrid
  real(kind(0d0)) :: test

  llgrid = ng0
  ulgrid = ng1

  if (rmax < zero) return

  do igrid = ng0, ng1
     test = abs(x(igrid)) - rmax
     if(test < zero .or. abs(test) < thrwfn) then
!     if(test < thrwfn .or. abs(test) < thrwfn) then
        llgrid = igrid
        exit
     end if
  end do
  do igrid = ng1, ng0, -1
     test = abs(x(igrid)) - rmax
     if(test < zero .or. abs(test) < thrwfn) then
!     if(test < thrwfn .or. abs(test) < thrwfn) then
        ulgrid = igrid
        exit
     end if
  end do

  llgrid = max(llgrid, ng0)
  ulgrid = min(ulgrid, ng1)

end subroutine get_irmax
!################################################################################
