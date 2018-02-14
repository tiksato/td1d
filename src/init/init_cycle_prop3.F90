!################################################################################
subroutine init_cycle_prop3(wfn, imethod)

  use const_mod, only : zero, one, two
  use thresh_mod, only : threne
  use init_mod, only : maxcyc
  use wfn_mod, only : energy

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
!nyi  integer icycle
!nyi  real(kind(0d0)) :: ene_ci, ene_mo, test
!nyi
!nyi  ene_ci = one
!nyi  ene_mo = two
!nyi  ! alternate mo-ci optimization
!nyi  do icycle = 1, maxcyc
!nyi
!nyi     call init_cycle_prop3_ci(wfn, imethod)
!nyi     ene_ci = energy
!nyi
!nyi     call init_cycle_prop3_mo(wfn, imethod)
!nyi     ene_mo = energy
!nyi
!nyi     test = ene_mo - ene_ci
!nyi     if(abs(test) < threne) exit
!nyi
!nyi  end do
!nyi
!nyi  energy = ene_mo
  energy = zero

end subroutine init_cycle_prop3
!################################################################################
