!################################################################################
subroutine init_cycle3(wfn, imethod)

  use const_mod, only : one, two
  use thresh_mod, only : threne
  use init_mod, only : maxcyc
  use wfn_mod, only : energy

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  integer icycle
  real(kind(0d0)) :: ene_ci, ene_mo, test

  ene_ci = one
  ene_mo = two
  ! alternate mo-ci optimization
  do icycle = 1, maxcyc

     call init_cycle1(wfn, imethod, .true., .false., .true.)
     ene_ci = energy
     call init_cycle1(wfn, imethod, .false., .true., .false.)
     ene_mo = energy

     test = ene_mo - ene_ci
     if(abs(test) < threne) exit

  end do

  energy = ene_mo

end subroutine init_cycle3
!################################################################################
