!################################################################################
subroutine input

  use io_mod, only : iostdi

  implicit none

!debug  call input_test
!debug  stop 'input_test.'

  call input_root(iostdi)
  call input_init(iostdi)
  call input_grid(iostdi)
  call input_mol(iostdi)
  call input_wfn(iostdi)
  call input_field(iostdi)
  call input_prop(iostdi)

!debug
!  stop 'for debug in input.'
!debug

end subroutine input
!################################################################################
!################################################################################
subroutine input_test


  complex(kind(0d0)) :: amat(3,3), uleft(3,3), uright(3,3)

  amat(1,1) = (1.d+0, 0.d+0); amat(1,2) = (4.d+0, 0.d+0); amat(1,3) = (7.d+0, 0.d+0)
  amat(2,1) = (2.d+0, 0.d+0); amat(2,2) = (5.d+0, 0.d+0); amat(2,3) = (8.d+0, 0.d+0)
  amat(3,1) = (3.d+0, 0.d+0); amat(3,2) = (6.d+0, 0.d+0); amat(3,3) = (9.d+0, 0.d+0)

!  amat(1,1) = (1.d+0, 0.d+0); amat(1,2) = (2.d+0, 0.d+0); amat(1,3) = (3.d+0, 0.d+0)
!  amat(2,1) = (2.d+0, 0.d+0); amat(2,2) = (4.d+0, 0.d+0); amat(2,3) = (5.d+0, 0.d+0)
!  amat(3,1) = (3.d+0, 0.d+0); amat(3,2) = (5.d+0, 0.d+0); amat(3,3) = (6.d+0, 0.d+0)

  write(6, "(3f20.10)") dble(amat(1,1)), dble(amat(1,2)), dble(amat(1,3))
  write(6, "(3f20.10)") dble(amat(2,1)), dble(amat(2,2)), dble(amat(2,3))
  write(6, "(3f20.10)") dble(amat(3,1)), dble(amat(3,2)), dble(amat(3,3))

  call util_gdiag_comp(.false., 3, amat, uleft, uright)
!  call util_diag_comp(.false., 3, amat, uleft)

  write(6, "(3f20.10)") dble(amat(1,1)), dble(amat(1,2)), dble(amat(1,3))
  write(6, "(3f20.10)") dble(amat(2,1)), dble(amat(2,2)), dble(amat(2,3))
  write(6, "(3f20.10)") dble(amat(3,1)), dble(amat(3,2)), dble(amat(3,3))

end subroutine input_test
!################################################################################
