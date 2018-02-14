!################################################################################
program td1d

  use root_mod, only : ical

  implicit none

!  call test
!  stop

  call input
  call grid_init
  call wfn_init
  call file_open

  call guess
  if (ical <= 0) call init
  if (ical >= 0) call tdse

  call file_close
  call wfn_final
  call grid_final

end program td1d
!################################################################################
!################################################################################
subroutine test

  integer :: i, j
  real(kind(0d0)) :: tmp(1:7, 1:7)
  real(kind(0d0)) :: exp_tmp(1:7, 1:7)

  tmp(1:7, 1:7) = 0.d+0
  exp_tmp(1:7, 1:7) = 0.d+0

  tmp(1,2) =   1.
  tmp(1,3) =  -1.
  tmp(1,4) =   2.
  tmp(1,5) =   3.
  tmp(1,6) =  -2.
  tmp(1,7) =  10.
  tmp(2,3) =   4.
  tmp(2,4) =  -4.
  tmp(2,5) =   5.
  tmp(2,6) =  -5.
  tmp(2,7) =  11.
  tmp(3,4) =   6.
  tmp(3,5) =  -6.
  tmp(3,6) =   7.
  tmp(3,7) =  12.
  tmp(4,5) =   8.
  tmp(4,6) =  -8.
  tmp(4,7) =  12.
  tmp(5,6) =   9.
  tmp(5,7) = -11.
  tmp(6,7) = -10.

  do i = 1, 7
     do j = 1, i -1
        tmp(i, j) = -tmp(j, i)
     end do
  end do

  call util_exp_asym(7, tmp, exp_tmp)

end subroutine test
!################################################################################
