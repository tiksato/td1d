!################################################################################
subroutine util_index_double(dsc, n, vec, index)

  use const_mod, only : one

  implicit none
  logical, intent(in) :: dsc
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: vec(1:n)
  integer, intent(out) :: index(1:n)

  integer :: i
  real(kind(0d0)) :: tmp
  real(kind(0d0)), allocatable :: vec2(:)

!  n = size(vec)
!  n = ubound(vec)
!debug
!write(6,"('util_index_double: n = ', i10)") n
!stop
!debug
  allocate(vec2(1:n))
  vec2(1:n) = vec(1:n)

  if (dsc) then
     tmp = minval(vec2) - one
     do i = 1, n
        index(i) = maxloc(vec2, dim = 1)
        vec2(index(i)) = tmp
     end do
  else
     tmp = maxval(vec2) + one
     do i = 1, n
        index(i) = minloc(vec2, dim = 1)
        vec2(index(i)) = tmp
     end do
  end if

  deallocate(vec2)

end subroutine util_index_double
!################################################################################
