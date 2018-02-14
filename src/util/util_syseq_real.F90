!################################################################################
subroutine util_syseq_real(n, a, b)

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  real(kind(0d0)), intent(inout) :: b(1:n)

  integer :: info
  integer, allocatable :: ipiv(:)

  allocate(ipiv(1:n))

  call dgesv(n, 1, a, n, ipiv, b, n, info)

  if (info /= 0) then
     write(6, "('util_syseq_real: info = ', i20)") info
     stop 'error in util_syseq_real.'
  end if
  
  deallocate(ipiv)

end subroutine util_syseq_real
!################################################################################
