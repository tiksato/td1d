!################################################################################
module const_mod

  implicit none

  ! numerical constants
  real(kind(0d0)), parameter :: zero  = 0.d+0
  real(kind(0d0)), parameter :: one   = 1.d+0
  real(kind(0d0)), parameter :: two   = 2.d+0
  real(kind(0d0)), parameter :: three = 3.d+0
  real(kind(0d0)), parameter :: four  = 4.d+0
  real(kind(0d0)), parameter :: five  = 5.d+0
  real(kind(0d0)), parameter :: six   = 6.d+0
  real(kind(0d0)), parameter :: eight = 8.d+0
  real(kind(0d0)), parameter :: half  = 5.d-1
  real(kind(0d0)), parameter :: quart = 2.5d-1
  real(kind(0d0)), parameter :: tenth = 1.d-1
  real(kind(0d0)), parameter :: pi = 3.141592653589793d+0
  real(kind(0d0)), parameter :: au2fs = 0.02418884254d+0
  complex(kind(0d0)), parameter :: czero = (zero, zero)
  complex(kind(0d0)), parameter :: runit = (one, zero)
  complex(kind(0d0)), parameter :: iunit = (zero, one)
  complex(kind(0d0)), parameter :: chalf = (half, zero)
  complex(kind(0d0)), parameter :: ctwo  = (two, zero)
  complex(kind(0d0)), parameter :: cfour = (four, zero)

  ! parameter
  integer, parameter :: maxproc = 256

end module const_mod
!################################################################################
