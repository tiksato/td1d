!################################################################################
subroutine prop_abm4(prop_type, time, dt, len, wfn0, wfn, hwfn, imethod)

  use const_mod, only : zero, one, two, three, six, half, iunit

  implicit none
  !--------------------------------------------------------------------
  character(len=*), intent(in) :: prop_type
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(in) :: wfn0(1:len)
  complex(kind(0d0)), intent(inout) :: wfn(1:len)
  complex(kind(0d0)), intent(inout) :: hwfn(1:len, 0:*)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: lfield, ene, error
  complex(kind(0d0)) :: dts
  real(kind(0d0)), external :: field
  real(kind(0d0)), external :: wfn_hprod
  real(kind(0d0)), external :: wfn_n2err
  complex(kind(0d0)), allocatable :: pwfn(:)
  complex(kind(0d0)), allocatable :: pder(:)
  complex(kind(0d0)), allocatable :: cwfn(:)

  dts = - iunit * dt / 24.d+0
  allocate(pwfn(1:len))
  allocate(pder(1:len))
  allocate(cwfn(1:len))

  ! predictor step: 4-th order adams-bashforth method
  pwfn(1:len) = -  9.d+0 * hwfn(1:len, 3) &
              & + 37.d+0 * hwfn(1:len, 2) &
              & - 59.d+0 * hwfn(1:len, 1) &
              & + 55.d+0 * hwfn(1:len, 0)
  pwfn(1:len) = pwfn(1:len) * dts
  pwfn(1:len) = pwfn(1:len) + wfn(1:len)

  ! estimated derivative
  lfield = field(time + dt)
  ene = wfn_hprod(.false., lfield, dt, wfn0, pwfn, pder, imethod)

  ! corrector step: 3-rd order adams-moulton method
  cwfn(1:len) =            hwfn(1:len, 2) &
              & -  5.d+0 * hwfn(1:len, 1) &
              & + 19.d+0 * hwfn(1:len, 0) &
              & +  9.d+0 * pder(1:len)
  cwfn(1:len) = cwfn(1:len) * dts
  cwfn(1:len) = cwfn(1:len) + wfn(1:len)

  ! error estimate
  if (trim(prop_type) == 'VABM') then
     pwfn(1:len) = cwfn(1:len) - pwfn(1:len)
     error = wfn_n2err(cwfn, pwfn, imethod)
     write(6, "('prop_abm4: error = ', f20.10)") error
  end if

  ! next value
  wfn(1:len) = cwfn(1:len)
  hwfn(1:len, 3) = hwfn(1:len, 2)
  hwfn(1:len, 2) = hwfn(1:len, 1)
  hwfn(1:len, 1) = hwfn(1:len, 0)

  deallocate(cwfn)
  deallocate(pder)
  deallocate(pwfn)

end subroutine prop_abm4
!################################################################################
