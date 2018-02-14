!######################################################################
real(kind(0d0)) function x2e_n2err(wfn1, wfn2)

  use const_mod, only : one
  use grid_mod, only : ngrid

  implicit none
  complex(kind(0d0)), intent(in) :: wfn1(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(in) :: wfn2(0:ngrid, 0:ngrid)
  complex(kind(0d0)), external :: x2e_iprod

  x2e_n2err = abs(x2e_iprod(-one, wfn2, wfn2))
  
end function x2e_n2err
!######################################################################
