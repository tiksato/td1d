!######################################################################
real(kind(0d0)) function wfn_norm(rmax, wfn, p0, p1, p2, imethod)

  use const_mod, only : zero
  use root_mod, only : icomp
  use grid_mod, only : norm_type
  use wfn_mod, only : ormas

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in)  :: wfn(*)
  real(kind(0d0)), intent(out) :: p0, p1, p2
  complex(kind(0d0)), allocatable  :: wfn2e(:)

  integer :: poscic, pospt1, lenx
  integer, external :: wfn_poscic, wfn_pospt1, wfn_size
  real(kind(0d0)), external :: x2e_norm, hf_norm
  real(kind(0d0)), external :: util_norm2e

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (trim(norm_type) == 'DIRECT') then
     lenx = wfn_size(-1)
     allocate(wfn2e(lenx))

     call wfn_wfn2e(wfn, wfn2e, imethod)
     wfn_norm = util_norm2e(rmax, wfn2e, p0, p1, p2)

     deallocate(wfn2e)
  else
     if (imethod == -1) then
        wfn_norm = x2e_norm(rmax, wfn, p0, p1, p2)
     else if (imethod == 0) then
        wfn_norm = hf_norm(rmax, wfn, p0, p1, p2)
     end if
  end if

end function wfn_norm
!######################################################################
