!######################################################################
subroutine wfn_wfn2e(wfn, wfn2e, imethod)

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn(*)
  complex(kind(0d0)), intent(out) :: wfn2e(*)

  integer :: poscic, pospt1, lenx
  integer, external :: wfn_size, wfn_poscic, wfn_pospt1

  lenx = wfn_size(-1)
  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
     call util_zcopy(lenx, wfn, 1, wfn2e, 1)
  else if (imethod == 0) then
     call hf_wfn2e(wfn, wfn2e)
  end if

end subroutine wfn_wfn2e
!######################################################################
