!######################################################################
subroutine wfn_printwfn2(iunit, wfn, imethod)

  implicit none
  integer, intent(in) :: iunit, imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: poscic, pospt1, lenx
  complex(kind(0d0)), allocatable  :: wfn2e(:)
  integer, external :: wfn_poscic, wfn_pospt1, wfn_size

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
     call x2e_print(iunit, wfn)
  else
     lenx = wfn_size(-1)
     allocate(wfn2e(lenx))

     if (imethod == 0) then
        call hf_wfn2e(wfn, wfn2e)
     end if

     call x2e_print(iunit, wfn2e)
     deallocate(wfn2e)
  end if

end subroutine wfn_printwfn2
!######################################################################
