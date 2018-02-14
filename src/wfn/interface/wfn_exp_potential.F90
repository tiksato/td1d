!######################################################################
subroutine wfn_exp_potential(dstep, lfield, wfn, expwfn, imethod)

  use fft_mod, only : fft_wfnr, fft_wfnk

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: dstep
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)  :: wfn(*)
  complex(kind(0d0)), intent(out) :: expwfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_exp_potential(dstep, lfield, wfn, expwfn, fft_wfnr, fft_wfnk)
  else if (imethod == 0) then
     call hf_exp_potential(dstep, lfield, wfn, expwfn)
  end if

end subroutine wfn_exp_potential
!######################################################################
