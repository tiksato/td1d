!######################################################################
subroutine wfn_exp_kinetic_fft(dstep, wfn, expwfn, imethod)

  use fft_mod, only : fft_wfnr, fft_wfnk
  use wfn_mod, only : nfun, nfcore, nspin

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: dstep
  complex(kind(0d0)), intent(in)  :: wfn(*)
  complex(kind(0d0)), intent(out) :: expwfn(*)

  integer :: poscic, lencic
  integer, external :: wfn_poscic, wfn_size

  if (imethod == -1) then
     call x2e_exp_kinetic_fft(dstep, wfn, expwfn, fft_wfnr, fft_wfnk)
  else if (imethod == 0) then
     call general_exp_kinetic_fft(dstep, nfun, nfcore, nspin, wfn, expwfn)
  end if

end subroutine wfn_exp_kinetic_fft
!######################################################################
