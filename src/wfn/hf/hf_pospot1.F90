!######################################################################
integer function hf_pospt1()

  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin

  implicit none

  hf_pospt1 = (ngrid + 1) * nfun * nspin + 1

end function hf_pospt1
!######################################################################
