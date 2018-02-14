!######################################################################
integer function hf_size()

  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin

  implicit none

  integer, external :: mp_size

  hf_size = (ngrid + 1) * nfun * nspin

end function hf_size
!######################################################################
