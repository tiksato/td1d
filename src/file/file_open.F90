!################################################################################
subroutine file_open

  use root_mod, only : name
  use io_mod, only : io_rho1

  implicit none
  character(len = 256) :: fname

!disabled  fname = trim(name)//".rho1"  
!disabled  open(unit = io_rho1, file = trim(fname), status = 'new', form = 'formatted')
!disabled  write(io_rho1, "('##### istep, time, rho1_fc, rho1_dc, rho1_a, rho1 #####')") 
!disabled  close(unit = io_rho1)

!nyi  fname = trim(name)//".den1"  
!nyi  open(unit = io_den1, file = trim(fname), status = 'new', form = 'formatted')
!nyi  fname = trim(name)//".ovlp0"  
!nyi  open(unit = io_ovlp0, file = trim(fname), status = 'new', form = 'formatted')
!nyi  fname = trim(name)//".ovlp1"  
!nyi  open(unit = io_ovlp1, file = trim(fname), status = 'new', form = 'formatted')
!nyi  fname = trim(name)//".dip"  
!nyi  open(unit = io_dip, file = trim(fname), status = 'new', form = 'formatted')
!nyi  fname = trim(name)//".acc"  
!nyi  open(unit = io_acc, file = trim(fname), status = 'new', form = 'formatted')

end subroutine file_open
!################################################################################
