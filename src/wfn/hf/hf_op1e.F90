!######################################################################
real(kind(0d0)) function hf_op1e(iop, rmax, wfn)

  use omp_mod
  use mol_mod, only : ne
  use const_mod, only : zero, one, two, czero, iunit
  use grid_mod, only : ngrid, dgrid, x, gv1, fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r
  use wfn_mod, only : nfun, nspin

  implicit none
  integer, intent(in) :: iop
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in)  :: wfn(0:ngrid, 1:nfun, 1:nspin)

  integer igrid, ifun, ispin, llgrid, ulgrid
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable  :: v(:)
  real(kind(0d0)) :: hffac
  complex(kind(0d0)), allocatable :: pwfn(:,:,:)

  if (nspin == 2 .or. ne(3) == 1) then
     hffac = one
  else
     hffac = two
  end if

  allocate(v(0:ngrid))
  allocate(pwfn(0:ngrid, 1:nfun, 1:nspin))

  hf_op1e = zero
  call get_irmax(rmax, llgrid, ulgrid, 0, ngrid)

  if (iop == 0) then
     ! dipole
     v(0:ngrid) = x(0:ngrid)
  else if (iop == 1) then
     ! acceleration
     v(0:ngrid) = gv1(0:ngrid)
  else if (iop == 2) then
     pwfn(0:ngrid, 1:nfun, 1:nspin) = czero
     !$omp parallel default(shared)
     call omp_mod_thread(llgrid, ulgrid)
     call general_diff(wfn(0,1,1), pwfn(0,1,1), fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r, ng0, ng1)
     if (nspin == 2) then
        call general_diff(wfn(0,1,2), pwfn(0,1,2), fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r, ng0, ng1)
     end if
     !$omp end parallel
  else
     stop 'hf_op1e: bad iop.'
  end if

!$omp parallel default(shared) private(tmp) reduction(+:hf_op1e)

  call omp_mod_thread(llgrid, ulgrid)

  if (iop /= 2) then
     tmp = czero
     do ispin = 1, nspin
        do ifun = 1, ne(ispin)
           do igrid = ng0, ng1
              tmp = tmp + conjg(wfn(igrid, ifun, ispin)) * v(igrid) &
                   &          * wfn(igrid, ifun, ispin)
           end do
        end do
     end do
     tmp = tmp * dgrid
     hf_op1e = hf_op1e + dble(tmp * hffac)
  else
     tmp = czero
     do ispin = 1, nspin
        do ifun = 1, ne(ispin)
           do igrid = ng0, ng1
              tmp = tmp + conjg(wfn(igrid, ifun, ispin)) &
                   &         * pwfn(igrid, ifun, ispin)
           end do
        end do
     end do
     tmp = - tmp * dgrid * iunit
     hf_op1e = hf_op1e + dble(tmp * hffac)
  end if
!$omp end parallel

  deallocate(pwfn)
  deallocate(v)

end function hf_op1e
!######################################################################
