!######################################################################
real(kind(0d0)) function x2e_op1e(iop, rmax, wfn)

  use omp_mod
  use const_mod, only : two, zero, czero, iunit
  use grid_mod, only : ngrid, dgrid, x, gv1, fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r

  implicit none
  integer, intent(in) :: iop
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)

  integer igrid, jgrid, llgrid, ulgrid
  complex(kind(0d0)) :: val, tmp, v12
  complex(kind(0d0)), allocatable :: v(:)
  complex(kind(0d0)), allocatable :: pwfn(:,:)

  allocate(v(0:ngrid))
  if(iop == 0) then
     v(0:ngrid) = x(0:ngrid)
  else if (iop == 1) then
     v(0:ngrid) = gv1(0:ngrid)
  end if

  x2e_op1e = zero
  call get_irmax(rmax, llgrid, ulgrid, 0, ngrid)

  if (iop == 2) then
     allocate(pwfn(0:ngrid, 0:ngrid))
     pwfn(0:ngrid, 0:ngrid) = czero
     !$omp parallel default(shared)
     call omp_mod_thread(llgrid, ulgrid)
     call general_diff_2d(wfn, pwfn, fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r, ng0, ng1)
     !$omp end parallel
  end if

!$omp parallel default(shared) private(igrid, jgrid, val, tmp, v12) reduction(+:x2e_op1e)

  call omp_mod_thread(llgrid, ulgrid)

! don't work here, but why?      if(iop == 0) then
! don't work here, but why?         v(ng0:ng1) = x(ng0:ng1)
! don't work here, but why?      end if

  val = czero

  if (iop /= 2) then
     tmp = czero
     do igrid = ng0, ng1
        v12 = v(igrid) + v(igrid)
        tmp = tmp + conjg(wfn(igrid, igrid)) * v12 * wfn(igrid, igrid)
     end do
     val = val + tmp
     
     tmp = czero
     do igrid = ng0, ng1
        do jgrid = llgrid, igrid - 1
           v12 = v(igrid) + v(jgrid)
           tmp = tmp + conjg(wfn(jgrid, igrid)) * v12 * wfn(jgrid, igrid)
        end do
     end do
  else
     tmp = czero
     do igrid = ng0, ng1
        tmp = tmp + conjg(wfn(igrid, igrid)) * pwfn(igrid, igrid)
     end do
     val = val + tmp
     
     tmp = czero
     do igrid = ng0, ng1
        do jgrid = llgrid, igrid - 1
           tmp = tmp + conjg(wfn(jgrid, igrid)) * pwfn(jgrid, igrid)
        end do
     end do
  end if

  val = val + tmp * two
  val = val * dgrid**two
  if (iop == 2) val = - val * iunit

  x2e_op1e = x2e_op1e + dble(val)

!$omp end parallel

  if (iop == 2) then
     deallocate(pwfn)
  end if
      
  deallocate(v)

end function x2e_op1e
!######################################################################
