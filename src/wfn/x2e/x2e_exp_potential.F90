!######################################################################
subroutine x2e_exp_potential(dstep, lfield, wfn, expwfn)

  use omp_mod
  use const_mod, only : iunit
  use root_mod, only : icomp, nocoulomb
  use grid_mod, only : ngrid, x, v1, v2, docap, mask

  implicit none
  complex(kind(0d0)), intent(in)  :: dstep
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)    :: wfn   (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: expwfn(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid
  complex(kind(0d0)) :: exp_veff
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(0:ngrid))

!$omp parallel default(shared) private(igrid, jgrid)

!old  call omp_mod_thread(1, ngrid - 1)
  call omp_mod_thread(0, ngrid)

  veff(0:ngrid) = v1(0:ngrid) - lfield * x(0:ngrid)
  if (docap .and. icomp == 1) then
     veff(0:ngrid) = veff(0:ngrid) - iunit * mask(0:ngrid)
  end if

  if (nocoulomb) then
     do igrid = ng0, ng1
        do jgrid = 1, igrid - 1
           exp_veff = veff(igrid) + veff(jgrid)
           exp_veff = exp(-iunit * dstep * exp_veff)
           expwfn(jgrid, igrid) = exp_veff * wfn(jgrid, igrid)
           expwfn(igrid, jgrid) = expwfn(jgrid, igrid)
        end do
        exp_veff = veff(igrid) + veff(igrid)
        exp_veff = exp(-iunit * dstep * exp_veff)
        expwfn(igrid, igrid) = exp_veff * wfn(igrid, igrid)
     end do
  else
     do igrid = ng0, ng1
        do jgrid = 1, igrid - 1
           exp_veff = veff(igrid) + veff(jgrid) + v2(jgrid, igrid)
           exp_veff = exp(-iunit * dstep * exp_veff)
           expwfn(jgrid, igrid) = exp_veff * wfn(jgrid, igrid)
           expwfn(igrid, jgrid) = expwfn(jgrid, igrid)
        end do
        exp_veff = veff(igrid) + veff(igrid) + v2(igrid, igrid)
        exp_veff = exp(-iunit * dstep * exp_veff)
        expwfn(igrid, igrid) = exp_veff * wfn(igrid, igrid)
     end do
  end if

!$omp end parallel

  deallocate(veff)

end subroutine x2e_exp_potential
!######################################################################
