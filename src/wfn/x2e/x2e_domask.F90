!######################################################################
subroutine x2e_domask(rmax, dt, wfn, q0, q1, q2)

  use omp_mod
  use grid_mod, only : ngrid, x, xmask, mask
  use const_mod, only : half

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: dt
  real(kind(0d0)), intent(inout) :: q0, q1, q2
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid
  real(kind(0d0)) :: absx, maski, maskij
  real(kind(0d0)) :: p0, p1, p2, ptot
  complex(kind(0d0)) :: dt2
  real(kind(0d0)), external :: x2e_norm

  dt2 = dt * half

  ptot = x2e_norm(rmax, wfn, p0, p1, p2)
  q0 = q0 + p0
  q1 = q1 + p1
  q2 = q2 + p2

!$omp parallel default(shared) private(igrid, jgrid, absx, maski, maskij)

  call omp_mod_thread(0, ngrid)

  do igrid = ng0, ng1
     absx = abs(x(igrid))
     if(absx > xmask) then

        maski = mask(igrid)

        do jgrid = 0, igrid - 1
           absx = abs(x(jgrid))
           if(absx > xmask) then
              maskij = maski * mask(jgrid)
           else
              maskij = maski
           end if
           wfn(jgrid, igrid) = wfn(jgrid, igrid) * maskij
           wfn(igrid, jgrid) = wfn(jgrid, igrid)
        end do

        maskij = maski * maski
        wfn(igrid, igrid) = wfn(igrid, igrid) * maskij

     else

        do jgrid = 0, igrid - 1
           absx = abs(x(jgrid))
           if(absx > xmask) then
              wfn(jgrid, igrid) = wfn(jgrid, igrid) * mask(jgrid)
              wfn(igrid, jgrid) = wfn(jgrid, igrid)
           end if
        end do

     end if
  end do

!$omp end parallel

  ptot = x2e_norm(rmax, wfn, p0, p1, p2)
  q0 = q0 - p0
  q1 = q1 - p1
  q2 = q2 - p2

end subroutine x2e_domask
!######################################################################
