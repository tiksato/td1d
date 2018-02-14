!######################################################################
subroutine x2e_dyson(wfn, wfn2e, dwfn)

  use omp_mod
  use const_mod, only : one, czero
  use grid_mod, only : ngrid, dgrid
  use wfn_mod, only : nfun

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(in) :: wfn2e(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: dwfn(0:ngrid, 1:nfun)

  integer :: nbas, nbnf
  integer :: ifun, igrid
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), external :: util_zdotc

  nbas = ngrid + 1
  nbnf = nbas * nfun

!$omp parallel default(shared) private(igrid, ifun, tmp)

  call omp_mod_thread(0, ngrid)

  do ifun = 1, nfun
     do igrid = ng0, ng1
        tmp = util_zdotc(nbas, wfn(0,ifun), 1, wfn2e(0, igrid), 1)
        dwfn(igrid, ifun) = tmp * dgrid
     end do
  end do

!$omp end parallel
  
end subroutine x2e_dyson
!######################################################################
