!######################################################################
subroutine general_int2e(nfun1, nfun2, nfun3, nfun4, wfn1, wfn2, veff, &
     & int2e, ng0, ng1)

  use const_mod, only : czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid

  implicit none
  integer, intent(in) :: nfun1, nfun2, nfun3, nfun4
  complex(kind(0d0)), intent(in) :: wfn1(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: wfn2(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: veff(0:ngrid, 1:nfun3, 1:nfun4)
  complex(kind(0d0)), intent(inout) :: int2e(1:nfun1, 1:nfun2, 1:nfun3, 1:nfun4)
  integer, intent(in) :: ng0, ng1

  integer :: ifun, jfun, kfun, lfun, igrid
  complex(kind(0d0)) :: rho
  complex(kind(0d0)), allocatable :: tmp(:,:)

  allocate(tmp(1:nfun3, 1:nfun4))

  do jfun = 1, nfun2
     do ifun = 1, nfun1
        call util_zcopy(nfun3*nfun4, czero, 0, tmp, 1)
        do igrid = ng0, ng1
           rho = conjg(wfn1(igrid, ifun)) &
                   & * wfn2(igrid, jfun) * dgrid
           if (abs(rho) > thrwfn) then
              do lfun = 1, nfun4
                 do kfun = 1, nfun3
                    tmp(kfun, lfun) = tmp(kfun, lfun) + rho * veff(igrid, kfun, lfun)
                 end do
              end do
           end if
        end do
        do lfun = 1, nfun4
           do kfun = 1, nfun3
              int2e(ifun, jfun, kfun, lfun) = &
            & int2e(ifun, jfun, kfun, lfun) + tmp(kfun, lfun)
           end do
        end do
     end do
  end do

  deallocate(tmp)

end subroutine general_int2e
!######################################################################
