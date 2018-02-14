!################################################################################
module omp_mod

  implicit none

  integer :: nproc
  integer, parameter :: maxproc = 24
  integer, external :: omp_get_num_threads
  integer, external :: omp_get_thread_num
  integer :: ng0, ng1, ngpp, iproc
  integer :: llstr_omp, ulstr_omp, nstr_omp
!$omp threadprivate(ng0, ng1, ngpp, llstr_omp, ulstr_omp, nstr_omp, iproc)

  contains

!###########################################################
    subroutine omp_mod_thread(n0, n1)

      implicit none
      integer, intent(in) :: n0, n1

      integer :: ndata

      ndata = n1 - n0 + 1
      ngpp = ndata / nproc
      iproc = omp_get_thread_num()

      if (iproc == nproc - 1) then
         ng0 = ngpp * iproc + n0
         ng1 = n1
      else
         ng0 = ngpp * iproc + n0
         ng1 = ng0 + ngpp - 1
      end if

    end subroutine omp_mod_thread
!###########################################################
!###########################################################
    subroutine omp_mod_thread_str(llstr, ulstr)

      implicit none
      integer, intent(in) :: llstr, ulstr

      integer :: ndata

      ndata = ulstr - llstr + 1
      nstr_omp = ndata / nproc
      iproc = omp_get_thread_num()

      if (iproc == nproc - 1) then
         llstr_omp = nstr_omp * iproc + llstr
         ulstr_omp = ulstr
      else
         llstr_omp = nstr_omp * iproc + llstr
         ulstr_omp = llstr_omp + nstr_omp - 1
      end if

    end subroutine omp_mod_thread_str
!###########################################################

end module omp_mod
!################################################################################
