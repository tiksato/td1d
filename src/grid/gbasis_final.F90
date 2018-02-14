!#########################################################
subroutine gbasis_final
  use grid_mod
  implicit none
  deallocate(govlp)
  deallocate(goinv)
  deallocate(gcore)
  deallocate(geris)
  deallocate(gmap2)
  deallocate(gmap4)
end subroutine gbasis_final
!#########################################################
