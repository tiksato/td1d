!################################################################################
subroutine input_mol(ioin)

  use const_mod, only : two
  use root_mod, only : softnuc
  use mol_mod, only : maxatom, natom, ne, nae, nbe, mult, c, z, enen

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr, iatom, jatom
  namelist /mol/ natom, nae, nbe, mult, c, z, ne, enen

  natom = 1
  nae = 1
  nbe = 0
  mult = 0

  rewind(ioin)
  read(unit=ioin, nml=mol, iostat=ioerr)
!  if(ioerr /= 0) stop "error in namelist mol."

  if(natom <= 0 .or. natom > maxatom) stop "error in mol.natom"
  if(nae <= 0) stop "error in mol.nae"
  if(nbe < 0) stop "error in mol.nbe"

  ne(1) = nae
  ne(2) = nbe
  ne(3) = nae + nbe

  enen = 0.D+0
  do iatom = 1, natom
     do jatom = 1, iatom - 1
        enen = enen + z(iatom) * z(jatom) / abs(c(iatom) - c(jatom))
!        enen = enen + z(iatom) * z(jatom) / sqrt((c(iatom) - c(jatom))**two + softnuc)
     end do
  end do

  write(6, nml=mol)

end subroutine input_mol
!################################################################################
