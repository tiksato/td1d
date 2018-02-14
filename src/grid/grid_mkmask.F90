!################################################################################
subroutine grid_mkmask

  use prop_mod, only : dstep
  use const_mod, only : zero, one, two, four, pi, half, eight
  use grid_mod, only : ngrid, x0, x, domask, docap, xmask, intmask, mask_type, mask, dgrid

  implicit none
  integer :: igrid
  real(kind(0d0)) :: x2, absx, pihalf, oo4, oo8

  pihalf = pi * half
  oo4 = one / four
  oo8 = one / eight

  ! as a real mask function
  if (domask) then
     mask(0:ngrid) = one
     do igrid = 0, ngrid
        absx = abs(x(igrid))
        if(absx > xmask+dgrid*0.1d0) then

           ! quadratic
           if (trim(mask_type) == 'QUADRATIC') then
              x2 = (absx - xmask)**two
!              mask(igrid) = exp(-intmask * x2)
              mask(igrid) = exp(-intmask * x2 * dstep)
            
           ! ln[cos**(1/4)]
           else if (trim(mask_type) == 'COS4') then
              if (igrid > 0 .and. igrid < ngrid) then
                 mask(igrid) = cos((absx - xmask) / (x0 - xmask) * pihalf) ** oo4
              else
                 mask(igrid) = zero
              end if

           ! ln[cos**(1/8)]
           else if (trim(mask_type) == 'COS8') then
              if (igrid > 0 .and. igrid < ngrid) then
                 mask(igrid) = cos((absx - xmask) / (x0 - xmask) * pihalf) ** oo8
              else
                 mask(igrid) = zero
              end if
           end if

        end if
     end do

  ! as a complex absorbing potential
  else if (docap) then
     mask(0:ngrid) = zero
     do igrid = 0, ngrid
        absx = abs(x(igrid))
        if(absx > xmask+dgrid*0.1d0) then

           ! quadratic
           if (trim(mask_type) == 'QUADRATIC') then
              x2 = (absx - xmask)**two
              mask(igrid) = intmask * x2

           ! cos**(1/4)
           else if (trim(mask_type) == 'COS4') then
              mask(igrid) = log(-cos((absx - xmask) / (x0 - xmask) * pihalf) ** oo4)

           ! cos**(1/8)
           else if (trim(mask_type) == 'COS8') then
              mask(igrid) = log(-cos((absx - xmask) / (x0 - xmask) * pihalf) ** oo8)
           end if

        end if
     end do
  end if

end subroutine grid_mkmask
!################################################################################
