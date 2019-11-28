! This file is part of xtb4stda.
!
! Copyright (C) 2015-2019 Stefan Grimme
!
! xtb4stda is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb4stda is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb4stda.  If not, see <https://www.gnu.org/licenses/>.

      subroutine shifteps(uhf,nmo,occa,occb,epsa,epsb,
     .                    lshift,lshifta,lshiftb)
      implicit none
      integer nmo
      logical uhf
      real*8 lshift,lshifta,lshiftb
      real*8 occa(nmo)
      real*8 occb(nmo)
      real*8 epsa(nmo)
      real*8 epsb(nmo)

      integer i
      real*8 dum

      if(uhf)then
! alpha
       do i=1,nmo
          dum=0.0d0
          if(occa(i).lt.1.d-6) dum=lshift              ! virt
          if(abs(occa(i)+occb(i)-1.0d0).lt.1.d-6) dum=dum-lshifta !  singly occ.
          epsa(i) = epsa(i) + dum
       enddo
! beta
       do i=1,nmo
          dum=0.0d0
          if(occb(i).lt.1.d-6) dum=lshift              ! virt
          if(abs(occa(i)+occb(i)-1.0d0).lt.1.d-6) dum=dum+lshiftb !  singly occ.
          epsb(i) = epsb(i) + dum
       enddo

      else

       do i=1,nmo
          if(occa(i)+occb(i).lt.1.d-6) epsa(i)=epsa(i)+lshift
       enddo

      endif

      end
