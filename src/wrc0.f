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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wrc0(fname,n1,at1,xyz1)
      implicit none
      character*(*) fname
      character*2 asym
      integer n1,at1(n1)
      real*8 xyz1(3,n1)
      integer j

      open(unit=1,file=fname)
      write(1,'(''$coord'')')
      do j=1,n1
         write(1,'(3F24.12,5x,a2)') xyz1(1:3,j),asym(at1(j))
      enddo
      write(1,'(''$end'')')
      close(1)

      end
