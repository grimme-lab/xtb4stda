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

      subroutine wrxyz(iu,n,at,xyz,e,e2)
      implicit none
      integer n,at(n),iu
      real*8 xyz(3,n),e,e2
      integer i,j
      character*2 asym
      real*8, parameter ::autoang=0.52917726d0

      write(iu,*) n
      write(iu,'('' SCF done '',2F16.8)') e,e2
      do j=1,n
         write(iu,'(a2,3F24.10)')asym(at(j)),xyz(1:3,j)*autoang
      enddo

      end
