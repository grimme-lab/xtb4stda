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

      subroutine warn(s)
      implicit none
      character*(*) s

      write(*,'(10x,65(''!''))')
      write(*,'(10x,65(''!''))')
      write(*,'(10x,''WARNING:'',a)')trim(s)
      write(*,'(10x,65(''!''))')
      write(*,'(10x,65(''!''))')

      end

      subroutine warn2(s1,s2)
      implicit none
      character*(*) s1,s2

      write(*,'(10x,45(''!''))')
      write(*,'(10x,45(''!''))')
      write(*,'(10x,a)')trim(s1)
      write(*,'(10x,a)')trim(s2)
      write(*,'(10x,45(''!''))')
      write(*,'(10x,45(''!''))')

      end
