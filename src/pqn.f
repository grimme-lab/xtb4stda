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


c principal quantum number of valence shell

      integer function pqn(at)
      integer at

      if(at.le.2)then
         pqn=1
      elseif(at.le.10)then
         pqn=2
      elseif(at.le.18)then
         pqn=3
      elseif(at.le.36)then
         pqn=4
      elseif(at.le.54)then
         pqn=5
      else
         pqn=6
      endif

      end

      integer function ncore(at)
      integer at

      if(at.le.2)then
         ncore=0
      elseif(at.le.10)then
         ncore=2
      elseif(at.le.18)then
         ncore=10
      elseif(at.le.29)then   !zn
         ncore=18
      elseif(at.le.36)then
         ncore=28
      elseif(at.le.47)then
         ncore=36
      elseif(at.le.54)then
         ncore=46
      elseif(at.le.71)then
         ncore=54
      elseif(at.le.79)then
         ncore=68
      elseif(at.le.86)then
         ncore=78
      endif

      end
