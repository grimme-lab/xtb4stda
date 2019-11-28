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

c aoat,lao,fila in CAO basis
c aoat2,lao2,fila2 in AO basis

      subroutine dtrafo2(n,nbf,nao)
      implicit none
      integer n,nbf,nao
      include 'ehtcommon.fh'
      integer i,j,be,en

c     do i=1,n
c        write(*,*) '6d ',fila(1:2,i)
c     enddo

      j=0
      do i=1,nbf
         if(lao(i).ne.5)then
            j=j+1
            lao2(j)=lao(i)
            if(lao(i).gt.4) lao2(j)=lao2(j)-1
            aoat2 (j)=aoat (i)
            valao2(j)=valao(i)
            hdiag2(j)=hdiag(i)
         endif
      enddo

      do j=1,n
         be=100000
         en=-1
         do i=1,nao
            if(aoat2(i).eq.j.and.i.lt.be)be=i
            if(aoat2(i).eq.j.and.i.gt.en)en=i
         enddo
         fila2(1,j)=be
         fila2(2,j)=en
      enddo

c     do i=1,n
c        write(*,*) '5d ',fila(1:2,i)
c     enddo
c     write(*,'(''6d '',40i3)') lao (1:nbf)
c     write(*,'(''5d '',40i3)') lao2(1:nbf)


      end
