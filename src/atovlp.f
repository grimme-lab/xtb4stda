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

      subroutine atovlp(l,npri,nprj,alpa,alpb,conta,contb,ss)
      implicit none
      integer l,npri,nprj
      real*8 alpa(*),alpb(*)
      real*8 conta(*),contb(*)
      real*8 ss

      integer ii,jj
      real*8 ab,s00,sss,pi,ab05
      data pi/3.1415926535897932384626433832795029d0/

              SS=0.0d0
              do ii=1,npri
                 do jj=1,nprj
                    ab =1./(alpa(ii)+alpb(jj))
                    s00=(pi*ab)**1.50d00
                    if(l.eq.0)then
                       sss=s00
                    endif
                    if(l.eq.1)then
                       ab05=ab*0.5
                       sss=s00*ab05
                    endif
                    SS=SS+SSS*conta(ii)*contb(jj)
                 enddo
              enddo

      end
