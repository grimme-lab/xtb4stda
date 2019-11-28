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

      subroutine outvip(eneut,ekat,shift)
      implicit none
      real*8 eneut,ekat,shift
      real*8 ip,dum
      logical ex

      write(*,*)
      write(*,*)'vertical deltaSCC IP calculation'
      ip=ekat-eneut-shift
      write(*,'(''empirical IP/EA shift (eV):'',f10.4)')
     .            27.21139570d0*shift
      write(*,'(''delta SCC IP (eV)'',f10.4)') 27.21139570d0*ip

c compare with reference for fit
      inquire(file='.ip',exist=ex)
      if(ex)then
         open(unit=11,file='.ip')
         read(11,*) dum
         close(11)
         open(unit=12,file='.IP')
         write(12,'(2F12.6)') dum,ip*27.21139570d0
         close(12)
      endif

      end

      subroutine outvea(eneut,eani,shift)
      implicit none
      real*8 eneut,eani,shift
      real*8 ea,dum
      logical ex

      write(*,*)
      write(*,*)'vertical deltaSCC EA calculation'
      ea=eneut-eani-shift
      write(*,'(''empirical IP/EA shift (eV):'',f10.4)')
     .            27.21139570d0*shift
      write(*,'(''delta SCC EA (eV):'',f10.4)') 27.21139570d0*ea

c compare with reference for fit
      inquire(file='.ea',exist=ex)
      if(ex)then
         open(unit=11,file='.ea')
         read(11,*) dum
         close(11)
         open(unit=12,file='.EA')
         write(12,'(2F12.6)') dum,ea*27.21139570d0
         close(12)
      endif

      end
