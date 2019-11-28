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

      subroutine wrcoord(iunit,n,xyz,iat,e,fname)
      use gbobc, only: lgbsa,lsalt,ionst,ion_rad
      implicit none
      include 'setcommon.fh'
      include 'fixcommon.fh'
      include 'scancommon.fh'
      include 'spherecommon.fh'
      include 'stuff.fh'

      integer n,iunit,iat(n)
      real*8 xyz(3,n),e
      character*(*) fname
      character*2 asym
      character*80 z1
      integer i,ii,nb,ierr,mo1(50000),nbond(3,3*n)
      real*8 bohr
      bohr=0.52917726d0

c     modflag=1
      mo1 = 0
      if(natomsf.gt.0)then
         do i=1,natomsf
            mo1(fixed(i))=1
         enddo
      endif

c take the input filename (and overwrite it!)
      if(index(fname,'none').ne.0)then
         fname = inputname
      endif

      if(iunit.ne.6) open(unit=iunit,file=fname)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc sdf
      if(index(inputname,'.sdf').ne.0)then
      call bondanal(n,iat,xyz,nb,nbond)
      write(iunit,'(a)') trim(molnameline)
      write(iunit,'('' xtb code, MCTC UBonn, Germany '')')
      write(iunit,'(a)') trim(commentline)
      write(iunit,'(2i3,''  0     0  0  0  0  0  0999 V2000'')')n,nb
      do i=1,n
         write(iunit,
     .   '(3F10.5,a3,''  0  0  0  0  0  0  0  0  0  0  0  0'')')
     .   xyz(1,i)*bohr,xyz(2,i)*bohr,xyz(3,i)*bohr,asym(iat(i))
      enddo
      do i=1,nb
         write(iunit,
     .   '(3i3,''  0  0  0  0'')')nbond(1:3,i)
      enddo
      write(iunit,'(''M  CHG  1   1'',i4)') ichrg
      write(iunit,'(''M  END'')')
      write(iunit,'(''> <total energy in Eh>'')')
      write(iunit,'(f18.8,/)') e
      if(gtot.ne.0)then
      write(iunit,'(''> <total free-energy in Eh>'')')
      write(iunit,'(f18.8,/)') gtot
      endif
      if(lgbsa)then
      write(iunit,'(''> <GBSA Gsolv in Eh>'')')
      write(iunit,'(f18.8,/)') gsolv
      endif
      write(iunit,'(''> <HL gap in eV>'')')
      write(iunit,'(f18.8,/)') hlgap
      write(iunit,'(''> <dipole moment in au>'')')
      write(iunit,'(f18.8,/)') dip
      write(iunit,'(''$$$$'')')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc xyz
      elseif(index(inputname,'.xyz').ne.0)then
      write(iunit,'(i6)')n
      write(iunit,'(F20.8)') e
      do i=1,n
         write(iunit,'(a2,6x,3F20.14)')asym(iat(i)),
     .   xyz(1,i)*bohr,xyz(2,i)*bohr,xyz(3,i)*bohr
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc coord with $set
      else

      write(iunit,'(a)')'$coord'
      do i=1,n
         write(iunit,'(3F24.14,6x,a2)')
     .   xyz(1,i),xyz(2,i),xyz(3,i),asym(iat(i))
      enddo
      write(iunit,'(a)')'$end'

      endif

      if(rdset) then
      write(iunit,'(a)')'$set'
      write(iunit,'(''chrg     '',i2)') ichrg
      write(iunit,'(''uhf      '',i2)') nalphabeta
      if(modflag(7) .eq.1)write(iunit,'('' etemp   '',f8.1)')etemp
      if(modflag(17).eq.1)write(iunit,'(''gbsa     '',a20)')solvent
      if(modflag(21).eq.1)then
                          write(iunit,'(''ion_st  '',f8.4)')ionst
                          write(iunit,'(''ion_rad '',f8.4)')ion_rad
      endif
      if(modflag(29).eq.1)write(iunit,'(''desy     '',f8.4)')desy
      write(iunit,'(a)')'$end'
      endif

      if(iunit.ne.6) close(unit=iunit)

      end

      subroutine bondanal(n,at,xyz,nb,nbond)
      implicit none
      integer n,at(n),nb,nbond(3,3*n)
      real*8 xyz(3,n)
      integer i,j
      real*8 w

      open(unit=22,file='wbo')
      nbond = 0
      nb=1
10    read(22,*,end=20) i,j,w
      nbond(1,nb)=i
      nbond(2,nb)=j
      if(w.lt.1.3)              nbond(3,nb)=1
      if(w.gt.1.3.and.w.lt.2.3) nbond(3,nb)=2
      if(w.gt.2.3             ) nbond(3,nb)=3
      if(w.gt.1.2.and.w.lt.1.5.and.at(i).eq.6.and.at(j).eq.6)
     .                          nbond(3,nb)=4
      nb = nb + 1
      goto 10
20    continue
      close(22)
      nb=nb-1
      end
