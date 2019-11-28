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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates in au or Ang. if its a xmol file
c redone by S.E. to avoid some input errors. Looks for $coord, ang, bohr or number (xmol) in the first line
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rdcoord(fname,n,xyz,iat,chrg)
      implicit none
      interface
        subroutine parse(str,delims,args,nargs)
          character(len=*),intent(inout) :: str
          character(len=*),intent(in)  :: delims
          character(len=*),dimension(:),intent(inout) :: args
          integer, intent(out) :: nargs
        end subroutine parse
      end interface

      real*8 xyz(3,*)
      integer iat(*),n,chrg
      character*(*) fname

      real*8 floats(3),f
      character*128 line
      character*80  strings(3)
      integer j,k,ich,cs,cf,ncheck

      if(index(fname,'.sdf').ne.0)then
         call rdsdf(fname,n,xyz,iat,chrg)
         return
      endif

      f=0.0d0
      ich=142
      open(unit=ich,file=fname)
      ncheck=0
      rewind(ich)
      DO
        read(ich,'(a)',end=200)line
        if(line.ne."") exit
      ENDDO

      call readline(line,floats,strings,cs,cf)
      if(cf.eq.1.and.floats(1).gt.0) then
         f=1./0.52917726d0
         read(ich,'(A)',end=200)line
      else if (index(line,'$coord').ne.0) then
         f=1.0d0
      else if (index(line,'ang').ne.0) then
         f=1./0.52917726d0
      else if (index(line,'bohr').ne.0) then
         f=1.0d0
      endif
      if(f.lt.1.0d0) then
       call stoprun('Coordinate format not recognized!')
      endif
      DO
         read(ich,'(a)',end=200)line
         if(index(line,'$redu').ne.0) exit
         if(index(line,'$user').ne.0) exit
         if(index(line,'$end' ).ne.0) exit
         if(index(line,'$set' ).ne.0) exit
         call readline(line,floats,strings,cs,cf)
         if(cf.ne.3) cycle
c        call readl(line,floats,k)
         call elem(strings(1),j)
         if(j.eq.0) then
          cycle !ignores dummies and unknown elements
         endif
         ncheck=ncheck+1
         xyz(1,ncheck)=floats(1)*f
         xyz(2,ncheck)=floats(2)*f
         xyz(3,ncheck)=floats(3)*f
         iat(ncheck)=j
      ENDDO

 200  continue

      if (n.ne.ncheck) then
          write(*,*)n,'/=',ncheck
          call stoprun('error reading coord file')
      endif
      close(ich)

c  321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
      end subroutine rdcoord


      subroutine rdatomnumber(fname,n)
      implicit none
      integer n
      character*(*) fname

      real*8 floats(3),f
      character*80 line
      character*80 strings(3)
      integer j,ich,cs,cf

      f=0.0d0
      ich=53
      open(unit=ich,file=fname)
      n=0
 300  read(ich,'(a)',end=200)line
      if(line.eq."") goto 300
      call readline(line,floats,strings,cs,cf)
      if(cf.eq.1.and.floats(1).gt.0.and.cs.eq.0) then
         f=1./0.52917726d0
!         write(*,*)floats(1)
         n=int(floats(1))
         close(ich)
         return
      else if (index(line,'$coord').ne.0) then
         f=1.0d0
      else if (index(line,'ang').ne.0) then
         f=1./0.52917726d0
      else if (index(line,'bohr').ne.0) then
         f=1.0d0
      endif
      if(f.lt.1.0d0) then
              write(*,*) f
       call stoprun('Coordinate format not recognized!')
      endif
      DO
         read(ich,'(a)',end=200)line
         if(index(line,'$redu').ne.0) exit
         if(index(line,'$user').ne.0) exit
         if(index(line,'$end' ).ne.0) exit
         call readline(line,floats,strings,cs,cf)
         if(cf.ne.3) exit
         call elem(strings(1),j)
         if(j.eq.0) cycle
         n=n+1
      ENDDO

 200  continue

      close(ich)

c  321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
      end subroutine rdatomnumber


!reads a line cuts the at blanks and tabstops and returns all floats and strings in order of occurence
      subroutine readline(line,floats,strings,cs,cf)
      implicit none
      real*8 floats(3)
      character*128 line
      character*80 strings(3)

      real*8 num
      character*80 stmp,str
      character*1 digit
      integer i,ty,cs,cf

      stmp=''
      cs=1
      cf=1
      strings(:)=''
      do i=1,len(trim(line))
       digit=line(i:i)
       if(digit.ne.' '.and.digit.ne.char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        stmp=trim(stmp)//trim(digit)
       elseif(stmp.ne.'')then
        call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
        if(ty.eq.0) then
         floats(cf)=num
         cf=cf+1
        elseif(ty.eq.1) then
         strings(cs)=str
         cs=cs+1
        else
          write(*,*)'Problem in checktype, must abort'
          exit
        endif
        stmp=''
       endif
       if(i.eq.len(trim(line))) then  !special case: end of line
        call checktype(stmp,num,str,ty)
        if(ty.eq.0) then
         floats(cf)=num
         cf=cf+1
        elseif(ty.eq.1) then
         strings(cs)=str
         cs=cs+1
        else
          write(*,*)'Problem in checktype, must abort'
          exit
        endif
        stmp=''
       endif
      enddo
      cs=cs-1
      cf=cf-1
      end subroutine readline


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine stoprun(s)
      character*(*) s
      write(*,*)'program stopped due to: ',s
      stop 'must stop!'
      end subroutine stoprun

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rd0(fname,n)
      implicit real*8 (a-h,o-z)
      dimension xx(10)
      character*128 line
      character*(*) fname
      logical ex

      n=0
      inquire(file=fname,exist=ex)
      if(.not.ex) return
      ich=142
      open(unit=ich,file=fname)

      if(index(fname,'.xyz').ne.0)then
      read(ich,*) n
      elseif(index(fname,'.sdf').ne.0)then
      read(ich,'(a)') line
      read(ich,'(a)') line
      read(ich,'(a)') line
      read(ich,'(i3)') n
      else
 100  read(ich,'(a)',end=200)line
         if(index(line,'$redu').ne.0)goto 200
         if(index(line,'$user').ne.0)goto 200
         if(index(line,'$end' ).ne.0)goto 200
         call readl(line,xx,nn)
         if(nn.ne.3) goto 100
         n=n+1
      goto 100
 200  continue
      endif

      if(n.gt.10000) stop 'too many (>10000) atoms'

      close(ich)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates, xyz only
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rdxyz(fname,n,xyz)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,n)
      character*128 line
      character*(*) fname
      real*8 xx(10)

      ich=142
      open(unit=ich,file=fname)

      if(index(fname,'xyz').ne.0)then
      read(ich,'(a)')line
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call readl(line,xx,nn)
         xyz(1:3,i)=xx(1:3)
      enddo
      xyz = xyz / 0.52917726
      else
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call readl(line,xx,nn)
         xyz(1:3,i)=xx(1:3)
      enddo
      endif

      close(ich)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read sdf coordinate file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rdsdf(fname,n,xyz,iat,chrg)
      implicit none
      include 'setcommon.fh'
      integer n, iat(n), chrg
      real*8 xyz(3,n)
      character*128 line,atmp
      character*(*) fname

      real*8 xx(20)
      integer nn,i,j,ich

      chrg=0
      ich=142
      open(unit=ich,file=fname)

      read(ich,'(a)')molnameline
      read(ich,'(a)')line
      read(ich,'(a)')commentline
      read(ich,'(i3)')n

      do i=1,n
         read(ich,'(a)')line
         call readl(line,xx,nn)
         call elem(line,j)
         iat(i) = j
         xyz(1:3,i)=xx(1:3)
      enddo

      xyz = xyz / 0.52917726

 100  read(ich,'(a)',end=200)line
         if(index(line,'$$$$').ne.0) goto 200
         if(index(line,'M  CHG').ne.0) then
            call readl(line,xx,nn)
            chrg=chrg+idint(xx(nn))
         endif
      goto 100
 200  continue

      close(ich)

      call execute_command_line('rm -f wbo') ! WBO are read for sdf output file

      end
