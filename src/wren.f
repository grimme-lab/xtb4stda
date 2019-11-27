cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wren(e)
      implicit none
      integer j,nl,nn,idum
      real*8 xx(10)
      real*8 e,dum
      character*128 a
      logical ex

      inquire(file='energy',exist=ex)
      if(.not.ex)then
        dum = -99
        idum= 1
        open(unit=42,file='energy')
        write(42,'(''$energy'')')
        write(42,'(i3,4F24.11)')idum,e,dum,dum,dum
        write(42,'(''$end'')')
        close(42)
        return
      endif

c write file energy
      nl=0
      j=1
      open(unit=42,file='energy')
      open(unit=43,file='energy.tmp')
 50   read(42,'(a)',end=100)a
      call readl(a,xx,nn)
      if(nn.gt.3)nl=j
      j=j+1
      goto 50
100   continue

      rewind 42
      j=0
  20  read(42,'(a)',end=200)a
      j=j+1
      if(j.lt.nl)then
         write(43,'(a)')trim(a)
         call readl(a,xx,nn)
      else
         call readl(a,xx,nn)
         xx(2)=e
         xx(3:10)=-99
         write(43,'(i6,4F24.11)')idint(xx(1)),xx(2:nn)
         a='$end'
         write(43,'(a)')trim(a)
         goto 200
      endif
      goto 20
 200  continue
      close(42)
      close(43)

      call execute_command_line('mv energy.tmp energy')

      end

