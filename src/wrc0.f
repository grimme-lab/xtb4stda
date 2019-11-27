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
