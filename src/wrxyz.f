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
