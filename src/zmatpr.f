      subroutine zmatpr(nat,at,geo,na,nb,nc,molnum)
      implicit none
      integer nat,na(nat),nb(nat),nc(nat),at(nat)
      real*8 geo(3,nat)
      character*2 asym
      character*20 filename
      logical ex
      integer i,l,m,n,molnum
      real*8 bl,ang,dihed,pi
      parameter (pi =  3.14159265358979D0)
         do 250 i=1,nat
            l=1
            m=1
            n=1
            if(i.eq.1)then
               l=0
               m=0
               n=0
            endif
            if(i.eq.2)then
               m=0
               n=0
            endif
            if(i.eq.3)n=0
         bl=geo(1,i)
         ang=geo(2,i)*180./pi
         dihed=geo(3,i)*180./pi
         if(dihed.gt.180.0d0)dihed=dihed-360.0d0
c250     write(6,'(2x,a2,f14.8,i4,f14.8,i5,f14.8,i5,i6,2i5)')
c    .   asym(at(i)),bl,l,ang,m,dihed,n,na(i),nb(i),nc(i)
 250     write(6,'(i4,2x,a2,f12.6,2x,f10.4,2x,f10.4,i6,2i5)')
     .   i,asym(at(i)),bl,ang,dihed,na(i),nb(i),nc(i)

         write(filename,'("zmatrix",i0,".zmat")') molnum
         open(unit=42, file=filename)

         write(42,'(a2)') asym(at(1))
         write(42,'(a2,x,i0,x,f8.3)') asym(at(2)), na(2), geo(1,2)
         write(42,'(a2,x,i0,x,f8.3,x,i0,x,f8.3)') asym(at(3)), na(3)
     .                      ,geo(1,3),nb(3), geo(2,3)*180./pi

         do i=4,nat
            bl=geo(1,i)
            ang=geo(2,i)*180./pi
            dihed=geo(3,i)*180./pi
            if(dihed.gt.180.0d0)dihed=dihed-360.0d0
            write(42,'(a2,x,i0,x,f8.3,x,i0,x,f8.3,x,i0,x,f8.3)')
     .                asym(at(i)),na(i),bl,nb(i),ang,nc(i),dihed
         enddo
         write(42,*)
         close(42)
      end
