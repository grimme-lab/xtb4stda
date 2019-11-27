      subroutine splitm(n,ic,xyz,cn)
      implicit none
      include 'aoelementcommon.fh'
      include 'splitcommon.fh'
      integer n,ic(n)
      real*8 xyz(3,n)
      real*8 cn(n)

      integer i,j,k
      real*8 r,f,rco,bond(n,n)

c determine covalent neighbours
      bond = 0
      do i=1,n
         do j=1,n
            r=sqrt((xyz(1,i)-xyz(1,j))**2+
     .             (xyz(2,i)-xyz(2,j))**2+
     .             (xyz(3,i)-xyz(3,j))**2)
            rco=rad(ic(i))+rad(ic(j))
            if(r.lt.2.5*rco) bond(j,i)=1
         enddo
         bond(i,i)=0
      enddo

      call mrec(i,xyz,cn,bond,n,ic,splitlist)

      iatf1=0
      iatf2=0
      do i=1,n
         if(splitlist(i).gt.2) splitlist(i)=2
         if(splitlist(i).eq.1) iatf1=iatf1+1
         if(splitlist(i).eq.2) iatf2=iatf2+1
      enddo

      end

      subroutine cmafrag(n,ic,xyz,r1,r2)
      implicit none
      include 'aoelementcommon.fh'
      include 'splitcommon.fh'
      real*8 xyz(3,n),r1(3),r2(3)
      integer n,ic(n)

      integer i
      real*8 sum1x,atmas,sum1y,sum1z,sum2x,sum2y,sum2z,dr(3)

      sum1x=0.d0
      sum1y=0.d0
      sum1z=0.d0
      sum2x=0.d0
      sum2y=0.d0
      sum2z=0.d0
      do i=1,n
         atmas=atmass(i)
         if(splitlist(i).eq.1)then
         sum1x=sum1x+atmas*xyz(1,i)
         sum1y=sum1y+atmas*xyz(2,i)
         sum1z=sum1z+atmas*xyz(3,i)
         else
         sum2x=sum2x+atmas*xyz(1,i)
         sum2y=sum2y+atmas*xyz(2,i)
         sum2z=sum2z+atmas*xyz(3,i)
         endif
      enddo
c
      r1(1)=sum1x/massf1
      r1(2)=sum1y/massf1
      r1(3)=sum1z/massf1
      r2(1)=sum2x/massf2
      r2(2)=sum2y/massf2
      r2(3)=sum2z/massf2

      dr(1:3)=r1(1:3)-r2(1:3)

      rcma=sqrt(sum(dr*dr))

      end

      subroutine splitprint(n,ic,xyz)
      implicit none
      include 'splitcommon.fh'
      include 'scancommon.fh'
      real*8 xyz(3,n)
      integer n,ic(n)

      real*8 ra(3),rb(3)
      integer i

      if(iatf1.eq.0.or.iatf2.eq.0) return

      massf1=0
      massf2=0
      do i=1,n
         if(splitlist(i).eq.1)then
            massf1=massf1+atmass(i)
         else
            massf2=massf2+atmass(i)
         endif
      enddo
      call cmafrag(n,ic,xyz,ra,rb)

      write(*,*)
      write(*,*) 'molecular fragmentation (1/2 indicates fragments):'
      write(*,'(80i1)')splitlist(1:n)
      write(*,'(''# atoms in fragment 1/2:'',2i6)')iatf1,iatf2
      write(*,'('' fragment masses (1/2) :'',2f12.2)')massf1,massf2
      write(*,'(''CMA distance (Bohr)    :'',f8.3)')rcma
      write(*,'(''constraining FC (au)   :'',f8.4)')fcconstr

      end

      subroutine cmaiface(n,at,xyz)
      implicit none
      include 'splitcommon.fh'
      real*8 xyz(3,n)
      integer n,at(n)

      integer i,j
      real*8 r,r0,rvdw(94),asum,bsum,xsum,ff
      rvdw(1:94)= (/
     .1.09155,0.86735,1.7478 ,1.5491 ,1.608  ,1.45515,1.31125,1.24085,
     .1.1498 ,1.0687 ,1.8541 ,1.74195,2.0053 ,1.89585,1.75085,1.65535,
     .1.5523 ,1.4574 ,2.12055,2.05175,1.94515,1.8821 ,1.86055,1.7207,
     .1.7731 ,1.72105,1.71635,1.6731 ,1.6504 ,1.61545,1.97895,1.93095,
     .1.83125,1.7634 ,1.6831 ,1.6048 ,2.3088 ,2.2382 ,2.1098 ,2.02985,
     .1.9298 ,1.87715,1.7845 ,1.73115,1.69875,1.67625,1.6654 ,1.731,
     .2.13115,2.0937 ,2.0075 ,1.94505,1.869  ,1.79445,2.52835,2.5907,
     .2.31305,2.31005,2.2851 ,2.26355,2.2448 ,2.22575,2.2117 ,2.06215,
     .2.12135,2.07705,2.1397 ,2.1225 ,2.1104 ,2.0993 ,2.0065 ,2.1225,
     .2.049  ,1.99275,1.94775,1.8745 ,1.7228 ,1.67625,1.6282 ,1.67995,
     .2.15635,2.1382 ,2.05875,2.0027 ,1.9322 ,1.8608 ,2.5398 ,2.4647,
     .2.35215,2.2126 ,2.2297 ,2.19785,2.17695,2.21705/)

      atmass=0
      ff=-5.d0

      asum=0
      do i=1,n
         if(splitlist(i).eq.1)then
            xsum=0
            do j=1,n
               if(splitlist(j).ne.1)then
               r=sqrt((xyz(1,i)-xyz(1,j))**2+
     .                (xyz(2,i)-xyz(2,j))**2+
     .                (xyz(3,i)-xyz(3,j))**2)
               r0=rvdw(at(i))+rvdw(at(j))
               xsum=xsum+exp(ff*r/r0**2)
               endif
            enddo
            asum=asum+xsum
            atmass(i)=xsum
         endif
      enddo
      do i=1,n
         if(splitlist(i).eq.1)atmass(i)=atmass(i)/asum
      enddo

      bsum=0
      do i=1,n
         if(splitlist(i).ne.1)then
            xsum=0
            do j=1,n
               if(splitlist(j).eq.1)then
               r=sqrt((xyz(1,i)-xyz(1,j))**2+
     .                (xyz(2,i)-xyz(2,j))**2+
     .                (xyz(3,i)-xyz(3,j))**2)
               r0=rvdw(at(i))+rvdw(at(j))
               xsum=xsum+exp(ff*r/r0**2)
               endif
            enddo
            bsum=bsum+xsum
            atmass(i)=xsum
         endif
      enddo
      do i=1,n
         if(splitlist(i).ne.1)atmass(i)=atmass(i)/bsum
      enddo
      ff=maxval(atmass(1:n))
      atmass(1:n)=atmass(1:n)/ff

c     do i=1,n
c        write(*,*) i,atmass(i)
c     enddo

      end
