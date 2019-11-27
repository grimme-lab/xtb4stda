c*****************************************************************
c transform to CMA for output on trj xyz file
c input coords remain unchanged
c*****************************************************************

      subroutine axis3(mode,numat,nat,coord,coordout,eig)
      implicit none
      integer numat,nat(numat),mode
      real*8  coord(3,numat),coordout(3,numat),eig(3)
      include 'atmass.fh'

      real*8 sumw,sumwx,sumwy,sumwz,atmass,xsum,eps
      real*8 t(6), evec(3,3)
      real*8 x(numat),y(numat),z(numat),coordtmp(3,numat)
      real*8 coord1(3,numat)
      integer i,j,k
      data t /6*0.d0/

      sumw=1.d-20
      sumwx=0.d0
      sumwy=0.d0
      sumwz=0.d0
      do 10 i=1,numat
         if(mode.eq.0)then
            atmass=ams(nat(i))
         else
            atmass=1./ams(nat(i))
         endif
         sumw=sumw+atmass
         sumwx=sumwx+atmass*coord(1,i)
         sumwy=sumwy+atmass*coord(2,i)
         sumwz=sumwz+atmass*coord(3,i)
   10 continue

      eps=1.d-3
      sumwx=sumwx/sumw
      sumwy=sumwy/sumw
      sumwz=sumwz/sumw
      do i=1,numat
         x(i)=coord(1,i)-sumwx
         y(i)=coord(2,i)-sumwy
         z(i)=coord(3,i)-sumwz
         coordtmp(1,i)=x(i)
         coordtmp(2,i)=y(i)
         coordtmp(3,i)=z(i)
      enddo

      do 40 i=1,6
   40 t(i)=dble(i)*1.0d-10

      do 50 i=1,numat
         atmass=ams(nat(i))
         t(1)=t(1)+atmass*(y(i)**2+z(i)**2)+eps
         t(2)=t(2)-atmass*x(i)*y(i)
         t(3)=t(3)+atmass*(z(i)**2+x(i)**2)+eps
         t(4)=t(4)-atmass*z(i)*x(i)
         t(5)=t(5)-atmass*y(i)*z(i)
         t(6)=t(6)+atmass*(x(i)**2+y(i)**2)+eps
   50 continue

      call rsp(t,3,3,eig,evec)

c   now to orient the molecule so the chirality is preserved
      xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) +
     1     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) +
     2     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))
      if( xsum .lt. 0) then
         do 80 j=1,3
   80    evec(j,1)=-evec(j,1)
      endif
c     call prmat(6,evec,3,3,'Rmat')

      do 120 i=1,numat
         do 120 j=1,3
            xsum=0.d0
            do 110 k=1,3
  110       xsum=xsum+coordtmp(k,i)*evec(k,j)
  120 coord1(j,i)=xsum

      do i=1,numat
         coordout(1,i)=coord1(1,i)
         coordout(2,i)=coord1(2,i)
         coordout(3,i)=coord1(3,i)
      enddo

      return
      end
