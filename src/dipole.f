

      subroutine dipolecent(n,nbf,nao,nocc,xyz,C)
      implicit none
      integer n,nbf,nao,nocc

      real*8 C(nao,nao),xyz(3,n)
      real*8,allocatable :: dip1(:),dip2(:),dip3(:)
      real*8,allocatable :: x(:),y(:),z(:),scr(:,:)

      integer i,j,k,m
      real*8 x1,x2,x3

      allocate(dip1(nbf*(nbf+1)/2),
     .         dip2(nbf*(nbf+1)/2),
     .         dip3(nbf*(nbf+1)/2))

      call Dints(n,nbf,xyz,dip1,dip2,dip3)
      call cao2saop(nbf,nao,dip1)
      call cao2saop(nbf,nao,dip2)
      call cao2saop(nbf,nao,dip3)

      do i=1,nocc
         call OTRONE3(I,C,dip1,dip2,dip3,X1,X2,X3,NAO)
      enddo

      deallocate(dip1,dip2,dip3)
      end

C***********************************************************************
C symmetric trafo for three integrals, mos i,j
C***********************************************************************

      SUBROUTINE OTRONE3(I,C,A1,A2,A3,X1,X2,X3,NAO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NAO,NAO)
      DIMENSION A1(*),A2(*),A3(*)

      K=0
      X1=0.0D0
      X2=0.0D0
      X3=0.0D0
      DO 110 M=1,NAO
         NEND=M-1
         DO 120 N=1,NEND
            K=K+1
            CC=2.*C(m,i)*C(n,i)
            X1=X1+A1(k)*CC
            X2=X2+A2(k)*CC
            X3=X3+A3(k)*CC
120      CONTINUE
         K=K+1
         CC=C(m,i)*C(m,i)
         X1=X1+A1(k)*CC
         X2=X2+A2(k)*CC
         X3=X3+A3(k)*CC
110   CONTINUE

      RETURN
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dipole(n,nbf,nao,xyz,z,P,dip,d,pr)
      implicit none
      integer n,nbf,nao
      real*8 xyz(3,n),z(n),P(nao,nao),dip
      logical pr
      real*8 d(3)

      real*8,allocatable :: dip1(:),dip2(:),dip3(:)

      integer i,j,k,m

      allocate(dip1(nbf*(nbf+1)/2),
     .         dip2(nbf*(nbf+1)/2),
     .         dip3(nbf*(nbf+1)/2))

c core part
      d = 0
      do i=1,n
         d(1)=d(1)+xyz(1,i)*z(i)
         d(2)=d(2)+xyz(2,i)*z(i)
         d(3)=d(3)+xyz(3,i)*z(i)
      enddo

      call Dints(n,nbf,xyz,dip1,dip2,dip3)
c     call prmat(6,dip1,nbf,0,'X')
c     call prmat(6,dip2,nbf,0,'Y')
c     call prmat(6,dip3,nbf,0,'Z')
      call cao2saop(nbf,nao,dip1)
      call cao2saop(nbf,nao,dip2)
      call cao2saop(nbf,nao,dip3)

c contraction with P
      k=0
      do i=1,nao
         do j=1,i-1
            k=k+1
            d(1)=d(1)+2.0d0*P(j,i)*dip1(k)
            d(2)=d(2)+2.0d0*P(j,i)*dip2(k)
            d(3)=d(3)+2.0d0*P(j,i)*dip3(k)
         enddo
         k=k+1
         d(1)=d(1)+P(i,i)*dip1(k)
         d(2)=d(2)+P(i,i)*dip2(k)
         d(3)=d(3)+P(i,i)*dip3(k)
      enddo

      dip=sqrt(d(1)**2+d(2)**2+d(3)**2)

      if(pr)then
      write(*,*)
      write(*,*)'dipole moment from electron density (au)'
      write(*,*)'    X       Y       Z   '
      write(*,'(3f9.4,''  total (Debye): '',f8.3)')
     .     d(1),   d(2),   d(3), dip*2.5418
      write(*,*)
      endif

      deallocate(dip1,dip2,dip3)
      end



      subroutine Dints(n,nbf,xyz,S1,S2,S3)
      use intpack
      implicit none
      include 'ehtcommon.fh'
      integer n,nbf
      real*8  S1(nbf*(nbf+1)/2)
      real*8  S2(nbf*(nbf+1)/2)
      real*8  S3(nbf*(nbf+1)/2)
      real*8  xyz(3,n)

      integer i,j,k,l,iprimcount,jprimcount
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn
      real*8  xyza(3),xyzb(3),rab,est,ss,sss,lmnfak(84),gama,pi,arg
      real*8  point(3),gm2,ttt(3),tt1,tt2,tt3,intcut
      data pi/3.1415926535897932384626433832795029d0/

      intcut=20.

      point=0
      s1=0
      s2=0
      s3=0

      k=0
      iprimcount=0
      do i=1,nbf
c aufpunkt i
         xyza(1:3)=xyz(1:3,aoat(i))
c #prims
         npri=nprim(i)
         jprimcount=0
         do j=1,i
            k=k+1
            nprj=nprim(j)
c aufpunkt j
              xyzb(1:3)=xyz(1:3,aoat(j))
              rab=(xyza(1)-xyzb(1))**2
     .           +(xyza(2)-xyzb(2))**2
     .           +(xyza(3)-xyzb(3))**2
              if(rab.gt.200) goto 42    ! cut-off gives crambin dipole accurate to 1d-3 Deb
c prim loop
              tt1=0.0d0
              tt2=0.0d0
              tt3=0.0d0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    gama=1.0d0/(alp(iii)+alp(jjj))
                    est=rab*alp(iii)*alp(jjj)*gama
c cutoff
                    if(est.lt.intcut)then
                    ttt=0
                    call propa(opab1,xyza,xyzb,point,
     .                         alp(iii),alp(jjj),
     .                         lao(i),lao(j),ttt,3)
                    tt1=tt1+ttt(1)*cont(iii)*cont(jjj)
                    tt2=tt2+ttt(2)*cont(iii)*cont(jjj)
                    tt3=tt3+ttt(3)*cont(iii)*cont(jjj)
                    endif
                 enddo
              enddo
            s1(k)=tt1
            s2(k)=tt2
            s3(k)=tt3
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine DQ3ints(n,nbf,xyz,S1,S2,S3,Q1,Q2,Q3)
      use intpack
      implicit none
      include 'ehtcommon.fh'
      integer n,nbf
      real*8  S1(nbf*(nbf+1)/2)
      real*8  S2(nbf*(nbf+1)/2)
      real*8  S3(nbf*(nbf+1)/2)
      real*8  Q1(nbf*(nbf+1)/2)
      real*8  Q2(nbf*(nbf+1)/2)
      real*8  Q3(nbf*(nbf+1)/2)
      real*8  xyz(3,n)

      integer i,j,k,l,iprimcount,jprimcount
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn
      real*8  xyza(3),xyzb(3),rab,est,ss,sss,lmnfak(84),gama,pi,arg,cc
      real*8  point(3),gm2,ttt(3),tt1,tt2,tt3,intcut,qqq(6),qt1,qt2,qt3
      data pi/3.1415926535897932384626433832795029d0/

      intcut=20.

      point=0
      s1=0
      s2=0
      s3=0
      q1=0
      q2=0
      q3=0

      k=0
      iprimcount=0
      do i=1,nbf
c aufpunkt i
         xyza(1:3)=xyz(1:3,aoat(i))
c #prims
         npri=nprim(i)
         jprimcount=0
         do j=1,i
            k=k+1
            nprj=nprim(j)
c aufpunkt j
              xyzb(1:3)=xyz(1:3,aoat(j))
              rab=(xyza(1)-xyzb(1))**2
     .           +(xyza(2)-xyzb(2))**2
     .           +(xyza(3)-xyzb(3))**2
              if(rab.gt.200) goto 42    ! cut-off gives crambin dipole accurate to 1d-3 Deb
c prim loop
              tt1=0.0d0
              tt2=0.0d0
              tt3=0.0d0
              qt1=0.0d0
              qt2=0.0d0
              qt3=0.0d0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    gama=1.0d0/(alp(iii)+alp(jjj))
                    est=rab*alp(iii)*alp(jjj)*gama
c cutoff
                    if(est.lt.intcut)then
                    ttt=0
                    call propa(opab1,xyza,xyzb,point,
     .                         alp(iii),alp(jjj),
     .                         lao(i),lao(j),ttt,3)
                    call propa(opab4,xyza,xyzb,point,
     .                         alp(iii),alp(jjj),
     .                         lao(i),lao(j),qqq,6)
                    cc=cont(iii)*cont(jjj)
                    tt1=tt1+ttt(1)*cc
                    tt2=tt2+ttt(2)*cc
                    tt3=tt3+ttt(3)*cc
                    qt1=qt1+qqq(1)*cc
                    qt2=qt2+qqq(2)*cc
                    qt3=qt3+qqq(3)*cc
                    endif
                 enddo
              enddo
            s1(k)=tt1
            s2(k)=tt2
            s3(k)=tt3
            q1(k)=qt1
            q2(k)=qt2
            q3(k)=qt3
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine DQints(n,nbf,xyz,D,Q)
      use intpack
      implicit none
      include 'ehtcommon.fh'
      integer n,nbf
      real*8  D(3,nbf*(nbf+1)/2)
      real*8  Q(6,nbf*(nbf+1)/2)
      real*8  xyz(3,n)

      integer i,j,k,l,iprimcount,jprimcount,ip
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn
      real*8  xyza(3),xyzb(3),rab,est,lmnfak(84),gama,pi,arg
      real*8  point(3),gm2,ttt(6),tt(6),intcut,cc,ss(3),sss(3)
      data pi/3.1415926535897932384626433832795029d0/

      intcut=20.

      point=0

      k=0
      iprimcount=0
      do i=1,nbf
c aufpunkt i
         xyza(1:3)=xyz(1:3,aoat(i))
c #prims
         npri=nprim(i)
         jprimcount=0
         do j=1,i
            k=k+1
            nprj=nprim(j)
c aufpunkt j
              xyzb(1:3)=xyz(1:3,aoat(j))
              rab=(xyza(1)-xyzb(1))**2
     .           +(xyza(2)-xyzb(2))**2
     .           +(xyza(3)-xyzb(3))**2
c prim loop
              ss=0
              tt=0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    gama=1.0d0/(alp(iii)+alp(jjj))
                    est=rab*alp(iii)*alp(jjj)*gama
c cutoff
                    if(est.lt.intcut)then
                    ttt=0
                    sss=0
                    call propa(opab4,xyza,xyzb,point,
     .                         alp(iii),alp(jjj),
     .                         lao(i),lao(j),ttt,6)
                    call propa(opab1,xyza,xyzb,point,
     .                         alp(iii),alp(jjj),
     .                         lao(i),lao(j),sss,3)
                    cc=cont(iii)*cont(jjj)
                    do ip=1,3
                       ss(ip)=ss(ip)+sss(ip)*cc
                       tt(ip)=tt(ip)+ttt(ip)*cc
                    enddo
                    do ip=4,6
                       tt(ip)=tt(ip)+ttt(ip)*cc
                    enddo
                    endif
                 enddo
              enddo
            D(1:3,k)=-ss(1:3)
            Q(1:6,k)=-tt(1:6)
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

      end
