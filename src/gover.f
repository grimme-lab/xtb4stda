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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computes s/t integrals in cao(6d) basis
c tm int version
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xstints(nat,nbf,xyz,s,t)
      use intpack
      implicit none
      integer nat,nbf
      real*8  xyz(3,nat)
      real*8  s(nbf,nbf)
      real*8  t(nbf*(nbf+1)/2)

      include 'ehtcommon.fh'

      integer i,j,k,l,iprimcount,jprimcount
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn,n
      integer lll(20),iall(4,4)
      data lll/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
      integer lin(84),min(84),nin(84)
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,
     8 1,1,4,2/
      integer lmnexp(84),ib,ie
      real*8  xyza(3),xyzb(3),rab,est,ss,sss,lmnfak(84),gama,pi,arg
      real*8  aa(20),bb(20),gm2,ttt(1),tt,intcut
      data pi/3.1415926535897932384626433832795029d0/
      external opap5

      intcut=25.

      s=0

      iall(1,1)=1
      iall(1,2)=4
      iall(2,1)=4
      iall(1,3)=10
      iall(3,1)=10
      iall(2,2)=10
      iall(2,3)=20
      iall(3,2)=20
      iall(1,4)=20
      iall(4,1)=20
      iall(2,4)=35
      iall(4,2)=35
      iall(3,3)=35
      iall(3,4)=56
      iall(4,3)=56
      iall(4,4)=84

      k=0
      iprimcount=0
      do i=1,nbf
c aufpunkt i
         xyza(1:3)=xyz(1:3,aoat(i))
c #prims
         npri=nprim(i)
         jprimcount=0
         li=lll(lao(i))
         do j=1,i
            lj=lll(lao(j))
            k=k+1
            nprj=nprim(j)
c aufpunkt j
              xyzb(1:3)=xyz(1:3,aoat(j))
              rab=(xyza(1)-xyzb(1))**2
     .           +(xyza(2)-xyzb(2))**2
     .           +(xyza(3)-xyzb(3))**2
c precalc some overlap terms that depend only on lm
              do ll=1,iall(li,lj)
              call lmnpre(lin(ll),min(ll),nin(ll),lmnexp(ll),lmnfak(ll))
              enddo
              aa=0
              bb=0
              aa(lao(i))=1.0d0
              bb(lao(j))=1.0d0
c prim loop
              ss=0.0d0
              tt=0.0d0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    gama=1.0d0/(alp(iii)+alp(jjj))
                    gm2 =0.5d0*gama
                    est=rab*alp(iii)*alp(jjj)*gama
c cutoff
                    if(est.lt.intcut)then
                    arg=(pi*gama)**1.50d0
                    call pola(xyza,xyzb,alp(iii),alp(jjj),
     .                        gama,gm2,lao(i),lao(j),iall(li,lj),
     .                        aa,bb,lmnexp,lmnfak,est,arg,sss)
                    ss=ss+sss   *cont(iii)*cont(jjj)
                    call prola(opap5,xyza,xyzb,
     .                         alp(iii),alp(jjj),
     .                         lao(i),lao(j),ttt,1)
                    tt=tt+ttt(1)*cont(iii)*cont(jjj)
                    endif
                 enddo
              enddo
            s(i,j)=ss
            s(j,i)=ss
            t(k)  =tt
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

c normalized?
      do i=1,nbf
c no diagonal contribution to H0
         if(abs(1.d0-1.0d0/sqrt(s(i,i))).gt.1.d-6)then
            write(*,*) i,s(i,i)
c           stop 'function not normalized inside stints'
         endif
      enddo

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computes s integrals in cao(6d) basis
c tm int version
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xsintsp(intcut,nat,nbf,xyz,s)
      use intpack
      implicit none
      integer nat,nbf
      real*8  xyz(3,nat)
      real*8  s(nbf*(nbf+1)/2)
      real*8  intcut

      include 'ehtcommon.fh'

      integer i,j,k,l,iprimcount,jprimcount
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn,n
      integer lll(20),iall(4,4)
      data lll/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
      integer lin(84),min(84),nin(84)
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,
     8 1,1,4,2/
      integer lmnexp(84),ib,ie
      real*8  xyza(3),xyzb(3),rab,est,ss,sss,lmnfak(84),gama,pi,arg
      real*8  aa(20),bb(20),gm2,ttt(1),tt
      data pi/3.1415926535897932384626433832795029d0/

      s=0

      iall(1,1)=1
      iall(1,2)=4
      iall(2,1)=4
      iall(1,3)=10
      iall(3,1)=10
      iall(2,2)=10
      iall(2,3)=20
      iall(3,2)=20
      iall(1,4)=20
      iall(4,1)=20
      iall(2,4)=35
      iall(4,2)=35
      iall(3,3)=35
      iall(3,4)=56
      iall(4,3)=56
      iall(4,4)=84

      k=0
      iprimcount=0
      do i=1,nbf
c aufpunkt i
         xyza(1:3)=xyz(1:3,aoat(i))
c #prims
         npri=nprim(i)
         jprimcount=0
         li=lll(lao(i))
         do j=1,i
            lj=lll(lao(j))
            k=k+1
            nprj=nprim(j)
c aufpunkt j
              xyzb(1:3)=xyz(1:3,aoat(j))
              rab=(xyza(1)-xyzb(1))**2
     .           +(xyza(2)-xyzb(2))**2
     .           +(xyza(3)-xyzb(3))**2
c precalc some overlap terms that depend only on lm
              do ll=1,iall(li,lj)
              call lmnpre(lin(ll),min(ll),nin(ll),lmnexp(ll),lmnfak(ll))
              enddo
              aa=0
              bb=0
              aa(lao(i))=1.0d0
              bb(lao(j))=1.0d0
c prim loop
              ss=0.0d0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    gama=1.0d0/(alp(iii)+alp(jjj))
                    gm2 =0.5d0*gama
                    est=rab*alp(iii)*alp(jjj)*gama
c cutoff
                    if(est.lt.intcut)then
                    arg=(pi*gama)**1.50d0
                    call pola(xyza,xyzb,alp(iii),alp(jjj),
     .                        gama,gm2,lao(i),lao(j),iall(li,lj),
     .                        aa,bb,lmnexp,lmnfak,est,arg,sss)
                    ss=ss+sss*cont(iii)*cont(jjj)
                    endif
                 enddo
              enddo
            s(k)=ss
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

c normalized?
      do i=1,nbf
c no diagonal contribution to H0
         ii=i+i*(i-1)/2
         if(abs(1.d0-1.0d0/sqrt(s(ii))).gt.1.d-6)then
            write(*,*) i,s(ii)
            stop 'function not normalized inside stints'
         endif
      enddo

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine opap5(ia,ib,ga,gb,ra,rb,l,m,n,gc,v,rc)
c
c     this subroutine calculates the expectation-values of
c     -p**2/2 (kinetic energy)
c     the program works correct only for s and p and the irreducible
c     linear combinations of the primitive d,f,... functions
c     written by s. brode in april,1984
c
      implicit real*8 (a-h,o-z)
      logical lodd,modd,nodd,leven,meven,neven
      dimension rc(3),rca(3),
     &lmni(35),o(25),dftr(7),ra(3),rb(3)
c ----- ra,rb centres ga,gb exponents ia,ib monom indices
      data a0/0.0d0/,a1/1.0d0/,a2/2.0d0/,a4/4.0d0/,a8/8.0d0/,
     &lmni/0,3*1,6*2,10*3,15*4/,
     &dftr/1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0,103935.0d0/
      data pi/3.1415926535897932384626433832795029d0/
c
c     define the overlap-integral (in line function definition)
c
      ola(ll,mm,nn)=gch**(ll/2+mm/2+nn/2)*dftr(ll/2+1)*dftr(mm/2+1)*
     &dftr(nn/2+1)
c
c     some previous calculations
c
      o(1:7)=0
      gch=a1/(a2*gc)
      ropigc=dsqrt(pi/gc)*pi/gc
      lodd=mod(l,2).eq.1
      leven=.not.lodd
      modd=mod(m,2).eq.1
      meven=.not.modd
      nodd=mod(n,2).eq.1
      neven=.not.nodd
      ga2=ga*a2
      prcaa=a0
      do 20 i=1,3
      rca(i)=rc(i)-ra(i)
      prcaa=prcaa+rca(i)*rca(i)
   20 continue
      zwa=ga2*prcaa-dble(2*lmni(ia)+3)
c
c     calculation of the overlap-integrals
c
      if(lodd.or.modd.or.nodd) go to 110
      o( 1)=ola(l   ,m  ,n  )
      o( 5)=ola(l+2,m   ,n  )
      o( 6)=ola(l   ,m+2,n  )
      o( 7)=ola(l   ,m  ,n+2)
  110 continue
c
      if(leven.or.modd.or.nodd) go to 120
      o( 2)=ola(l+1,m   ,n  )
  120 continue
c
      if(lodd.or.meven.or.nodd) go to 130
      o( 3)=ola(l   ,m+1,n  )
  130 continue
c
      if(lodd.or.modd.or.neven) go to 140
      o( 4)=ola(l   ,m  ,n+1)
  140 continue
c
      p2=-ga*(zwa*o(1)
     &       +ga2*(a2*(rca(1)*o(2)+rca(2)*o(3)+rca(3)*o(4))
     &       +o(5)+o(6)+o(7)))


      v=p2*ropigc
      return
      end

