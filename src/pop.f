cccccccccccccccccccccccccccccccccccccccccccccc
c density matrix
c C: MO coefficient
c X: scratch
c P  dmat
cccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dmat(ndim,focc,C,P)
      implicit none
      integer ndim
      real*8 focc(*)
      real*8 C(ndim,ndim)
      real*8 P(ndim,ndim)
      integer i,m
      real*8,allocatable ::Ptmp(:,:)

      allocate(Ptmp(ndim,ndim))
      do m=1,ndim
         do i=1,ndim
            Ptmp(i,m)=C(i,m)*focc(m)
         enddo
      enddo
      call DGEMM('N','T',ndim,ndim,ndim,1.0d0,C,
     .                   ndim,Ptmp,ndim,0.0d0,P,ndim)
      deallocate(Ptmp)

      end

cccccccccccccccccccccccccccccccccccccccccccccc
c            Wiberg BOs
cccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wiberg(n,ndim,at,xyz,P,S,wb,ex,fit)
      implicit none
      integer n,ndim,at(n)
      real*8 xyz(3,n)
      real*8 P(ndim,ndim)
      real*8 S(ndim,ndim)
      real*8 wb (n,n)
      logical ex,fit
      include 'ehtcommon.fh'

      real*8,allocatable ::Ptmp(:,:)
      real*8 xsum,rab
      real*8 wbr(n,n)
      real*8 ddum(3,n,n),cdum(n)
      integer i,j,k,m,ibmax,imem(n)
      character*2 asym

      allocate(Ptmp(ndim,ndim))
      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,P,
     .                   ndim,S,ndim,0.0d0,Ptmp,ndim)
      open(unit=43,file='wbo')
      wb=0
      do i=1,n
         do j=1,i-1
         xsum=0.0
         rab=(xyz(1,i)-xyz(1,j))**2
     .      +(xyz(2,i)-xyz(2,j))**2
     .      +(xyz(3,i)-xyz(3,j))**2
         if(rab.lt.100.0)then
            do k=fila2(1,i),fila2(2,i)   ! AOs on atom i
               do m=fila2(1,j),fila2(2,j) ! AOs on atom j
                  xsum=xsum+Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wb(i,j)=xsum
         wb(j,i)=xsum
         if(abs(xsum).gt.0.3) write(43,*)i,j,xsum
         enddo
      enddo
      close(43)
      deallocate(Ptmp)

      wbr = wb

c fit
c     if(ex)then
c        open(unit=43,file='wbo3')
c        do i=1,n
c           do j=1,i-1
c              if(abs(wbr(i,j)).gt.0.1.and.
c    .            abs( wb(i,j)).gt.0.1) write(43,*) wbr(j,i),wb(j,i)
c           enddo
c        enddo
c        close(43)
c     endif
c     if(fit) then
c        call geobo(n,at,xyz,cdum,ddum)
c        open(unit=112,file='wbofit')
c     endif

      write(*,*)
      write(*,*)'Wiberg/Mayer (AO) data. WBOs > 0.1 written file <wbo>'
      write(*,*)'largest (>0.05) Wiberg bond orders for each atom'
      write(*,*)'          total WBO             WBO to atom ...'
      do i=1,n
         do j=1,n
            imem(j)=j
         enddo
         call wibsort(n,i,imem,wbr)
         ibmax=0
         xsum =0
         do j=1,n
            if(wbr(i,j).gt.0.05)ibmax=j
            xsum=xsum+wbr(i,j)
         enddo
         write(*,'(i6,a4,1x,f6.3,9(4x,a2,i4,f6.3))')
     .   i,asym(at(i)),xsum,
     .   (asym(at(imem(j))),imem(j),wbr(i,j),j=1,ibmax)
c        if(fit) write(112,'(2F22.12)') xsum,cdum(i)
      enddo

c     if(fit) close(112)

      end

      subroutine wiberg_nosort(n,ndim,at,xyz,P,S,wb,ex)
      implicit none
      integer n,ndim,at(n)
      real*8 xyz(3,n)
      real*8 P(ndim,ndim)
      real*8 S(ndim,ndim)
      real*8 wb (n,n)
      logical ex
      include 'ehtcommon.fh'

      real*8,allocatable ::Ptmp(:,:)
      real*8 xsum,rab
      real*8 wbr(n,n)
      integer i,j,k,m,ibmax,imem(n)
      character*2 asym

      allocate(Ptmp(ndim,ndim))
      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,P,
     .                   ndim,S,ndim,0.0d0,Ptmp,ndim)
      open(unit=43,file='wbo')
      wb=0
      do i=1,n
         do j=1,i-1
         xsum=0.0
         rab=(xyz(1,i)-xyz(1,j))**2
     .      +(xyz(2,i)-xyz(2,j))**2
     .      +(xyz(3,i)-xyz(3,j))**2
         if(rab.lt.100.0)then
            do k=fila2(1,i),fila2(2,i)
               do m=fila2(1,j),fila2(2,j)
                  xsum=xsum+Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wb(i,j)=xsum
         wb(j,i)=xsum
         if(abs(xsum).gt.0.3) write(43,*)i,j,xsum
         enddo
      enddo
      close(43)
      deallocate(Ptmp)

      end

cccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE wibsort(ncent,imo,imem,qmo)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension qmo(ncent,ncent)
      dimension imem(ncent)

      do 140   ii = 2,ncent
         i = ii - 1
         k = i
         pp= qmo(imo,i)
         do 120   j = ii, ncent
            if (qmo(imo,j) .lt. pp) go to 120
            k = j
            pp=qmo(imo,j)
  120    continue
         if (k .eq. i) go to 140
         qmo(imo,k) = qmo(imo,i)
         qmo(imo,i) = pp

         ihilf=imem(i)
         imem(i)=imem(k)
         imem(k)=ihilf
  140 continue

      end

cccccccccccccccccccccccccccccccccccccccc
c Mulliken pop + AO pop
cccccccccccccccccccccccccccccccccccccccc

      subroutine mpopall(n,nao,aoat,S,P,qao,q)
      implicit none
      real*8  S (nao,nao)
      real*8  P (nao,nao)
      real*8  qao(nao),q(n),ps
      integer nao,n,aoat(nao)

      integer i,j,ii,jj,ij,is,js

      q  = 0
      qao= 0
      do i=1,nao
         ii=aoat(i)
         do j=1,i-1
            jj=aoat(j)
            ps=p(j,i)*s(j,i)
            q(ii)=q(ii)+ps
            q(jj)=q(jj)+ps
            qao(i)=qao(i)+ps
            qao(j)=qao(j)+ps
         enddo
         ps=p(i,i)*s(i,i)
         q(ii)=q(ii)+ps
         qao(i)=qao(i)+ps
      enddo

      end

cccccccccccccccccccccccccccccccccccccccc
c Mulliken pop
cccccccccccccccccccccccccccccccccccccccc

      subroutine mpop0(n,nao,aoat,S,P,q)
      implicit none
      real*8  S (nao,nao)
      real*8  P (nao,nao)
      real*8  q(n),ps
      integer nao,n,aoat(nao)

      integer i,j,ii,jj,ij,is,js

      q = 0
      do i=1,nao
         ii=aoat(i)
         do j=1,i-1
            jj=aoat(j)
            ps=p(j,i)*s(j,i)
            q(ii)=q(ii)+ps
            q(jj)=q(jj)+ps
         enddo
         ps=p(i,i)*s(i,i)
         q(ii)=q(ii)+ps
      enddo

      end

cccccccccccccccccccccccccccccccccccccccc
c Mulliken AO pop
cccccccccccccccccccccccccccccccccccccccc

      subroutine mpopao(n,nao,S,P,qao)
      implicit none
      real*8  S (nao,nao)
      real*8  P (nao,nao)
      real*8  qao(nao),ps
      integer nao,n

      integer i,j

      qao = 0
      do i=1,nao
         do j=1,i-1
            ps=p(j,i)*s(j,i)
            qao(i)=qao(i)+ps
            qao(j)=qao(j)+ps
         enddo
         ps=p(i,i)*s(i,i)
         qao(i)=qao(i)+ps
      enddo

      end

cccccccccccccccccccccccccccccccccccccccc
c Mulliken pop
cccccccccccccccccccccccccccccccccccccccc

      subroutine mpop(n,nao,aoat,lao,S,P,q,ql)
      implicit none
      real*8  S (nao,nao)
      real*8  P (nao,nao)
      real*8  q(n),ps
      real*8  ql(3,n)
      integer nao,n,aoat(nao),lao(nao)

      integer i,j,ii,jj,ij,is,js,mmm(20)
      data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/

      ql= 0
      q = 0
      do i=1,nao
         ii=aoat(i)
         is=mmm(lao(i))
         do j=1,i-1
            jj=aoat(j)
            js=mmm(lao(j))
            ps=p(j,i)*s(j,i)
            q(ii)=q(ii)+ps
            q(jj)=q(jj)+ps
            ql(is,ii)=ql(is,ii)+ps
            ql(js,jj)=ql(js,jj)+ps
         enddo
         ps=p(i,i)*s(i,i)
         q(ii)=q(ii)+ps
         ql(is,ii)=ql(is,ii)+ps
      enddo

      end

cccccccccccccccccccccccccccccccccccccccc
c Mulliken pop shell wise
cccccccccccccccccccccccccccccccccccccccc

      subroutine mpopsh(n,nao,nshell,ao2sh,S,P,qsh)
      implicit none
      real*8  S (nao,nao)
      real*8  P (nao,nao)
      real*8  qsh(nshell),ps
      integer nao,n,nshell,ao2sh(nao)

      integer i,j,ii,jj,ij

      qsh=0
      do i=1,nao
         ii =ao2sh(i)
         do j=1,i-1
            jj =ao2sh(j)
            ps=p(j,i)*s(j,i)
            qsh(ii)=qsh(ii)+ps
            qsh(jj)=qsh(jj)+ps
         enddo
         ps=p(i,i)*s(i,i)
         qsh(ii)=qsh(ii)+ps
      enddo

      end

      subroutine qsh2qat(n,at,nshell,qsh,q)
      implicit none
      integer n,nshell,at(n)
      real*8  qsh(nshell),q(n)
      include 'aoelementcommon.fh'

      integer i,mi

      nshell=0
      do i=1,n
         q(i)=0
         do mi=1,ao_n(at(i))
            nshell=nshell+1
            q(i)=q(i)+qsh(nshell)
         enddo
      enddo

      end


cccccccccccccccccccccccccccccccccccccccc
c Loewdin pop
cccccccccccccccccccccccccccccccccccccccc

      subroutine lpop(n,nao,aoat,lao,occ,C,f,q,ql)
      implicit none
      real*8  C (nao,nao)
      real*8  occ(nao)
      real*8  q(n)
      real*8  ql(3,n)
      real*8  f
      integer nao,n,aoat(nao),lao(nao)

      integer i,j,ii,jj,js,mmm(20)
      data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
      real*8  cc

      do i=1,nao
         if(occ(i).lt.1.d-8) cycle
         do j=1,nao
            cc=f*C(j,i)*C(j,i)*occ(i)
            jj=aoat(j)
            js=mmm(lao(j))
            q(jj)=q(jj)+cc
            ql(js,jj)=ql(js,jj)+cc
         enddo
      enddo

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c atomic valence shell pops and total atomic energy
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setzshell(n,at,z,nshell,zshell,lshell,ashell,
     .                     q,qsh,ini,e,gfnver)
      implicit none
      include 'aoelementcommon.fh'
      include 'ehtcommon.fh'
      integer n,nshell,at(n),lshell(nshell),ashell(nshell)
      integer gfnver
      real*8 zshell(nshell),z(n),qsh(nshell),q(n),e
      logical ini

      real*8  ntot,fracz
      real*8  iox(86,0:2,2) ! GFN1 on 1
      integer i,j,k,m,l,ll(0:3),iat,lll,iver
      data ll /1,3,5,7/

!     H           Initial     Orbital Occupancies                     He
!     Li Be                                            B  C  N  O  F  Ne
!     Na Mg                                            Al Si P  S  Cl Ar
!     K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!     Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!     Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!                                      spd shell

c GFN1
      data iox /
     &1,                                                          2, !He
     &1,2,                                         2, 2, 2, 2, 2, 2, !Ne
     &1,2,                                         2, 2, 2, 2, 2, 2, !Ar
     &1,2,2,             2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, !Kr
     &1,2,2,             2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, !Xe
     &1,2,2, 14*2       ,2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, !Rn

     &0,                                                          0, !He
     &0,0,                                         1, 2, 3, 4, 5, 6, !Ne
     &0,0,                                         1, 2, 3, 4, 5, 6, !Ar
     &0,0,0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, !Kr
     &0,0,0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, !Xe
     &0,0,0,  14*0,   0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, !Rn

     & 0,                                                         0, !He
     & 0,0,                                        0, 0, 0, 0, 0, 0, !Ne
     & 0,0,                                        0, 0, 0, 0, 0, 0, !Ar
     & 0,0,1,          2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, !Kr
     & 0,0,1,          2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, !Xe
     & 0,0,1, 14*1,    2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, !Rn
c GFN2
     &1,                                                            2, !He
     &1,2,                                        2  ,1  ,1.5,  2,2,2, !Ne
     &1,2,                                        2  ,1.5,1.5,  2,2,2, !Ar
     &1,2,1,             1, 1, 1, 1, 1, 1, 1,1,2, 2,  1.5,1.5,  2,2,2, !Kr
     &1,2,1,             1, 1, 1, 1, 1, 1, 1,1,2, 2,  2,  2 ,   2,2,2, !Xe
     &1,2,1, 14*1       ,1, 1, 1, 1, 1, 1, 1,1,2, 2,  2,  2 ,   2,2,2, !Rn

     &0,                                                            0, !He
     &0,0,                                        1  , 3  ,3.5, 4,5,6, !Ne
     &0,0,                                        1  , 2.5,3.5, 4,5,6, !Ar
     &0,0,1,            1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2.5,3.5, 4,5,6, !Kr
     &0,0,1,            1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2,  3,   4,5,6, !Xe
     &0,0,1,  14*1,     1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2,  3,   4,5,6, !Rn

     & 0,                                                          0, !He
     & 0,0,                                        0, 0 , 0, 0, 0, 0, !Ne
     & 0,0,                                        0, 0 , 0, 0, 0, 0, !Ar
     & 0,0,1,          2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0, !Kr
     & 0,0,1,          2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0, !Xe
     & 0,0,1, 14*1,    2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0  !Rn
     & /

      iver=1
      if(gfnver.gt.1) iver = 2

      k=0
      e=0
      do i=1,n
         iat=at(i)
         ntot=-1.d-6
         do m=1,ao_n(iat)
            l=ao_l(m,iat)
            k=k+1
            zshell(k)=iox(iat,l,iver)
            lshell(k)=l
            ashell(k)=i
            ntot=ntot+zshell(k)
            if(ntot.gt.z(i)) zshell(k)=0
            fracz=zshell(k)/z(i)
            if(ini) qsh(k)=fracz*q(i)
            e=e+ao_lev(m,iat)*zshell(k)
         enddo
      enddo
      if(k.ne.nshell) stop 'internal setzshell error 1'
      if(abs(sum(z)-sum(zshell)).gt.1.d-4) then
         write(*,*) i,sum(z),sum(zshell)
         stop 'internal setzshell error 2'
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc
