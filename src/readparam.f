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

c testcode
c     implicit none
c     include 'aoelementcommon.fh'
c     character*128  atmp
c     integer i,j
c     open(unit=1,file='~/.paramx.xtb')
c     read(1,'(a)') head1
c     read(1,'(a)') head2
c     do j=1,86
c        read(1,'(i2,7F12.6,2x,a30)')
c    .   i,ao_lev(1:3,i),ao_exp(1:3,i),en(i),timestp(i)
c        ao_n(i)=3
c        ao_l(1,i)=0
c        ao_l(2,i)=1
c        ao_l(3,i)=2
c     enddo
c     close(1)
c     call rdelparam('paramv')
c     call wrelparam('test')
c     end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read global parameters (as string only) and element parameters +basis
c to common block
c basis set string layout:
c <pqn1><l1><pqn2><l2>....   e.g.
c 2s2p3d
c rules:
c 1. sp and spd denote Rydberg diffuse shells and they must follow
c    the valence shells (e.g. 3s3p4sp3d)
c 2. polarization exponents are taken as scaled of the larges valence l
c
c additionally sets default values for some parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rdelparam(ini,fname,globpar,gfnver)
      implicit none
      character*200 fname
      real*8 globpar(25)
      integer ini,gfnver
      include 'aoelementcommon.fh'
      real*8 readaa

c local
      integer ich
      integer iz
      integer nn
      integer i,j,k,l
      integer nao
      integer lmax
      real*8  xx(20),emax
      character*128 atmp,btmp,ctmp
      character*2 asym

      timestp='------------------------------'
      ao_pqn=0
      ao_l  =0
      ao_n  =0
      ao_lev=0
      ao_exp=0
      ao_typ=0
      polyr =0
      cxb   =0
      rep   =0
      mc    =0
      lpar  =0
      gam3  =0
      kcnat =0
      kqat  =0
      kpair =1

c     new GFN2 stuff AES related
      dpolc =0 ! read values are scaled by 0.01
      qpolc =0 !  "     "     "    "    "   "
      radaes=5 ! default atom radius
c     exceptions from default H,C-F,Si-Cl,Se
      radaes(1) =1.4
      radaes(2) =3.0
      radaes(6) =3.0
      radaes(7) =1.9
      radaes(8) =1.8
      radaes(9) =2.4
      radaes(14)=3.9
      radaes(15)=2.1
      radaes(16)=3.1
      radaes(17)=2.5
      radaes(34)=3.9
      radaes(35)=4.0

      if(ini.eq.1) globpar = 0.0d0

      ich=11
      open(unit=ich,file=fname)

 10   read(ich,'(a)',end=100) atmp

      if(index(atmp,'$globpar').ne.0)then
 12      read(ich,'(a)',end=100) atmp
         if(index(atmp,'ks').ne.0)     globpar(1) =readaa(atmp,1,i,j)
         if(index(atmp,'kp').ne.0)     globpar(2) =readaa(atmp,1,i,j)
         if(index(atmp,'kd ').ne.0)    globpar(3) =readaa(atmp,1,i,j)
         if(index(atmp,'kf').ne.0)     globpar(4) =readaa(atmp,1,i,j)
         if(index(atmp,'kdiffa').ne.0) globpar(5) =readaa(atmp,1,i,j)
         if(index(atmp,'kdiffb').ne.0) globpar(6) =readaa(atmp,1,i,j)
         if(index(atmp,'wllscal').ne.0)globpar(7) =readaa(atmp,1,i,j)
         if(index(atmp,'gscal'  ).ne.0)globpar(8) =readaa(atmp,1,i,j)
         if(index(atmp,'zcnf'   ).ne.0)globpar(9) =readaa(atmp,1,i,j)
         if(index(atmp,'tscal'  ).ne.0)globpar(10)=readaa(atmp,1,i,j)
         if(index(atmp,'kcn'    ).ne.0)globpar(11)=readaa(atmp,1,i,j)
         if(index(atmp,'fpol'   ).ne.0)globpar(12)=readaa(atmp,1,i,j)
         if(index(atmp,'ken'    ).ne.0)globpar(13)=readaa(atmp,1,i,j)
         if(index(atmp,'lshift ').ne.0)globpar(14)=readaa(atmp,1,i,j)
         if(index(atmp,'lshifta').ne.0)globpar(15)=readaa(atmp,1,i,j)
         if(index(atmp,'split'  ).ne.0)globpar(16)=readaa(atmp,1,i,j)
         if(index(atmp,'zqf'    ).ne.0)globpar(17)=readaa(atmp,1,i,j)
         if(index(atmp,'alphaj' ).ne.0)globpar(18)=readaa(atmp,1,i,j)
         if(index(atmp,'kexpo'  ).ne.0)globpar(19)=readaa(atmp,1,i,j)
         if(index(atmp,'dispa ' ).ne.0)globpar(20)=readaa(atmp,1,i,j)
         if(index(atmp,'dispb' ).ne.0) globpar(21)=readaa(atmp,1,i,j)
         if(index(atmp,'dispc' ).ne.0) globpar(22)=readaa(atmp,1,i,j)
         if(index(atmp,'dispatm').ne.0)globpar(23)=readaa(atmp,1,i,j)
         if(index(atmp,'xbdamp').ne.0) globpar(24)=readaa(atmp,1,i,j)
         if(index(atmp,'xbrad' ).ne.0) globpar(25)=readaa(atmp,1,i,j)
         if(index(atmp,'$end').eq.0)goto 12
      endif

      if(index(atmp,'$pairpar').ne.0)then
         write(*,*)
 13      read(ich,'(a)',end=100) atmp
         call readl(atmp,xx,nn)
         if(int(xx(1)).gt.0.and.int(xx(1)).le.86.and.
     .      int(xx(2)).gt.0.and.int(xx(1)).le.86.and.xx(3).gt.0)then
                kpair(int(xx(1)),int(xx(2)))=xx(3)
                kpair(int(xx(2)),int(xx(1)))=xx(3)
                write(*,'(''KAB for pair '',a2,'' - '',a2,'' :'',f8.4)')
     .          asym(int(xx(1))),asym(int(xx(2))),xx(3)
         endif
         if(index(atmp,'$end').eq.0)goto 13
      endif

      if(index(atmp,'Z=').ne.0)then
         call readl(atmp,xx,nn)

c        ordinal number
         iz=idint(xx(1))

c        timestamp
         timestp(iz)=atmp(7:35)

  20     read(ich,'(a)',end=100) atmp
         if(index(atmp,'$end').ne.0) goto 10

c        basis
         if(index(atmp,'ao=').ne.0)then
            btmp=atmp
            ctmp=atmp
            do i=1,len(atmp)
               if(atmp(i:i).eq.'=')j=i
            enddo
            j=j+1
            do i=j,len(atmp)
               if(btmp(i:i).eq.'s') btmp(i:i)=' '
               if(btmp(i:i).eq.'p') btmp(i:i)=' '
               if(btmp(i:i).eq.'d') btmp(i:i)=' '
               if(btmp(i:i).eq.'f') btmp(i:i)=' '
               if(ctmp(i:i).eq.'1') ctmp(i:i)=' '
               if(ctmp(i:i).eq.'2') ctmp(i:i)=' '
               if(ctmp(i:i).eq.'3') ctmp(i:i)=' '
               if(ctmp(i:i).eq.'4') ctmp(i:i)=' '
               if(ctmp(i:i).eq.'5') ctmp(i:i)=' '
               if(ctmp(i:i).eq.'6') ctmp(i:i)=' '
               if(ctmp(i:i).eq.'7') ctmp(i:i)=' '
               if(btmp(i:i).eq.'S') btmp(i:i)=' '
            enddo
            call readl(btmp,xx,ao_n(iz))
            k=j
            do i=1,ao_n(iz)
               ao_pqn(i,iz)=idint(xx(i))
               if(ctmp(k:k+1).eq.' s') ao_l(i,iz)=0
               if(ctmp(k:k+1).eq.' S') ao_l(i,iz)=11
c              if(ctmp(k:k+1).eq.' S') ao_l(i,iz)=0
               if(ctmp(k:k+1).eq.' p') ao_l(i,iz)=1
               if(ctmp(k:k+1).eq.' d') ao_l(i,iz)=2
               if(ctmp(k:k+1).eq.' f') ao_l(i,iz)=3
               if(ctmp(k:k+1).eq.' g') ao_l(i,iz)=4
               if(ctmp(k:k+3).eq.' spd')then
                                       ao_l(i,iz)=13
                                       k=k+2
                                       goto 30
                                       endif
               if(ctmp(k:k+2).eq.' sp')then
                                       ao_l(i,iz)=12
                                       k=k+1
                                       endif
  30           k=k+2
            enddo
         endif

c        level
         if(index(atmp,'lev=').ne.0)then
            call readl(atmp,xx,nn)
            if(nn.ne.ao_n(iz)) then
               write(*,*) nn,ao_n(iz),iz
               stop 'basis level read error'
            endif
            do i=1,nn
               ao_lev(i,iz)=xx(i)
            enddo
         endif

c        exponent
         if(index(atmp,'exp=').ne.0)then
            call readl(atmp,xx,nn)
            if(nn.ne.ao_n(iz)) then
               write(*,*) nn,ao_n(iz),iz
               stop 'basis exponent read error'
            endif
            do i=1,nn
               ao_exp(i,iz)=xx(i)
            enddo
c           polarization
            do i=1,nn
               if(ao_exp(i,iz).le.0) then
                 lmax=-1
                 do j=1,nn
                    if(ao_exp(j,iz).gt.0.and.ao_exp(j,iz).lt.4)then
                       if(ao_l(j,iz).gt.lmax)then
                          lmax=ao_l(j,iz)
                          emax=ao_exp(j,iz)
                       endif
                    endif
                 enddo
c              scaled valence exponent for polarization fct
               ao_exp(i,iz)=emax*globpar(12)
               ao_typ(i,iz)=-1
               endif
            enddo
         endif

c        EN
         if(index(atmp,'EN=').ne.0)then
            call readl(atmp,xx,nn)
            en(iz)=xx(1)
         endif
c        GAMMA
         if(index(atmp,'GAM=').ne.0)then
            call readl(atmp,xx,nn)
            gam(iz)=xx(1)
         endif
c        xi in Goed
         if(index(atmp,'XI=').ne.0)then  !NEWSG
            call readl(atmp,xx,nn)
            dpolc(iz)=xx(1)
         endif
c        alpha in Goed
         if(index(atmp,'ALPG=').ne.0)then
            call readl(atmp,xx,nn)
            radaes(iz)=xx(1)
         endif
c        third order
         if(index(atmp,'GAM3=').ne.0)then
            call readl(atmp,xx,nn)
            gam3(iz)=xx(2)*0.1
            if(gfnver.eq.0) gam3(iz)=xx(2)
         endif
c        XB potential
         if(index(atmp,'CXB=').ne.0)then
            call readl(atmp,xx,nn)
            cxb(iz)=xx(1)*0.1
            if(gfnver.eq.0) cxb(iz)=xx(1)
         endif
c        ONE CENTER AES DIPOLE SCAL
         if(index(atmp,'DPOL=').ne.0)then
            call readl(atmp,xx,nn)
            dpolc(iz)=xx(1)*0.01
            if(gfnver.eq.0) dpolc(iz)=xx(1)
         endif
c        ONE CENTER AES QPOLE SCAL
         if(index(atmp,'QPOL=').ne.0)then
            call readl(atmp,xx,nn)
            qpolc(iz)=xx(1)*0.01
            if(gfnver.eq.0) qpolc(iz)=xx(1)
         endif
c        REP
         if(index(atmp,'REPA=').ne.0)then
            call readl(atmp,xx,nn)
            rep(1,iz)=xx(nn)
         endif
         if(index(atmp,'REPB=').ne.0)then
            call readl(atmp,xx,nn)
            rep(2,iz)=xx(nn)
         endif
c        S*F(R)
         if(index(atmp,'POLYS=').ne.0)then
            call readl(atmp,xx,nn)
            polyr(1,iz)=xx(nn)
         endif
         if(index(atmp,'POLYP=').ne.0)then
            call readl(atmp,xx,nn)
            polyr(2,iz)=xx(nn)
         endif
         if(index(atmp,'POLYD=').ne.0)then
            call readl(atmp,xx,nn)
            polyr(3,iz)=xx(nn)
         endif
         if(index(atmp,'POLYF=').ne.0)then
            call readl(atmp,xx,nn)
            polyr(4,iz)=xx(nn)
         endif
         if(index(atmp,'LPARS=').ne.0)then
            call readl(atmp,xx,nn)
            lpar(0,iz)=xx(nn)*0.1
         endif
         if(index(atmp,'LPARP=').ne.0)then
            call readl(atmp,xx,nn)
            lpar(1,iz)=xx(nn)*0.1
         endif
         if(index(atmp,'LPARD=').ne.0)then
            call readl(atmp,xx,nn)
            lpar(2,iz)=xx(nn)*0.1
         endif
         if(index(atmp,'KCNS=').ne.0)then
            call readl(atmp,xx,nn)
            kcnat(0,iz)=xx(nn)*0.1
         endif
         if(index(atmp,'KCNP=').ne.0)then
            call readl(atmp,xx,nn)
            kcnat(1,iz)=xx(nn)*0.1
         endif
         if(index(atmp,'KCND=').ne.0)then
            call readl(atmp,xx,nn)
            kcnat(2,iz)=xx(nn)*0.1
         endif
         if(index(atmp,'KQS=').ne.0)then
            call readl(atmp,xx,nn)
            kqat(1,iz)=xx(nn)
         endif
         if(index(atmp,'KQP=').ne.0)then
            call readl(atmp,xx,nn)
            kqat(2,iz)=xx(nn)
         endif
         if(index(atmp,'KQD=').ne.0)then
            call readl(atmp,xx,nn)
            kqat(3,iz)=xx(nn)
         endif

         goto 20
      endif
      goto 10
100   continue
      close(ich)

      end

      logical function maingroup(i)
      integer i
      logical main_group(107)
      data main_group /
     &  2*.true.,                               ! H  - He
     &  8*.true.,                               ! Li - Ne
     &  8*.true.,                               ! Na - Ar
     &  2*.true., 9*.false., 7*.true.,          ! K  - Kr
     &  2*.true., 9*.false., 7*.true.,          ! Rb - Xe
     &  2*.true.,23*.false., 7*.true.,          ! Cs - Rn
     & 21*.true.                               / ! Fr - Tv

      maingroup = main_group(i)

      end

c global, predefined pair parameters
      subroutine setpair(gfnver)
      implicit none
      include 'aoelementcommon.fh'
      integer gfnver
      integer i,j,ii,jj
      integer tmmetal
      real*8  kp(3)
      logical notset

      if(gfnver.eq.1)then
      kp(1)=1.1    ! 3d
      kp(2)=1.2    ! 4d
      kp(3)=1.2    ! 5d or 4f
      write(*,'(''KAB for pair M(3d)-M(3d) :'',f8.4)')kp(1)
      write(*,'(''KAB for pair M(4d)-M(4d) :'',f8.4)')kp(2)
      write(*,'(''KAB for pair M(5d)-M(5d) :'',f8.4)')kp(3)
      elseif(gfnver.eq.0)then
      kp(1)=1.1
      kp(2)=1.1
      kp(3)=1.1
      write(*,'(''KAB for pair M(3d)-M(3d) :'',f8.4)')kp(1)
      write(*,'(''KAB for pair M(4d)-M(4d) :'',f8.4)')kp(2)
      write(*,'(''KAB for pair M(5d)-M(5d) :'',f8.4)')kp(3)
      elseif(gfnver.gt.1)then
      kp(1)=1.   ! 3d
      kp(2)=1.   ! 4d
      kp(3)=1.   ! 5d or 4f
c     write(*,'(''KAB for pair M(3d)-M(3d) :'',f8.4)')kp(1)
c     write(*,'(''KAB for pair M(4d)-M(4d) :'',f8.4)')kp(2)
c     write(*,'(''KAB for pair M(5d)-M(5d) :'',f8.4)')kp(3)
      endif

      do i=21,79
         do j=21,i
            ii=tmmetal(i)
            jj=tmmetal(j)
c           metal-metal interaction
            notset=abs(kpair(i,j)-1.0d0).lt.1.d-6 .and.
     .             abs(kpair(j,i)-1.0d0).lt.1.d-6
            if(ii.gt.0.and.jj.gt.0.and.notset) then
               kpair(i,j)=0.5*(kp(ii)+kp(jj))
               kpair(j,i)=0.5*(kp(ii)+kp(jj))
            endif
         enddo
      enddo

      end

      integer function tmmetal(i)
      integer i,j

      j=0
      if(i.gt.20.and.i.lt.30) j=1
      if(i.gt.38.and.i.lt.48) j=2
      if(i.gt.56.and.i.lt.80) j=3

      tmmetal=j

      end
