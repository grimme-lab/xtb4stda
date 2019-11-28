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
c setup TB basis
c runtyp > 1 means GFN0 or GFN1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xbasis(n,at,nbf,coord,q,cn,
     .           zqf,zcnf,split,ok,diff,runtyp,gfnver)
      implicit none
      integer elem,n,nbf,runtyp,gfnver
      integer at(n)
      logical ok,diff
      real*8  q(n),cn(n),zqf,split,zcnf
      real*8 coord(3,n)
      include 'ehtcommon.fh'
      include 'aoelementcommon.fh'

      integer i,j,m,l,ibf,ipr,p,thisprim,thisprimR,idum,npq,npqR,pqn
      real*8  a(10),c(10),zeta,k1,k2,split1,pp,zqfR,zcnfR,qi
      real*8  aR(10),cR(10),ss
      real*8  as(10),cs(10)
      real*8  ap(10),cp(10)

      if(n.gt.maxao/4) stop 'system too large. recompile code!'

      hdiag(1:nbf)=1.d+42
c note: Rydbergs are identified by valao(*)=0
c       polarization by            valao(*)=-1
      valao(1:nbf)=1

      ibf=0
      ipr=0
      ok=.true.

c     atom loop
      do i=1,n
c     charge only for non metals
      qi=0
      if(cn(i).gt.1.d-6) qi=q(i)
c     AO=shell loop
      fila(1,i)=ibf+1
      do m=1,ao_n(at(i))
c        principle QN
         npq=ao_pqn(m,at(i))
         l=ao_l(m,at(i))
cccccccccccccccccccccccccccccccccc
c        H-He
cccccccccccccccccccccccccccccccccc
         if(l.eq.0.and.at(i).le.2.and.npq.eq.1)then
c  s
            ibf =ibf+1
            zeta=ao_exp(m,at(i))*(1.0d0+zqf*qi)*(1.d0+0.1*zcnf*cn(i))
            if(runtyp.gt.1)then
               if(gfnver.eq.0)then !GFN0
               call setsto3(thisprim,npq,1,zeta,a,c)
               else
               call setsto4(thisprim,npq,1,zeta,a,c)
               endif
            else
               call setsto4(thisprim,npq,1,zeta,a,c)
            endif
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=c(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=1
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))
         endif

         if(l.eq.12.and.at(i).le.2)then
c diff sp
            ibf =ibf+1
            valao(ibf)=0
            zeta=ao_exp(m,at(i))
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,thisprim,thisprimR,a,aR,c,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               alp (ipr)=aR(p)
               cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=-ss*c(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=1
            nprim(ibf)=thisprimR+thisprim
            call atovlp(0,nprim(ibf),nprim(ibf),alp(idum),alp(idum),
     .                                        cont(idum),cont(idum),ss)
            do p=1,nprim(ibf)
               cont(idum-1+p)=cont(idum-1+p)/sqrt(ss)
            enddo

            split1=0
            if(at(i).eq.2)split1=2.25  ! set to SCS-CC2 value

            hdiag(ibf)=ao_lev(m,at(i))-split1
            zeta=ao_exp(m,at(i))
            call setsto4(thisprim,2,2,zeta,a,c)
            do j=2,4
            ibf=ibf+1
            valao(ibf)=0
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=c(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))+split1
            enddo
         endif

         if(l.eq.0.and.at(i).le.2.and.npq.eq.2)then
c diff s
            ibf =ibf+1
            valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(gfnver.eq.0)then !GFN0
               call setsto2(thisprimR,npq,1,zeta,aR,cR)
            else
               call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,a,aR,c,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               alp (ipr)=aR(p)
               cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=-ss*c(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=1
            nprim(ibf)=thisprimR+thisprim
            call atovlp(0,nprim(ibf),nprim(ibf),alp(idum),alp(idum),
     .                                        cont(idum),cont(idum),ss)
            do p=1,nprim(ibf)
               cont(idum-1+p)=cont(idum-1+p)/sqrt(ss)
            enddo
            hdiag(ibf)=ao_lev(m,at(i))
         endif

c p polarization
         if(l.eq.1.and.at(i).le.2)then
            zeta=ao_exp(m,at(i))*(1.0d0+zqf*qi)*(1.d0+0.1*zcnf*cn(i))
            call setsto3(thisprim,npq,2,zeta,ap,cp)
            do j=2,4
            ibf=ibf+1
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=ap(p)
               cont(ipr)=cp(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            valao(ibf)=-1
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

cccccccccccccccccccccccccccccccccc
c general sp
cccccccccccccccccccccccccccccccccc
         if(l.eq.0.and.at(i).gt.2)then
c   s
            ibf=ibf+1
            zeta=ao_exp(m,at(i))*(1.0d0+zqf*qi)*(1.d0+0.1*zcnf*cn(i))
            if(npq.gt.5) then
            call setsto6(thisprim,npq,1,zeta,as,cs)
            else
            if(runtyp.gt.1)then
               if(gfnver.eq.0)then !GFN0
               call setsto4(thisprim,npq,1,zeta,as,cs) !SGNEW
               else
               call setsto6(thisprim,npq,1,zeta,as,cs)
               endif
            else
               call setsto4(thisprim,npq,1,zeta,as,cs)
            endif
            endif
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=as(p)
               cont(ipr)=cs(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=1
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))
         endif
c   p
         if(l.eq.1.and.at(i).gt.2)then
            zeta=ao_exp(m,at(i))*(1.0d0+zqf*qi)*(1.d0+0.1*zcnf*cn(i))
            if(npq.gt.5) then
            call setsto6(thisprim,npq,2,zeta,ap,cp)
            else
            if(runtyp.gt.1)then
               if(gfnver.eq.0)then !GFN0
               call setsto3(thisprim,npq,2,zeta,ap,cp)
               else
               call setsto6(thisprim,npq,2,zeta,ap,cp)
               endif
            else
               call setsto4(thisprim,npq,2,zeta,ap,cp)
            endif
            endif
            do j=2,4
            ibf=ibf+1
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=ap(p)
               cont(ipr)=cp(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

cccccccccccccccccccccccccccccccccc
c diffuse. important: the code assumes that the previous shells are
c          valence on which the functions are orthogonalized
c          i.e. prims of val. are in as,ap, cs,cp
cccccccccccccccccccccccccccccccccc
         if(l.ge.12.and.at(i).gt.2)then
c   s
            ibf=ibf+1
            valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,as,aR,cs,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               alp (ipr)=aR(p)
               cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=as(p)
               cont(ipr)=-ss*cs(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=1
            nprim(ibf)=thisprimR+thisprim
            hdiag(ibf)=ao_lev(m,at(i))-split
            call atovlp(0,nprim(ibf),nprim(ibf),alp(idum),alp(idum),
     .                                        cont(idum),cont(idum),ss)
            do p=1,nprim(ibf)
               cont(idum-1+p)=cont(idum-1+p)/sqrt(ss)
            enddo
c  p
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprimR,npq,2,zeta,aR,cR)
            else
            call setsto3(thisprimR,npq,2,zeta,aR,cR)
            endif
            call atovlp(1,thisprim,thisprimR,ap,aR,cp,cR,pp)
            do j=2,4
            ibf=ibf+1
            valao(ibf)=0
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               alp (ipr)=aR(p)
               cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=ap(p)
               cont(ipr)=-pp*cp(p)
            enddo
            nprim(ibf)=thisprimR+thisprim
            call atovlp(1,nprim(ibf),nprim(ibf),alp(idum),alp(idum),
     .                                        cont(idum),cont(idum),ss)
            do p=1,nprim(ibf)
               cont(idum-1+p)=cont(idum-1+p)/sqrt(ss)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            hdiag(ibf)=ao_lev(m,at(i))+split
            enddo
c  d if its an spd Ryd shell
            if(l.eq.13)then
            zeta=ao_exp(m,at(i))
            if(npq.lt.5) then
            call setsto3(thisprim,npq,3,zeta,a,c)
            else
            call setsto4(thisprim,5  ,3,zeta,a,c)
            endif
            do j=5,10
            ibf=ibf+1
            valao(ibf)=0
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=c(p)
               if(j.gt.7)cont(ipr)=cont(ipr)*sqrt(3.0d0)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))+2.0d0*split
            enddo
            endif
         endif

c DZ s
         if(l.eq.11)then
            ibf=ibf+1
            valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,as,aR,cs,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               alp (ipr)=aR(p)
               cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=as(p)
               cont(ipr)=-ss*cs(p)
            enddo
            aoat (ibf)=i
            lao  (ibf)=1
            nprim(ibf)=thisprimR+thisprim
            hdiag(ibf)=ao_lev(m,at(i))
            call atovlp(0,nprim(ibf),nprim(ibf),alp(idum),alp(idum),
     .                                        cont(idum),cont(idum),ss)
            do p=1,nprim(ibf)
               cont(idum-1+p)=cont(idum-1+p)/sqrt(ss)
            enddo
         endif

cccccccccccccccccccccccccccccccccc
c d
cccccccccccccccccccccccccccccccccc
         if(l.eq.2)then
            zeta=ao_exp(m,at(i))*(1.0d0+zqf*qi)*(1.d0+0.1*zcnf*cn(i))
c           valence
            if(ao_typ(m,at(i)).ne.-1)then
            call setsto4(thisprim,npq,3,zeta,a,c)
            else
c           polarization
            call setsto4(thisprim,npq,3,zeta,a,c)
            endif
            do j=5,10
            ibf=ibf+1
            if(ao_typ(m,at(i)).eq.-1)valao(ibf)=-1
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=c(p)
               if(j.gt.7)cont(ipr)=cont(ipr)*sqrt(3.0d0)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

cccccccccccccccccccccccccccccccccc
c f
cccccccccccccccccccccccccccccccccc
         if(l.eq.3)then
            zeta=ao_exp(m,at(i))*(1.0d0+zqf*qi)*(1.d0+0.1*zcnf*cn(i))
            do j=11,20
            ibf=ibf+1
c           valence
            call setsto4(thisprim,npq,4,zeta,a,c)
            do p=1,thisprim
               ipr=ipr+1
               alp (ipr)=a(p)
               cont(ipr)=c(p)
               if(j.gt.13.and.j.lt.20)cont(ipr)=cont(ipr)*sqrt(5.0d0)
               if(j.eq.20            )cont(ipr)=cont(ipr)*sqrt(15.0d0)
            enddo
            aoat (ibf)=i
            lao  (ibf)=j
            nprim(ibf)=thisprim
            hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif
c next shell
         enddo
c next atom
         fila(2,i)=ibf
      enddo


      diff=.false.
      do i=1,ibf
         if(valao(i).eq.0) diff=.true.
         if(hdiag(i).gt.1.d+10)then
           write(*,*)'Hii not defined for',i,aoat(i)
           ok=.false.
         endif
      enddo
      do i=1,ipr
         if(alp(i).eq.0) then
          ok=.false.
          write(*,*)'alp=0 for',i
         endif
      enddo

      if(nbf.ne.ibf) then
      write(*,*) ibf,nbf
      stop 'internal error'
      endif

c     write(*,*) alp(1:ipr)
c     write(*,*) cont(1:ipr)

      do i=1,n
         do j=1,ao_n(at(i))
            l = ao_l(j,at(i))
            if(l.eq.11)ao_l(j,at(i))=0
         enddo
      enddo

      end

cccccccccccccccccccccccccccccccccccccccccc
c determine bf limits
cccccccccccccccccccccccccccccccccccccccccc

      subroutine xbasis0(n,at,nbf,nshell)
      implicit none
      integer n,nbf,nshell
      integer at(n)

      include 'ehtcommon.fh'
      include 'aoelementcommon.fh'

      integer i,j,k,l

      nbf=0
      nshell=0

      do i=1,n
         k=0
         do j=1,ao_n(at(i))
            l = ao_l(j,at(i))
            k = k + 1
            nshell=nshell+1
            if(l.eq.0) nbf=nbf+1
            if(l.eq.1) nbf=nbf+3
            if(l.eq.2) nbf=nbf+6
            if(l.eq.3) nbf=nbf+10
            if(l.eq.3) stop 'f-functions not implemented'
            if(l.eq.4) stop 'g-functions not implemented'
            if(l.eq.11)nbf=nbf+1
            if(l.eq.12)nbf=nbf+4
            if(l.eq.13)nbf=nbf+10
         enddo
         if(k.eq.0) then
             write(*,*) 'no basis found for atom', i,' Z= ',at(i)
             stop
         endif
      enddo

      if(nbf.gt.maxao) stop 'TB basis too large'

      end


