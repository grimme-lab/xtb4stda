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
c non-iterative, distance dependent Gasteiger-Grimme type charges
c in addition computes CN
c n: number of atoms
c nel: number of electrons
c at: ordinal numbers of atoms
c xyz: coordinates in Bohr
c z: nuclear charges
c en: atomic EN
c q: partial charges (output)
c cn: D3 CN (output)
c kchrg1: parameter
c requires block.f as well as a call to setmetal.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine iniqcn(n,nel,at,z,xyz,chrg,q,cn,kchrg1,runtyp,gfnver)
      implicit none
      include 'aoelementcommon.fh'
      include 'd3common.fh'
      integer n,i,j,nel,chrg,runtyp,gfnver
      integer at(n)
      real*8  kchrg1
      real*8  z(n)
      real*8  q(n)
      real*8 cn(n)
      real*8  xyz(3,n)

      real*8 rij,dum,norm,r0,d(3),den,den2,ena(n),rexp
      real*8 rav,xij,yij,zij,q2(n)
      integer lin,kk,k,idum

      call setrcov(rcov)

      if(n.eq.1)then
         q(1)=float(chrg)
         cn(1)=0
         return
      endif

      write(*,*)
      write(*,*) 'doing EN charges ...'

      call ncoord2(.false.,n,at,xyz,cn,4000.0d0) ! CN needed for printout

c the effective CN corrected ENs
      do i=1,n
         ena(i)=en(at(i))
         if(runtyp.gt.1.and.metal(at(i)).eq.1) ena(i)=0
         if(gfnver.gt.1.and.runtyp.gt.1) then
            ena(i)=ena(i)+0.2*cn(i)**0.5 ! avoids too big start q
         else
            ena(i)=ena(i)-kchrg1*cn(i)**0.5
         endif
      enddo

c the neutral part
      q = z

      do i=1,n
         dum=0.0d0
         do j=1,n
            if(i.eq.j) cycle
            xij=xyz(1,i)-xyz(1,j)
            yij=xyz(2,i)-xyz(2,j)
            zij=xyz(3,i)-xyz(3,j)
            rav=0.5*(rcov(at(i))+rcov(at(j)))
            rij=sqrt(xij**2+yij**2+zij**2)/rav
            den =ena(i)-ena(j)
            rexp=1.0d0/rij**6
c           rexp=exp(-0.95*rij**2)   ! SIMILAR THAN 1/R^6
            dum=dum+den*rexp
         enddo
         q(i)=q(i)+dum
      enddo

c normalize to Nel
      q = q*(sum(z)-float(chrg))/sum(z)

c partial charges
      q = z - q


c distribute the extra molecular charge seperately
c     if(abs(chrg).gt.1.d-6)then
c     q2 = chrg/float(n)
c     do i=1,n
c        dum=0.0d0
c        do j=1,n
c           if(i.eq.j) cycle
c           xij=xyz(1,i)-xyz(1,j)
c           yij=xyz(2,i)-xyz(2,j)
c           zij=xyz(3,i)-xyz(3,j)
c           rav=0.5*(rcov(at(i))+rcov(at(j)))
c           rij=sqrt(xij**2+yij**2+zij**2)/rav
c           den = ena(i)-ena(j)
c           den =den+10.0*kchrg1*(gam(at(i))-gam(at(j)))**3
c           rexp=1.0d0/(rij**kchrg2+1.0d0)
c           dum=dum+den*rexp
c        enddo
c        q2(i)=q2(i)-dum
c     enddo
c     q = q2 + q
c     endif

c dipole moment of PCs
      d = 0
      do i=1,n
         d(1)=d(1)+xyz(1,i)*q(i)
         d(2)=d(2)+xyz(2,i)*q(i)
         d(3)=d(3)+xyz(3,i)*q(i)
      enddo
      write(*,'('' sum q : '',d14.7)')sum(q)
      write(*,*)'point charge moment (au)'
      write(*,*)'    X       Y       Z   '
      write(*,'(3f9.4,''  total (Debye): '',f8.3)')
     .     d(1),   d(2),   d(3),
     .sqrt(d(1)**2+d(2)**2+d(3)**2)*2.5418

c fitting stuff (EN parameters)
c     open(unit=43,file='~/.fitpar.xtb')
c     read(43,*) idum
c     read(43,'(20i1)') idum
c     if(idum.eq.2)then
c     open(unit=22,file='charges')
c     do i=1,n
c        read(22,*) ena(i)
c     enddo
c     close(22)
c     open(unit=43,file='charges2')
c     do i=1,n
c        write(* ,'(2F12.6)') ena(i)   ,q(i)
c        write(43,'(2F12.6)') ena(i)*10,q(i)*10
c     enddo
c     close(43)
c     stop
c     endif

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c intial guess for P in GFN2 runs
c assuming canonical occupation
c corrected for partial EN-based charges
c C. Bannwarth : modified from S. Grimme's setzshell routine
c n: number of atoms
c nao: # of MOs
c at: ordinal numbers of atoms
c z: nuclear charges
c q: partial charges
c P: guess for density matrix (output)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! for inital guess for P, add this after P allocation in main.f
!      if(gfnver.gt.1) call inip(n,nao,at,z,q,P)

      subroutine inip(n,nao,at,z,q,P)
      implicit none
      include 'aoelementcommon.fh'
      include 'ehtcommon.fh'

      integer n,nao
      integer at(n)
      real*8  z(n),P0(nao)
      real*8  q(n),P(nao,nao),qtot,fred

      integer lin,kk,k,idum,i,j,iat,ntot,m,l,lll,ll(0:3)
      data ll /1,3,5,7/
      call inip0(n,nao,at,P0)
      ! canonical orbital occupation
      P=0.0d0
      k=0
      j=0
      do i=1,n
         iat=at(i)
         qtot=-1.d-6
         ntot=0
         do m=1,ao_n(iat)
            l=ao_l(m,iat)
            do lll=1,ll(l)
               k=k+1
               j=j+1
               ntot=ntot+1
               P(k,k)=P0(k) ! ref densmat is diagonal
            enddo
         enddo
         j=j-ntot
         fred=1.0d0/dble(ntot)
         ! correct for partial charge
         do m=1,ao_n(iat)
            l=ao_l(m,iat)
            do lll=1,ll(l)
               j=j+1
               P(j,j)=P(j,j)-q(i)*fred
            enddo
         enddo
      enddo
      if(k.ne.nao) stop 'internal error'
!      call prmat(6,p,nao,nao,'Pguess')
      end subroutine inip


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c intial guess for P in GFN2 runs
c assuming canonical occupation
c corrected for partial EN-based charges
c C. Bannwarth : modified from S. Grimme's setzshell routine
c n: number of atoms
c nao: # of MOs
c at: ordinal numbers of atoms
c z: nuclear charges
c q: partial charges
c P: guess for density matrix (output)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! for inital guess for P, add this after P allocation in main.f
!      if(gfnver.gt.1) call inip(n,nao,at,z,q,P)

      subroutine inip0(n,nao,at,P0)
      implicit none
      include 'aoelementcommon.fh'
      include 'ehtcommon.fh'

      integer n,nao
      integer at(n)
      real*8  P0(nao),qtot,fred

      integer lin,kk,k,idum,i,j,iat,ntot
      integer iox(86,0:3),m,l,lll,ll(0:3)
      data ll /1,3,5,7/
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

     & 0,                                                         0, !He
     & 0,0,                                        0, 0, 0, 0, 0, 0, !Ne
     & 0,0,                                        0, 0, 0, 0, 0, 0, !Ar
     & 0,0,0,          0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, !Kr
     & 0,0,0,          0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, !Xe
     & 0,0,0,
c    &       1,2,3,4,5,6,7,8,9,10,11,12,13,14,                       !f block
     &       14*0,                                                   !f block
     &                 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0/ !Rn

      ! canonical orbital occupation
      P0=0.0d0
      k=0
      j=0
      do i=1,n
         iat=at(i)
         qtot=-1.d-6
         ntot=0
         do m=1,ao_n(iat)
            l=ao_l(m,iat)
            do lll=1,ll(l)
               k=k+1
               P0(k)=dble(iox(iat,l))/dble(ll(l))
               qtot=qtot+P0(k)
            enddo
         enddo
      enddo
      if(k.ne.nao) stop 'internal error'
!      call prmat(6,p,nao,nao,'Pguess')
      end subroutine inip0

