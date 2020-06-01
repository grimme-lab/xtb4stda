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

      module gbobc
      implicit none
c     cutoffs
c     Born radii
      real*8 :: lrcut_a=35.d0
c     SASA = 2*(w+maxrasasa) + srcut_add
      real*8 :: srcut_add=2.d0
c     flag for gbsa
      logical :: lgbsa
c     number of particles
      integer :: nat
c     number of pairs
      integer :: ntpair

c     van der Waals radii of the particles
      real*8, allocatable :: vdwr(:)
c     greatest van der Waals radius
      real*8 :: maxvdwr

c     pair descreening approximation radii
      real*8, allocatable :: rho(:)
c     offset van der Waals radii
      real*8, allocatable :: svdw(:)

c     Neighbor list:
c     cut-off radius for the Born radius NN list
      real*8 :: lrcut
c     cut-off radius for the SASA NN list
      real*8 :: srcut
c     number of neighbors for Born radii
      integer :: nnrad
c     number of neighbors for SASA computation
      integer, allocatable :: nnsas(:)
c     neighbors of an atom for Born radii
      integer, allocatable :: nnlistr(:,:)
c     neighbors of an atom for SASA
      integer, allocatable :: nnlists(:,:)
c     all pairs indeces array
      integer, allocatable :: ppind(:,:)
c     all pairs vector differences and magnitudes array
      real*8, allocatable ::  ddpair(:,:)

c     angstroem to atomic unit
      real*8, parameter :: atoau=1.d0/0.52917726d0
c     atomic units to kcal
      real*8, parameter :: autokcal=627.509541d0
c     4*pi
      real*8, parameter :: pi4=12.56637061435916D0
c     a.u. to eV
      real*8, parameter :: autoeV=27.21138505d0

c     GBOBC parameters
c     offset parameter (fitted)
c     real*8 :: soset=0.09d0*atoau
      real*8 :: soset
c     van der Waals to Lee-Richard's surface correction
      real*8 :: alp=1.d0
      real*8 :: bet=0.8d0
      real*8 :: gam=4.85d0

c     Smoothing dielectric function parameters
      real*8, parameter :: w=0.3d0*atoau
      real*8, parameter :: w3=w*(w*w)
      real*8, parameter :: ah0=0.5d0
      real*8, parameter :: ah1=3.d0/(4.d0*w)
      real*8, parameter :: ah3=-1.d0/(4.d0*w3)

c     Angular grid. 38 points lead rot.inv.errors for c60 of 2 kcal
      integer, parameter :: nangsa=230
      real*8 :: grida(4,nangsa)
      include 'grida230.fh'
c     integer, parameter :: nangsa=86
c     include 'grida86.fh'
c     integer, parameter :: nangsa=110
c     include 'grida110.fh'

c     real space cut-offs
      real*8, parameter :: tolsesp=1.d-6

c     Atom specific surface data
      real*8, allocatable :: vdwsa(:)
      real*8, allocatable :: wrp(:)
      real*8, allocatable :: trj2(:,:)

c     Dielectric data
      real*8 :: gborn
      real*8 :: epsv
      real*8 :: epsu
      real*8 :: keps

c     Surface tension (mN/m=dyn/cm)
      real*8, parameter :: mNmkcal=4.0305201015221386d-4
      real*8 :: gammas
      real*8 :: gamscale(94)

c     Solvent density (g/cm^3) and molar mass (g/mol)
      real*8, parameter :: molcm3au=8.92388d-2
      real*8 :: smass
      real*8 :: rhos

c     Born radii
      real*8 :: c1
      real*8, allocatable :: brad(:)
      real*8, allocatable :: brt(:)

c     Salt screening
      logical :: lsalt=.false.
      real*8  :: ionst=0.d0
      real*8  :: kappa_const=0.7897d-3
      real*8  :: ion_rad=0.d0
      real*8  :: kappa=0.d0
      real*8, allocatable :: ionscr(:)
      real*8, allocatable :: discr(:)

c     Atomic surfaces
      real*8 :: rprobe
      real*8 :: sasamol
      real*8 :: gsasa
      real*8 :: sasagam
      real*8, allocatable :: gamsasa(:)
      real*8, allocatable :: sasa(:)

c     Hydrogen bond contribution
      logical :: lhb=.true.
      real*8 :: ghb
      real*8, allocatable :: hbw(:)

c     Gradient:
c     Born radii gradient
      real*8, allocatable :: brdr(:,:,:)

c     Molecular Surface gradient
      real*8, allocatable :: dsdr(:,:)
      real*8, allocatable :: dsdrp(:,:)
      real*8, allocatable :: dsdrt(:,:,:)

c     Hydrogen bond gradient
      real*8, allocatable :: dhbdw(:)

c     GB energy gradient
      real*8, allocatable :: grdgb(:,:)
      real*8, allocatable :: grddb(:)
      real*8, allocatable :: dbrdp(:)
      real*8, allocatable :: dgbta(:,:)
      real*8, allocatable :: dgbtb(:)

c     Parameters:
c     van der Waals radii
      real*8 :: rvdw(94)
c     dielectric descreening parameters
      real*8 :: sx(94)
c     solvent accesible surface radii
      real*8 :: rasasa(94)
c     HB correction present if zero no HB correction
      integer :: at_hb(94)
c     solvent HB donor or acceptor strength
      real*8 :: hb_mag(94)
c     Gshift (gsolv=reference vs. gsolv)
      real*8 :: gshift

      integer, allocatable :: at(:)
      save

      contains

      subroutine rd_rvdw(nat,at,sname,mode,temp)
      implicit none
      include 'setcommon.fh'
      character*(*) sname
      integer nat,at(nat),mode

      integer :: i,fix,inum
      real*8 :: rad
      real*8 :: gamma_in, rvdwscal, tmp(94), gstate, dum, temp
      character*80 fname
      logical ex
      character*200 a200

c     D3 cut-off radii
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

c     hydrogen bonding parameters
      lhb=.false.

      at_hb=0
      at_hb(1)=1
      at_hb(6)=1
      at_hb(7)=1
      at_hb(8)=1
      at_hb(9)=1
      at_hb(15)=1
      at_hb(16)=1
      at_hb(17)=1
      at_hb(34)=1
      at_hb(35)=1
      at_hb(53)=1

      rvdwscal=1.0d0

      write(fname,'(''.param_gbsa_'',a)')trim(sname)
      a200=trim(XTB4STDAHOME) // trim(fname)
      fname=a200
      write(*,*) 'Solvent             : ', trim(sname)
      write(*,*) 'GBSA parameter file : ', trim(fname)

      inquire(file=fname,exist=ex)
      if(.not.ex)then
        write(*,*) 'solvent :',trim(sname),' not implemented'
        stop 'init_gbsa'
      endif

      open(unit=1,file=fname)
      read(1,*)epsv
      read(1,*)smass
      read(1,*)rhos
      read(1,*)c1
      read(1,*)rprobe
      read(1,*)gshift
      read(1,*)soset
      read(1,*)dum

      if(mode.eq.1) then ! gsolv=reference option in COSMOTHERM
c               RT*(ln(ideal gas mol volume)+ln(rho/M))
         gstate=(temp*8.31451/1000./4.184)*
     .      (log(24.79d0*temp/298.15)+
     .       log(1000.0d0*rhos/smass))
         gshift=(gshift+gstate)/autokcal
         write(*,*) 'Gsolv state corr. (kcal):',gstate
         a200='gsolv=reference [X=1]'
      elseif(mode.eq.0)then !gsolv option in COSMOTHERM to which it was fitted
         gshift=gshift/autokcal
         a200='gsolv [1 M gas/solution]'
      elseif(mode.eq.2)then ! 1 bar gas/ 1 M solution is not implemented in COSMOTHERM although its the canonical choice
         gstate=(temp*8.31451/1000./4.184)*log(24.79d0*temp/298.15)
         gshift=(gshift+gstate)/autokcal
         write(*,*) 'Gsolv state corr. (kcal):',gstate
         a200='gsolv [1 bar gas/ 1 M solution]'
      endif

      do i=1,94
         read(1,*)gamscale(i),sx(i),tmp(i)
         if(abs(tmp(i)).gt.1.d-3) lhb=.true.
      enddo
      close(1)

c     if(fit)then !penalty to avoid small sx which lead to numerical instabs
c     dum=0
c     do i=1,nat
c        dum=dum+2.*(sx(at(i))-0.8)**4
c     enddo
c     gshift=gshift+dum/autokcal
c     endif

c     hydrogen bonding magnitude
      hb_mag = -(tmp**2)/autokcal

c     scaling of the van der Waals radius
      rvdw = rvdw * rvdwscal

c     add the probe radius to the molecular surface
      rasasa=rvdw+rprobe

c     surface tension scaling
      gamma_in=1.0d0
      gammas=gamma_in*(1.0d-5)*autokcal/mNmkcal

c     dielectric scaling
      epsu=1.d0
      keps=((1.d0/epsv)-(1.d0/epsu))

c     set the salt term
      if(lsalt) then
c      convert to au
       ion_rad=ion_rad*atoau
c      inverse Debye screening length
       kappa=sqrt(epsv*temp*kappa_const/ionst)*atoau
       kappa=1.d0/kappa
      endif

c     print parameters
      write(*,*) 'Gsolv ref. state (COSMO-RS): ',trim(a200)
      write(*,*) 'temperature (mdtemp)       : ',temp
      write(*,*) 'dielectric constant        : ',epsv
      write(*,*) 'rho                        : ',rhos
      write(*,*) 'mass                       : ',smass
      write(*,*) 'surface tension            : ',gammas
      write(*,*) 'probe radius               : ',rprobe
      write(*,*) 'vdW radii scaling          : ',rvdwscal
      write(*,*) 'Gshift (Eh)                : ',gshift
      write(*,*) 'c1                         : ',c1
      write(*,*) 'soset                      : ',soset
      write(*,*) 'HB correction              : ',lhb
      if(lsalt) then
      write(*,*) 'Debye screening length     : ',1.d0/kappa/atoau
      endif

      soset=0.1*soset/0.52917726d0

      rhos=rhos*molcm3au/smass

      return
      end subroutine rd_rvdw

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine init_gbsa(n12,at12,sname,mode,temp)
       implicit none
       character(20) :: sname

       integer n12,mode
       integer at12(n12)
       real*8 temp

       integer i,j,k
       integer ierr
       real*8 minvdwr
       real*8 maxrasasa
       real*8 r

        call rd_rvdw(n12,at12,sname,mode,temp)

        nat=n12

c       initialize the vdw radii array
        allocate(vdwr(nat))
        allocate(rho(nat))
        allocate(svdw(nat))
        allocate(at(nat))
        at=at12
        maxvdwr=0.d0
        minvdwr=1000.d0
        do i=1,nat
         vdwr(i)=rvdw(at12(i))*atoau
         rho(i)=vdwr(i)*sx(at12(i))
         svdw(i)=vdwr(i)-soset
         maxvdwr=max(maxvdwr,vdwr(i))
         minvdwr=min(minvdwr,vdwr(i))
        enddo

c       initialize Born radii
        allocate(brad(nat),brt(nat))
        allocate(brdr(3,nat,nat))
        allocate(grdgb(3,nat),grddb(nat))
        allocate(dbrdp(nat))
        allocate(dgbta(3,nat),dgbtb(nat))

c       nearest-neighbor list preparation
        lrcut = lrcut_a*atoau

        ntpair=nat*(nat-1)/2

        allocate(nnlistr(3,ntpair),nnsas(nat),nnlists(nat,nat),
     .           ppind(2,ntpair),stat=ierr)
        if(ierr.ne.0) stop 'Memory allocation failed!'
        allocate(ddpair(4,ntpair),stat=ierr)
        if(ierr.ne.0) stop 'Memory allocation failed!'

        k=0
        do i=1,nat
         do j = 1,i-1
           k=k+1
           ppind(1,k)=i
           ppind(2,k)=j
         enddo
        enddo

c       initialize solvent-accessible atomic surface area computation (SASA)
        allocate(vdwsa(nat))
        allocate(wrp(nat))
        allocate(trj2(2,nat))

        allocate(sasa(nat),gamsasa(nat))
        allocate(dsdr(3,nat),dsdrp(3,nat),dsdrt(3,nat,nat))

        maxrasasa=0.d0
        do i = 1, nat
         vdwsa(i) = rasasa(at12(i))*atoau
         maxrasasa=max(maxrasasa,vdwsa(i))
         trj2(1,i) = (vdwsa(i)-w)**2
         trj2(2,i) = (vdwsa(i)+w)**2
         r=vdwsa(i)+w
         wrp(i)=(0.25d0/w+
     .            3.d0*ah3*(0.2d0*r*r-0.5*r*vdwsa(i)+
     .            vdwsa(i)*vdwsa(i)/3.))*r*r*r
         r=vdwsa(i)-w
         wrp(i)=wrp(i)-(0.25/w+
     .    3.d0*ah3*(0.2d0*r*r-0.5*r*vdwsa(i)+
     .            vdwsa(i)*vdwsa(i)/3.))*r*r*r
        enddo

        srcut = 2.d0*(w + maxrasasa) + srcut_add*atoau

        gammas=gammas*mNmkcal/autokcal
        sasagam=pi4*gammas
        do i = 1, nat
         gamsasa(i)=gamscale(at12(i))*pi4*gammas
        enddo

c       initialize the hydrogen bonding contribution
        ghb=0.d0
        if(lhb) then
         allocate(hbw(nat),dhbdw(nat))
        endif

c       initialize the salt term
        if(lsalt) then
         allocate(ionscr(nat),discr(nat))
        endif

       return
       end subroutine init_gbsa

       subroutine compute_fgb(n,xyz,fgb,fhb)
       implicit none

       integer n
       real*8 xyz(3,n)
       real*8 fgb(n,n)
       real*8 fhb(n)

       integer i,j,nnj
       integer kk
       real*8, parameter :: a13=1.d0/3.d0
       real*8, parameter :: a4=0.25d0
       real*8 aa,r2,gg,iepsu
       real*8 dd,edd,dfgb,hkeps

c      initialize
       fgb=0.d0

       hkeps=keps*autoeV

c      compute Born radii
       call compute_brad_sasa(n,xyz)

c      compute the Debye-Hueckel ion exclusion term
       if(lsalt) then
        aa=0.5d0/epsv
        do i = 1, n
         gg=kappa*(brad(i)+ion_rad)
         ionscr(i)=aa*exp(gg)/(1.d0+gg)
         discr(i)=ionscr(i)*kappa*gg/(1.d0+gg)
        enddo
       endif

       if(lsalt) then

       iepsu=1.d0/epsu

c      compute energy and fgb direct and radii derivatives
!$OMP PARALLEL PRIVATE(i,j,r2,aa,dd,edd,dfgb)
!$OMP DO
       do kk = 1, ntpair
         r2=ddpair(1,kk)
         r2=r2*r2

         i=ppind(1,kk)
         j=ppind(2,kk)

         aa=brad(i)*brad(j)
         dd=a4*r2/aa
         edd=exp(-dd)
         dfgb=sqrt(r2+aa*edd)
         gg=ionscr(i)+ionscr(j)
         fgb(i,j)=autoeV*(exp(-kappa*dfgb)*gg-iepsu)/dfgb
         fgb(j,i)=fgb(i,j)
       enddo
!$OMP ENDDO
!$OMP END PARALLEL

c      self-energy part
       do i = 1, n
        gg=ionscr(i)*2.d0
        fgb(i,i)=autoeV*(exp(-kappa*brad(i))*gg-iepsu)/brad(i)
       enddo

       else

c      compute energy and fgb direct and radii derivatives
!$OMP PARALLEL PRIVATE(i,j,r2,aa,dd,edd,dfgb)
!$OMP DO
       do kk = 1, ntpair
         r2=ddpair(1,kk)
         r2=r2*r2

         i=ppind(1,kk)
         j=ppind(2,kk)

         aa=brad(i)*brad(j)
         dd=a4*r2/aa
         edd=exp(-dd)
         dfgb=1.d0/(r2+aa*edd)
         fgb(i,j)=hkeps*sqrt(dfgb)
         fgb(j,i)=fgb(i,j)
       enddo
!$OMP ENDDO
!$OMP END PARALLEL

c      self-energy part
       do i = 1, n
        fgb(i,i)=hkeps/brad(i)
       enddo

       endif

c      compute the HB term
       if(lhb) then
         call compute_fhb(n,xyz,fhb)
       else
         fhb=0.d0
       endif

       return
       end subroutine compute_fgb

       subroutine compute_gb_egrad(n,xyz,q,gbornh,grd,lpr)
       implicit none

       integer n
       real*8 xyz(3,n)
       real*8 q(n)
       real*8 gbornh
       real*8 grd(3,n)
       logical lpr

       integer i,j,nnj
       integer kk
       real*8, parameter :: a13=1.d0/3.d0
       real*8, parameter :: a4=0.25d0
       real*8 aa,r2,fgb,br3
       real*8 qq,dd,edd,dfgb,egb,ap,bp,qfg
       real*8 gg,efg,epu
       real*8 r0vdw,r01,r02,ar02
       real*8 r(3)
       real*8 kq(n)

c      GB energy and gradient
       if(.not.lsalt) then

c      dielectric scaling of the charges
       kq=keps*q

       gborn=0.d0
       grdgb=0.d0
       grddb=0.d0

c      compute energy and fgb direct and radii derivatives
!$OMP PARALLEL PRIVATE(i,j,r2,r,qq,aa,dd,edd,fgb,dfgb,egb,ap,bp),
!$OMP& PRIVATE(dgbta,dgbtb),REDUCTION(+:gborn,grdgb,grddb)
       egb=0.d0
       dgbta=0.d0
       dgbtb=0.d0
!$OMP DO
       do kk = 1, ntpair
         r2=ddpair(1,kk)
         r2=r2*r2

         i=ppind(1,kk)
         j=ppind(2,kk)

         qq=q(i)*kq(j)
         aa=brad(i)*brad(j)
         dd=a4*r2/aa
         edd=exp(-dd)
         dfgb=1.d0/(r2+aa*edd)
         fgb=qq*sqrt(dfgb)
         dfgb=dfgb*fgb

         egb=egb+fgb

         ap=(1.d0-a4*edd)*dfgb
         r(1)=ap*ddpair(2,kk)
         r(2)=ap*ddpair(3,kk)
         r(3)=ap*ddpair(4,kk)
         dgbta(1,i)=dgbta(1,i)-r(1)
         dgbta(2,i)=dgbta(2,i)-r(2)
         dgbta(3,i)=dgbta(3,i)-r(3)
         dgbta(1,j)=dgbta(1,j)+r(1)
         dgbta(2,j)=dgbta(2,j)+r(2)
         dgbta(3,j)=dgbta(3,j)+r(3)

         bp=-0.5d0*edd*(1.d0+dd)*dfgb
         dgbtb(i)=dgbtb(i)+brad(j)*bp
         dgbtb(j)=dgbtb(j)+brad(i)*bp
       enddo
!$OMP ENDDO
       gborn=gborn+egb
       grdgb=grdgb+dgbta
       grddb=grddb+dgbtb
!$OMP END PARALLEL

c      self-energy part
       do i = 1, n
        qq=q(i)/brad(i)
        gborn=gborn+0.5d0*kq(i)*qq
        grddb(i)=dbrdp(i)*(grddb(i)-0.5d0*keps*qq*qq)
       enddo

c      contract with the Born radii derivatives
!$OMP PARALLEL PRIVATE(dgbta,i,j), REDUCTION(+:grdgb)
       dgbta=0.d0
!$OMP DO
       do kk = 1, nnrad
         i=nnlistr(1,kk)
         j=nnlistr(2,kk)
         dgbta(1,i)=dgbta(1,i)+grddb(j)*brdr(1,j,i)
         dgbta(2,i)=dgbta(2,i)+grddb(j)*brdr(2,j,i)
         dgbta(3,i)=dgbta(3,i)+grddb(j)*brdr(3,j,i)
         dgbta(1,j)=dgbta(1,j)+grddb(i)*brdr(1,i,j)
         dgbta(2,j)=dgbta(2,j)+grddb(i)*brdr(2,i,j)
         dgbta(3,j)=dgbta(3,j)+grddb(i)*brdr(3,i,j)
       enddo
!$OMP ENDDO
       grdgb=grdgb+dgbta
!$OMP END PARALLEL

c      self-energy part
       do i = 1, n
        grdgb(1,i)=grdgb(1,i)+grddb(i)*brdr(1,i,i)
        grdgb(2,i)=grdgb(2,i)+grddb(i)*brdr(2,i,i)
        grdgb(3,i)=grdgb(3,i)+grddb(i)*brdr(3,i,i)
       enddo

       else
c      GB-SE energy and gradient
       gborn=0.d0
       grdgb=0.d0
       grddb=0.d0

       epu=1.d0/epsu

c      compute energy and fgb direct and radii derivatives
!$OMP PARALLEL PRIVATE(i,j,r2,r,qq,aa,dd,edd,fgb,dfgb,egb,ap,bp,qfg),
!$OMP& PRIVATE(gg,efg,dgbta,dgbtb),REDUCTION(+:gborn,grdgb,grddb)
       egb=0.d0
       dgbta=0.d0
       dgbtb=0.d0
!$OMP DO
       do kk = 1, ntpair
         r2=ddpair(1,kk)
         r2=r2*r2

         i=ppind(1,kk)
         j=ppind(2,kk)

         qq=q(i)*q(j)
         aa=brad(i)*brad(j)
         dd=a4*r2/aa
         edd=exp(-dd)
         dfgb=r2+aa*edd
         fgb=sqrt(dfgb)
         aa=kappa*fgb
         efg=exp(-aa)
         gg=(ionscr(i)+ionscr(j))*efg
         qfg=qq/fgb

         egb=egb+qfg*(gg-epu)

         dfgb=qfg*(gg*(1.d0+aa)-epu)/dfgb

         ap=(1.d0-a4*edd)*dfgb
         r(1)=ap*ddpair(2,kk)
         r(2)=ap*ddpair(3,kk)
         r(3)=ap*ddpair(4,kk)
         dgbta(1,i)=dgbta(1,i)-r(1)
         dgbta(2,i)=dgbta(2,i)-r(2)
         dgbta(3,i)=dgbta(3,i)-r(3)
         dgbta(1,j)=dgbta(1,j)+r(1)
         dgbta(2,j)=dgbta(2,j)+r(2)
         dgbta(3,j)=dgbta(3,j)+r(3)

         qfg=qfg*efg
         bp=-0.5d0*edd*(1.d0+dd)*dfgb
         dgbtb(i)=dgbtb(i)+brad(j)*bp+qfg*discr(i)
         dgbtb(j)=dgbtb(j)+brad(i)*bp+qfg*discr(j)
       enddo
!$OMP ENDDO
       gborn=gborn+egb
       grdgb=grdgb+dgbta
       grddb=grddb+dgbtb
!$OMP END PARALLEL

c      self-energy part
       do i = 1, n
        gg=exp(-kappa*brad(i))
        aa=2.d0*ionscr(i)*gg-epu
        qq=q(i)/brad(i)
        gborn=gborn+0.5d0*qq*q(i)*aa
        ap=aa-brad(i)*2.d0*(discr(i)+ionscr(i)*kappa)*gg
        grddb(i)=dbrdp(i)*(grddb(i)-0.5d0*qq*qq*ap)
       enddo

c      contract with the Born radii derivatives
!$OMP PARALLEL PRIVATE(dgbta,i,j), REDUCTION(+:grdgb)
       dgbta=0.d0
!$OMP DO
       do kk = 1, nnrad
         i=nnlistr(1,kk)
         j=nnlistr(2,kk)
         dgbta(1,i)=dgbta(1,i)+grddb(j)*brdr(1,j,i)
         dgbta(2,i)=dgbta(2,i)+grddb(j)*brdr(2,j,i)
         dgbta(3,i)=dgbta(3,i)+grddb(j)*brdr(3,j,i)
         dgbta(1,j)=dgbta(1,j)+grddb(i)*brdr(1,i,j)
         dgbta(2,j)=dgbta(2,j)+grddb(i)*brdr(2,i,j)
         dgbta(3,j)=dgbta(3,j)+grddb(i)*brdr(3,i,j)
       enddo
!$OMP ENDDO
       grdgb=grdgb+dgbta
!$OMP END PARALLEL

c      self-energy part
       do i = 1, n
        grdgb(1,i)=grdgb(1,i)+grddb(i)*brdr(1,i,i)
        grdgb(2,i)=grdgb(2,i)+grddb(i)*brdr(2,i,i)
        grdgb(3,i)=grdgb(3,i)+grddb(i)*brdr(3,i,i)
       enddo
       endif

       gbornh = gborn
       grd = grd + grdgb + dsdr

       if(lhb) then
        call compute_hb(n,q,grd)
       endif

c      if(lopt.and.lpr) then
c       write(*,'(/,a)') 'Results GBOBC:'
c       write(*,*) 'At #,  Z , GBOBC (A), RVDW (A)'
c       do i = 1, nat
c        write(*,'(I5,2x,I2,6F12.4)') i,at(i),brad(i)/atoau,
c    .   rvdw(at(i)),sx(at(i)),xyz(1:3,i)/atoau
c       enddo
c       write(*,'(/,a)') 'Free Energy (kcal/mol):'
c       write(*,'(''G-EL  = '',F8.3)') gborn*autokcal
c       write(*,'(''GCAV  = '',F8.3)') gsasa*autokcal
c       write(*,'(''G-HB  = '',F8.3)') ghb*autokcal
c       write(*,'(''GSOL  = '',F8.3)') (gborn+gsasa+ghb)*autokcal
c      endif

       return
       end subroutine compute_gb_egrad

       subroutine compute_fhb(n,xyz,fhb)
       implicit none

       integer n
       real*8 :: xyz(3,n)
       real*8 :: fhb(n)

       integer :: i

        hbw=0.d0
        dhbdw=0.d0

!$OMP PARALLEL DO
        do i = 1, n
          call compute_fhb_i(i,n,xyz)
        enddo
!$OMP END PARALLELDO

        fhb=hbw*autoeV

       return
       end subroutine compute_fhb

       subroutine compute_fhb_i(i,n,xyz)
       implicit none

       integer :: i,n
       real*8 :: xyz(3,n)

       integer :: iz,nhb
       real*8 :: hbed,dhbed
       real*8 :: smaxd,sasad,sasaw
       real*8 :: sfw,dsfw,w3,w2,w1
       integer :: j
       real*8 :: wbh,wah

c       atomic Z
        iz=at(i)
c       number of HB
        nhb=at_hb(iz)
        if(nhb.gt.0) then
c        SASA-D for HB
         smaxd=1.d0/(vdwsa(i)*vdwsa(i)*gamsasa(i))
         sasad=sasa(i)*smaxd
         hbw(i)=hb_mag(iz)*sasad
         dhbdw(i)=hb_mag(iz)*smaxd
        endif

       return
       end subroutine compute_fhb_i

       subroutine compute_hb(n,q,grd)
       implicit none

       integer :: n
       real*8 :: q(n)
       real*8 :: ghbh
       real*8 :: grd(3,n)

       integer :: i,j
       real*8 :: dhbed
       real*8 :: qq

       ghb=0.d0
       do i = 1, n
        qq = q(i)*q(i)
        ghb = ghb + hbw(i)*qq
       enddo

!$OMP PARALLEL PRIVATE(dsdrp,j,dhbed), REDUCTION(+:grd)
       dsdrp=0.d0
!$OMP DO
       do i = 1, n
         dhbed=dhbdw(i)
         if(abs(dhbed).gt.0.d0) then
          dhbed=dhbed*(q(i)*q(i))
          do j = 1, n
           dsdrp(1,j) = dsdrp(1,j) + dsdrt(1,j,i)*dhbed
           dsdrp(2,j) = dsdrp(2,j) + dsdrt(2,j,i)*dhbed
           dsdrp(3,j) = dsdrp(3,j) + dsdrt(3,j,i)*dhbed
          enddo
         endif
       enddo
!$OMP ENDDO
       grd=grd+dsdrp
!$OMP END PARALLEL

       return
       end subroutine compute_hb

       subroutine compute_brad_sasa(n,xyz)
       implicit none

       integer n
       real*8 :: xyz(3,n)

       integer i,j,kk
       real*8 brdrd(3,n)
       real*8 brdrt(3,n)

       brad=0.d0
       dsdr=0.d0
       dsdrt=0.d0
       brdr=0.d0
       brdrd=0.d0

c      compute Born radii and their derivatives
!$OMP PARALLEL PRIVATE(brt,brdrt), REDUCTION(+:brad,brdrd)
       brt=0.d0
       brdrt=0.d0
!$OMP DO SCHEDULE(DYNAMIC,1)
       do kk = 1, nnrad
         call compute_psi(kk,brt,brdr,brdrt)
       enddo
!$OMP ENDDO
       brad=brad+brt
       brdrd=brdrd+brdrt
!$OMP END PARALLEL

       do i = 1, nat
        brdr(1:3,i,i)=brdrd(1:3,i)
       enddo

!$OMP PARALLELDO
       do i = 1, nat
        call compute_bornr(i,brad(i),dbrdp(i))
       enddo
!$OMP END PARALLELDO

c      compute solvent accessible surface and its derivatives
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
       do i = 1, nat
        call compute_numsa(n,xyz,i,vdwsa(i),sasa(i),dsdrt(:,:,i))
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL PRIVATE(dsdrp,j), REDUCTION(+:dsdr)
      dsdrp=0.d0
!$OMP DO
       do i = 1, nat
        do j = 1, nat
         dsdrp(1,j) = dsdrp(1,j) + dsdrt(1,j,i)
         dsdrp(2,j) = dsdrp(2,j) + dsdrt(2,j,i)
         dsdrp(3,j) = dsdrp(3,j) + dsdrt(3,j,i)
        enddo
       enddo
!$OMP ENDDO
       dsdr=dsdr+dsdrp
!$OMP END PARALLEL

       gsasa = sum(sasa)

       return
       end subroutine compute_brad_sasa

       subroutine compute_bornr(iat,br,dpsi)
       implicit none

       integer iat
       real*8 br
       real*8 dpsi
       real*8 svdwi,vdwri
       real*8 s1,v1,s2
       real*8 arg,arg2,th,ch
       real*8 alpi,beti,gami

        svdwi=svdw(iat)
        vdwri=vdwr(iat)
        s1=1.d0/svdwi
        v1=1.d0/vdwri
        s2=0.5d0*svdwi

        br=br*s2

        arg2=br*(gam*br-bet)
        arg=br*(alp+arg2)
        arg2=2.d0*arg2+alp+gam*br*br

        th=tanh(arg)
        ch=cosh(arg)

        br=1.d0/(s1-v1*th)
c       Include GBMV2-like scaling
        br=c1*br

        dpsi=ch*(s1-v1*th)
        dpsi=s2*v1*arg2/(dpsi*dpsi)
        dpsi=c1*dpsi

       return
       end subroutine compute_bornr

       subroutine compute_psi(kk,br,gbr,gbrd)
       implicit none

       integer kk
       real*8 br(nat),gbr(3,nat,nat),gbrd(3,nat)
       integer ii,jj,nn
       real*8 dr(3),r,rhoi,rhoj
       real*8 gi,gj,ap,am,lnab,rhab,ab,dgi,dgj
       real*8 drjj(3)
       real*8 rh1,rhr1,r24,rh2,r1,aprh1,r12
       real*8 rvdwi,rvdwj
       integer ovij,ovji,ov

        ii=nnlistr(1,kk)
        jj=nnlistr(2,kk)
        nn=nnlistr(3,kk)

        r=ddpair(1,nn)
        dr(1)=ddpair(2,nn)
        dr(2)=ddpair(3,nn)
        dr(3)=ddpair(4,nn)

        rhoi=rho(ii)
        rhoj=rho(jj)
        rvdwi=vdwr(ii)
        rvdwj=vdwr(jj)

        ovij=1
        ovji=1
        if(r.ge.(rvdwi+rhoj)) ovij=0
        if(r.ge.(rhoi+rvdwj)) ovji=0
        ov=ovij+10*ovji

        select case(ov)
c        ij do not overlap; ji do not overlap
         case(0)
c        nonoverlaping spheres
         if(abs(rhoi-rhoj).lt.1.d-8) then
c         equal reduced radii
          r1=1.d0/r
          ap=r+rhoj
          am=r-rhoj
          ab=ap*am
          rhab=rhoj/ab
          lnab=0.5d0*log(am/ap)*r1
          gi=rhab+lnab
          dgi=-2.d0*rhab/ab+(rhab-lnab)*r1*r1
c         accumulate psi
          br(ii)=br(ii)+gi
          br(jj)=br(jj)+gi
c         accumulate psi gradient
          drjj(1)=dgi*dr(1)
          drjj(2)=dgi*dr(2)
          drjj(3)=dgi*dr(3)
          gbrd(1,ii)=gbrd(1,ii)+drjj(1)
          gbrd(2,ii)=gbrd(2,ii)+drjj(2)
          gbrd(3,ii)=gbrd(3,ii)+drjj(3)
          gbr(1,ii,jj)=gbr(1,ii,jj)-drjj(1)
          gbr(2,ii,jj)=gbr(2,ii,jj)-drjj(2)
          gbr(3,ii,jj)=gbr(3,ii,jj)-drjj(3)
          gbrd(1,jj)=gbrd(1,jj)-drjj(1)
          gbrd(2,jj)=gbrd(2,jj)-drjj(2)
          gbrd(3,jj)=gbrd(3,jj)-drjj(3)
          gbr(1,jj,ii)=gbr(1,jj,ii)+drjj(1)
          gbr(2,jj,ii)=gbr(2,jj,ii)+drjj(2)
          gbr(3,jj,ii)=gbr(3,jj,ii)+drjj(3)
         else
c         unequal reduced radii
c         ij contribution
          r1=1.d0/r
          ap=r+rhoj
          am=r-rhoj
          ab=ap*am
          rhab=rhoj/ab
          lnab=0.5d0*log(am/ap)*r1
          gi=rhab+lnab
          dgi=-2.d0*rhab/ab+(rhab-lnab)*r1*r1
c         ji contribution
          ap=r+rhoi
          am=r-rhoi
          ab=ap*am
          rhab=rhoi/ab
          lnab=0.5d0*log(am/ap)*r1
          gj=rhab+lnab
          dgj=-2.d0*rhab/ab+(rhab-lnab)*r1*r1
c         accumulate psi
          br(ii)=br(ii)+gi
          br(jj)=br(jj)+gj
c         accumulate psi gradient
          drjj(1)=dgi*dr(1)
          drjj(2)=dgi*dr(2)
          drjj(3)=dgi*dr(3)
          gbrd(1,ii)=gbrd(1,ii)+drjj(1)
          gbrd(2,ii)=gbrd(2,ii)+drjj(2)
          gbrd(3,ii)=gbrd(3,ii)+drjj(3)
          gbr(1,ii,jj)=gbr(1,ii,jj)-drjj(1)
          gbr(2,ii,jj)=gbr(2,ii,jj)-drjj(2)
          gbr(3,ii,jj)=gbr(3,ii,jj)-drjj(3)

          drjj(1)=dgj*dr(1)
          drjj(2)=dgj*dr(2)
          drjj(3)=dgj*dr(3)
          gbrd(1,jj)=gbrd(1,jj)-drjj(1)
          gbrd(2,jj)=gbrd(2,jj)-drjj(2)
          gbrd(3,jj)=gbrd(3,jj)-drjj(3)
          gbr(1,jj,ii)=gbr(1,jj,ii)+drjj(1)
          gbr(2,jj,ii)=gbr(2,jj,ii)+drjj(2)
          gbr(3,jj,ii)=gbr(3,jj,ii)+drjj(3)
         endif
c        ij do not overlap; ji overlap
         case(10)

c         ij contribution
          r1=1.d0/r
          ap=r+rhoj
          am=r-rhoj
          ab=ap*am
          rhab=rhoj/ab
          lnab=0.5d0*log(am/ap)*r1
          gi=rhab+lnab
          dgi=-2.d0*rhab/ab+(rhab-lnab)*r1*r1
c         accumulate psi
          br(ii)=br(ii)+gi
c         accumulate psi gradient
          drjj(1)=dgi*dr(1)
          drjj(2)=dgi*dr(2)
          drjj(3)=dgi*dr(3)
          gbrd(1,ii)=gbrd(1,ii)+drjj(1)
          gbrd(2,ii)=gbrd(2,ii)+drjj(2)
          gbrd(3,ii)=gbrd(3,ii)+drjj(3)
          gbr(1,ii,jj)=gbr(1,ii,jj)-drjj(1)
          gbr(2,ii,jj)=gbr(2,ii,jj)-drjj(2)
          gbr(3,ii,jj)=gbr(3,ii,jj)-drjj(3)

          if((r+rhoi).gt.rvdwj) then
c          ji contribution
           r1=1.d0/r
           r12=0.5d0*r1
           r24=r12*r12

           ap=r+rhoi
           am=r-rhoi
           rh1=1.d0/rvdwj
           rhr1=1.d0/ap
           aprh1=ap*rh1
           lnab=log(aprh1)

           gj=rh1-rhr1+r12*(0.5d0*am*(rhr1-rh1*aprh1)-lnab)

           dgj=rhr1*rhr1*(1.d0-0.25d0*am*r1*(1.d0+aprh1*aprh1))+
     .         rhoi*r24*(rhr1-rh1*aprh1)+
     .         r12*(r1*lnab-rhr1)
           dgj=dgj*r1
c          accumulate psi
           br(jj)=br(jj)+gj
c          accumulate psi gradient
           drjj(1)=dgj*dr(1)
           drjj(2)=dgj*dr(2)
           drjj(3)=dgj*dr(3)
           gbrd(1,jj)=gbrd(1,jj)-drjj(1)
           gbrd(2,jj)=gbrd(2,jj)-drjj(2)
           gbrd(3,jj)=gbrd(3,jj)-drjj(3)
           gbr(1,jj,ii)=gbr(1,jj,ii)+drjj(1)
           gbr(2,jj,ii)=gbr(2,jj,ii)+drjj(2)
           gbr(3,jj,ii)=gbr(3,jj,ii)+drjj(3)
          endif

c        ij overlap; ji do not overlap
         case(1)

          if((r+rhoj).gt.rvdwi) then
c          ij contribution
           r1=1.d0/r
           r12=0.5d0*r1
           r24=r12*r12

           ap=r+rhoj
           am=r-rhoj
           rh1=1.d0/rvdwi
           rhr1=1.d0/ap
           aprh1=ap*rh1
           lnab=log(aprh1)

           gi=rh1-rhr1+r12*(0.5d0*am*(rhr1-rh1*aprh1)-lnab)

           dgi=rhr1*rhr1*(1.d0-0.25d0*am*r1*(1.d0+aprh1*aprh1))+
     .         rhoj*r24*(rhr1-rh1*aprh1)+
     .         r12*(r1*lnab-rhr1)
           dgi=dgi*r1
c          accumulate psi
           br(ii)=br(ii)+gi
c          accumulate psi gradient
           drjj(1)=dgi*dr(1)
           drjj(2)=dgi*dr(2)
           drjj(3)=dgi*dr(3)
           gbrd(1,ii)=gbrd(1,ii)+drjj(1)
           gbrd(2,ii)=gbrd(2,ii)+drjj(2)
           gbrd(3,ii)=gbrd(3,ii)+drjj(3)
           gbr(1,ii,jj)=gbr(1,ii,jj)-drjj(1)
           gbr(2,ii,jj)=gbr(2,ii,jj)-drjj(2)
           gbr(3,ii,jj)=gbr(3,ii,jj)-drjj(3)
          endif

c         ji contribution
          ap=r+rhoi
          am=r-rhoi
          ab=ap*am
          rhab=rhoi/ab
          lnab=0.5d0*log(am/ap)*r1
          gj=rhab+lnab
          dgj=-2.d0*rhab/ab+(rhab-lnab)*r1*r1
c         accumulate psi
          br(jj)=br(jj)+gj
c         accumulate psi gradient
          drjj(1)=dgj*dr(1)
          drjj(2)=dgj*dr(2)
          drjj(3)=dgj*dr(3)
          gbrd(1,jj)=gbrd(1,jj)-drjj(1)
          gbrd(2,jj)=gbrd(2,jj)-drjj(2)
          gbrd(3,jj)=gbrd(3,jj)-drjj(3)
          gbr(1,jj,ii)=gbr(1,jj,ii)+drjj(1)
          gbr(2,jj,ii)=gbr(2,jj,ii)+drjj(2)
          gbr(3,jj,ii)=gbr(3,jj,ii)+drjj(3)

c        ij and ji overlap
         case(11)
c         overlaping spheres
          if((r+rhoj).gt.rvdwi) then
c          ij contribution
           r1=1.d0/r
           r12=0.5d0*r1
           r24=r12*r12

           ap=r+rhoj
           am=r-rhoj
           rh1=1.d0/rvdwi
           rhr1=1.d0/ap
           aprh1=ap*rh1
           lnab=log(aprh1)

           gi=rh1-rhr1+r12*(0.5d0*am*(rhr1-rh1*aprh1)-lnab)

           dgi=rhr1*rhr1*(1.d0-0.25d0*am*r1*(1.d0+aprh1*aprh1))+
     .         rhoj*r24*(rhr1-rh1*aprh1)+
     .         r12*(r1*lnab-rhr1)
           dgi=dgi*r1
c          accumulate psi
           br(ii)=br(ii)+gi
c          accumulate psi gradient
           drjj(1)=dgi*dr(1)
           drjj(2)=dgi*dr(2)
           drjj(3)=dgi*dr(3)
           gbrd(1,ii)=gbrd(1,ii)+drjj(1)
           gbrd(2,ii)=gbrd(2,ii)+drjj(2)
           gbrd(3,ii)=gbrd(3,ii)+drjj(3)
           gbr(1,ii,jj)=gbr(1,ii,jj)-drjj(1)
           gbr(2,ii,jj)=gbr(2,ii,jj)-drjj(2)
           gbr(3,ii,jj)=gbr(3,ii,jj)-drjj(3)
          endif

          if((r+rhoi).gt.rvdwj) then
c          ji contribution
           r1=1.d0/r
           r12=0.5d0*r1
           r24=r12*r12

           ap=r+rhoi
           am=r-rhoi
           rh1=1.d0/rvdwj
           rhr1=1.d0/ap
           aprh1=ap*rh1
           lnab=log(aprh1)

           gj=rh1-rhr1+r12*(0.5d0*am*(rhr1-rh1*aprh1)-lnab)

           dgj=rhr1*rhr1*(1.d0-0.25d0*am*r1*(1.d0+aprh1*aprh1))+
     .         rhoi*r24*(rhr1-rh1*aprh1)+
     .         r12*(r1*lnab-rhr1)
           dgj=dgj*r1
c          accumulate psi
           br(jj)=br(jj)+gj
c          accumulate psi gradient
           drjj(1)=dgj*dr(1)
           drjj(2)=dgj*dr(2)
           drjj(3)=dgj*dr(3)
           gbrd(1,jj)=gbrd(1,jj)-drjj(1)
           gbrd(2,jj)=gbrd(2,jj)-drjj(2)
           gbrd(3,jj)=gbrd(3,jj)-drjj(3)
           gbr(1,jj,ii)=gbr(1,jj,ii)+drjj(1)
           gbr(2,jj,ii)=gbr(2,jj,ii)+drjj(2)
           gbr(3,jj,ii)=gbr(3,jj,ii)+drjj(3)
          endif

        end select

       return
       end subroutine compute_psi

       subroutine compute_numsa(n,xyz,iat,rsas,sasai,grads)
       implicit none

       integer n,iat
       real*8 xyz(3,n),rsas,sasai,grads(3,n)

       integer ip,jj,nnj,nnk
       real*8 xyza(3),xyzp(3),sasap
       real*8 wr,wsa,drjj(3)
       integer :: nno
       real*8, allocatable :: grds(:,:)
       integer :: nni
       integer, allocatable :: grdi(:)

c       allocate space for the gradient storage
        nno=nnsas(iat)
        allocate(grds(3,nno))
        allocate(grdi(nno))

c       initialize storage
        sasai=0.d0

c       atomic position
        xyza(1)=xyz(1,iat)
        xyza(2)=xyz(2,iat)
        xyza(3)=xyz(3,iat)
c       radial atomic weight
        wr=wrp(iat)*gamsasa(iat)

c       loop over grid points
        do ip=1,nangsa
c        grid point position
         xyzp(1) = xyza(1) + rsas*grida(1,ip)
         xyzp(2) = xyza(2) + rsas*grida(2,ip)
         xyzp(3) = xyza(3) + rsas*grida(3,ip)
c        atomic surface function at the grid point
         call compute_w_sp(n,xyz,iat,nno,xyzp,sasap,grds,nni,grdi)

         if(sasap.gt.tolsesp) then
c         numerical quadrature weight
          wsa = grida(4,ip)*wr*sasap
c         accumulate the surface area
          sasai = sasai + wsa
c         accumulate the surface gradient
          do jj = 1, nni
           nnj = grdi(jj)
           drjj(1) = wsa*grds(1,jj)
           drjj(2) = wsa*grds(2,jj)
           drjj(3) = wsa*grds(3,jj)
           grads(1,iat)=grads(1,iat)+drjj(1)
           grads(2,iat)=grads(2,iat)+drjj(2)
           grads(3,iat)=grads(3,iat)+drjj(3)
           grads(1,nnj)=grads(1,nnj)-drjj(1)
           grads(2,nnj)=grads(2,nnj)-drjj(2)
           grads(3,nnj)=grads(3,nnj)-drjj(3)
          enddo
         endif
        enddo

       return
       end subroutine compute_numsa

       subroutine compute_w_sp(n,xyza,iat,nno,xyzp,sasap,grds,nni,grdi)
       implicit none

       integer n,iat,nno,nni
       real*8 xyza(3,n)
       real*8 xyzp(3)
       real*8 sasap
       real*8 grds(3,nno)
       integer grdi(nno)

       integer i,ia
       real*8 tj(3),tj2,sqtj
       real*8 uj,uj3,ah3uj2
       real*8 sasaij,dsasaij

c      initialize storage
       nni=0
       sasap=1.d0
       do i = 1, nno
         ia = nnlists(i,iat)
c        compute the distance to the atom
         tj(1)=xyzp(1)-xyza(1,ia)
         tj(2)=xyzp(2)-xyza(2,ia)
         tj(3)=xyzp(3)-xyza(3,ia)
         tj2=tj(1)*tj(1)+tj(2)*tj(2)+tj(3)*tj(3)
c        if within the outer cut-off compute
         if(tj2.lt.trj2(2,ia)) then
           if(tj2.le.trj2(1,ia)) then
            sasap=0.d0
            return
           else
            sqtj = sqrt(tj2)
            uj = sqtj - vdwsa(ia)
            ah3uj2 = ah3*uj*uj
            dsasaij = ah1+3.d0*ah3uj2
            sasaij =  ah0+(ah1+ah3uj2)*uj

c           accumulate the molecular surface
            sasap = sasap*sasaij
c           compute the gradient wrt the neighbor
            dsasaij = dsasaij/(sasaij*sqtj)
            nni=nni+1
            grdi(nni)=ia
            grds(1,nni) = dsasaij*tj(1)
            grds(2,nni) = dsasaij*tj(2)
            grds(3,nni) = dsasaij*tj(3)
           endif
c          check if the point is already buried
c          if(sasap.lt.tolsesp) return
         endif
       enddo

       return
       end subroutine compute_w_sp

       subroutine update_nnlist_gbsa(n12,xyz)
       implicit none

       integer n12
       real*8 xyz(3,n12)

       integer kk,i1,i2
       real*8 rcutn2,lrcut2,srcut2
       real*8 x,y,z,dr2
       integer ip,ip2,thrid,nproc
       integer omp_get_thread_num,omp_get_max_threads
       integer, allocatable :: npid(:)
       integer, allocatable :: plisttr(:,:,:)
       integer, allocatable :: nntmp(:)
       integer, allocatable :: nnls(:,:)

c      setup a pairlist and compute pair distances of all neighbors
c      within thresholds lrcut and srcut
       nproc=OMP_GET_MAX_THREADS()

       allocate(plisttr(3,ntpair,nproc),nnls(nat,nat))
       allocate(nntmp(nat),npid(nproc))
       npid = 0

       lrcut2 = lrcut*lrcut
       srcut2 = srcut*srcut

       nnsas=0
       nnlists=0
!$OMP PARALLEL PRIVATE ( i1,i2,x,y,z,dr2,ip,ip2,thrid,nntmp,nnls )
       ip=0
       ip2=0
       nntmp=0
       nnls=0
       thrid=omp_get_thread_num() + 1
!$OMP DO
       do kk=1,ntpair
         i1=ppind(1,kk)
         i2=ppind(2,kk)
         x=xyz(1,i1)-xyz(1,i2)
         y=xyz(2,i1)-xyz(2,i2)
         z=xyz(3,i1)-xyz(3,i2)
         dr2=x**2+y**2+z**2
         ddpair(2,kk)=x
         ddpair(3,kk)=y
         ddpair(4,kk)=z
         ddpair(1,kk)=sqrt(dr2)
         if(dr2.lt.lrcut2) then
          ip = ip + 1
          plisttr(1,ip,thrid)=i1
          plisttr(2,ip,thrid)=i2
          plisttr(3,ip,thrid)=kk
          if(dr2.lt.srcut2) then
           nntmp(i1) = nntmp(i1) + 1
           nntmp(i2) = nntmp(i2) + 1
           nnls(nntmp(i1),i1)=i2
           nnls(nntmp(i2),i2)=i1
          endif
         endif
       enddo
!$OMP END DO
       npid(thrid)=ip
!$OMP CRITICAL
       do i1=1,nat
        do i2=1,nntmp(i1)
         nnlists(nnsas(i1)+i2,i1)=nnls(i2,i1)
        enddo
        nnsas(i1)=nnsas(i1)+nntmp(i1)
       enddo
!$OMP END CRITICAL
!$OMP END PARALLEL

       nnrad=0
       do thrid=1,nproc
        do kk = nnrad+1,nnrad+npid(thrid)
         nnlistr(1:3,kk)=plisttr(1:3,kk-nnrad,thrid)
        enddo
        nnrad = nnrad + npid(thrid)
       enddo

       deallocate(nntmp,nnls)

       return
       end subroutine update_nnlist_gbsa

       subroutine update_dist_gbsa(n,xyz)
       implicit none

       integer n
       real*8 xyz(3,n)

       integer i1,i2,kk

!$OMP PARALLEL PRIVATE(i1,i2)
!$OMP DO
       do kk = 1, ntpair
         i1=ppind(1,kk)
         i2=ppind(2,kk)
         ddpair(2:4,kk)=xyz(1:3,i1)-xyz(1:3,i2)
         ddpair(1,kk)=sqrt(ddpair(2,kk)**2+
     .                     ddpair(3,kk)**2+
     .                     ddpair(4,kk)**2)
       enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end subroutine update_dist_gbsa

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end module gbobc
