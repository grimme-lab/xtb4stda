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

c filen: coordinate filename
c arg(): command line arguments
c n: # of atoms
c chrg: total charge
c spin: # S*2
c mel : # of electrons
c at(): their ordinal numbers
c z() : nuclear charges
c xyz(3,): cartesian coordinates

      subroutine molbld(filen,arg,n,chrg,spin,nel,at,z,xyz,lgbsa,
     .                  nfrozh,frozh)
      use gbobc, only: lsalt,ionst,ion_rad
      implicit none
      include 'setcommon.fh'
      include 'splitcommon.fh'
      include 'atmass.fh'
      include 'fixcommon.fh'

c input
      character*(*) filen
      character*80  arg(20)
      integer n,chrg,spin,nel,at(n),nfrozh,frozh(n)
      integer, allocatable :: iseed(:)
      real*8 z(n),xyz(3,n)
      logical lgbsa
c local
      integer i,j,ncore,rohf
      real*8 xx(50)
      real   xxx
      logical ex,inchrg,inspin

      chrg=0
      spin=0
      rohf=1 ! HS default
      pgroup='C1  '
      inchrg=.false.
      inspin=.false.
      molnameline=''
      commentline=''

      do i=1,20
         if(index(arg(i),'-chrg').ne.0.or.
     .      index(arg(i),'charge ').ne.0)then
            call readl(arg(i+1),xx,j)
            chrg=idint(xx(j))
            inchrg=.true.
         endif
         if(index(arg(i),'-uhf').ne.0)then
            call readl(arg(i+1),xx,j)
            spin=idint(xx(1))
            nalphabeta=spin
            inspin=.true.
         endif
         if(index(arg(i),'-LS').ne.0) rohf=0
      enddo

      inquire(file='.CHRG',exist=ex)
      if(ex.and.(.not.inchrg))then
         open(unit=1,file='.CHRG')
         read(1,*)chrg
         close(1)
         ichrg=int(chrg)
      endif

      inquire(file='.UHF',exist=ex)
      if(ex.and.(.not.inspin))then
         open(unit=1,file='.UHF')
         read(1,*)spin
         close(1)
         nalphabeta=int(spin)
      endif

      inquire(file='.SOLVENT',exist=ex)
      if(ex)then
        open(unit=1,file='.SOLVENT')
        read(1,'(a)') solvent
        write(*,*) 'file <.SOLVENT> found and consequently GBSA is on.'
        lgbsa=.true.
c       read(1,*) ionst  ! M=mol/l
c       read(1,*) ion_rad ! ang
        close(1)
        if(ionst.gt.0.d0) lsalt=.true.
      endif

c     inquire(file='.SYM',exist=ex)
c     if(ex)then
c       open(unit=1,file='.SYM')
c       read(1,'(a)') pgroup
c       close(1)
c     endif

      call rdcoord(filen,n,xyz,at,chrg)

      do i=1,n
         z(i) = at(i) - ncore( at(i) )
         if(at(i).gt.57.and.at(i).lt.72) z(i)=3   ! lanthanides without f are treated as La
         atmass(i)=ams(at(i))
      enddo

      call readset(filen,at,n,xyz,chrg,spin,lgbsa,
     .             nfrozh,frozh,rohf)

      nel = idint(sum(z))

      nel=nel-chrg

      if(spin.eq.0.and.mod(nel,2).ne.0) spin=1

      ichrg=chrg

      write(*,'('' name of molecule           :'',a)') trim(molnameline)
      write(*,'('' comment line               :'',a)') trim(commentline)
      write(*,'('' number of atoms            :'',i6)') n
      write(*,'('' number of electrons        :'',i6)') nel
      write(*,'('' charge                     :'',i3)') chrg
      write(*,'('' spin                       :'',f4.1)') 0.5*spin
      if(spin.gt.1)
     .write(*,'('' low(0)/high(1) spin ROHF   :'',i2  )') rohf

      if(mod(spin,2).ne.0.and.mod(nel,2).eq.0) stop 'Nel even, 2*S odd!'
      if(mod(spin,2).eq.0.and.mod(nel,2).ne.0) stop 'Nel odd, 2*S even!'

      if(samerand)then
        call random_seed(size=j)
        allocate(iseed(j))
        iseed(1:j)=41                    ! start random number generator for same sequenz
        do i=1,j
           iseed(i)=iseed(i)+j
        enddo
        call random_seed(put=iseed)
        deallocate(iseed)
      else
        call random_seed()
      endif
      call random_number(xxx)
      write(*,'('' first test random number   :'',f8.4)') xxx

      if(runtyp.le.3.and.extcode.gt.0)
     .stop 'SP/grad and extcode > 0 make no sense!'

      end




      subroutine readset(filename,at,n,xyz,chrg,spin,lgbsa,
     .                   nfrozh,frozh,rohf)
      use gbobc, only: lsalt,ionst,ion_rad
      implicit none
      integer chrg, spin, n, at(n), nfrozh, frozh(n),rohf
      real*8 xyz(3,n)
      character*(*) filename
      logical lgbsa
      include 'setcommon.fh'
      include 'atmass.fh'
      include 'fixcommon.fh'
      include 'scancommon.fh'
      include 'splitcommon.fh'
      include 'spherecommon.fh'
      include 'atomlistcommon.fh'

      character*128 line
      character*80  str
      character*2  asym
      real*8 xx(100),minmass,pi,phi,valijkl,ex_open_HS,ex_open_LS,rij(3)
      real*8 ra(3),rb(3)
      data pi/3.1415926535897932384626433832795029d0/
      integer nn,i,j,idum,jdum(n),k
      logical comline1,comline2,comline3,comline4,ex

      rdset=.false.
      comline1=.true.
      comline2=.true.
      comline3=.true.
      comline4=.true.

      inputname   = filename
      restart_md  = .false.
      honly       = .false.
      enan_siman  = .false.
      fit         = .false. ! write fit data in scf.f
      maxscciter  =250
      maxoptcycle =0        ! det. in ancopt routine if not read in
      time_md     =10.0
      temp_md     =298.15
      dump_md     = 250     ! scoord
      dump_md2    = 10      ! molden xyzfile
      skip_md     = 1000    ! mdopt, mdav
      accu_md     = 2.0
      tstep_md    =-1
      md_hmass    =1
      nvt_md      =.true.
      Tend_siman  =600.0
      ntemp_siman= 3
      accu_hess   = 0.25
      step_hess   = 0.005
      micro_opt   = 25
      hlow_opt    = 0.005
      maxdispl_opt= 1.000
      s6_opt      = 20.0
      fixfc       = 0.05
      natomsf     = 0
      fixed       = 0
      thermotemp  = 0
      thermotemp(1)=298.15
      thermo_sthr = 50.0d0
      broydamp    = 0.25
      path_nrun =3
      path_nopt =50
      path_anopt=3
      path_kpush= 0.05
      path_kpull=-0.15
      path_alp  = 0.80
      path_dens  =0.01
      cube_step  =0.50d0
      cube_pthr  =0.01d0
      ewin_conf  =4.0
      fcconstr   =0.05
      nconstr    =0
      iconstr    =0
      zconstr    =0
      atconstr   =0
      nscan      =0
      valscan    =0
      splitlist  =1
      iatf1      =0
      iatf2      =0
      desy       =0.1
      maxatdesy  =500
      mode_nscan =21
      mode_step  =0.5
      mode_vthr  =500.
      mode_updat =0.1
      mode_local =0
      mode_prj   =0
      ex_open_HS =0.3  ! exchange scaling factor for highspin
      ex_open_LS =-1.4 ! exchange scaling factor for lowspin
      samerand   =.false. ! initialize at each start the random number generator if .FALSE.
      mdrtrconstr=.false. ! not used
      check_rmsd =.true.
      sphere     =-1      ! off, 1=sphere, 2=ellipsoid, 4=logfermi
      boxr       =-1
      springexpo =2.0d0
      velodump   =.false.
      orcaexe    =''
      orcaline   =''
      stm_alp    =1.5
      stm_targ   =1.d-4
      stm_grid   =0.5
      stm_pot    =0.0
      stm_thr    =1.0
      metadyn_k  =0
      metadyn_a  =1.0
      metadyn_n  =0
      wallexpo   =10.
      natlist    =n ! ini common atomlist
      do i=1,n
         atlist(i)=i
      enddo

c those can be set by command line
      if(eTemp.lt.0) then
         comline1=.false.
         etemp=300.0
      endif
      if(index(solvent,'none').ne.0)then
         comline3=.false.
         solvent ='h2o'
      endif
      modflag     =0

c read defaults
      line=trim(XTB4STDAHOME) // trim('.xtb4stdarc')
      inquire(file=line,exist=ex)
      if(ex)then
      open(unit=1,file=line)
 100  read(1,'(a)',end=200)line
         if(index(line,'$set').ne.0)then
 110        read(1,'(a)',end=200)line
            if(index(line,'etemp ').ne.0)then
               call readl(line,xx,nn)
               if(.not.comline1)etemp=xx(1)
            endif
            if(index(line,'gbsa ').ne.0.and.(.not.comline3))then
               i=index(line,'gbsa')+5
               j=index(line,'#')
               if(index(line,'none').eq.0)then
                  solvent(1:j-i+1)=line(i:j-1)
                  lgbsa=.true.
               endif
            endif
            if(index(line,'samerand').ne.0)samerand=.true.
            if(index(line,'ion_st ').ne.0)then
               call readl(line,xx,nn)
               ionst=xx(1)
               if(ionst.gt.0.d0) lsalt=.true.
            endif
            if(index(line,'ion_rad ').ne.0)then
               call readl(line,xx,nn)
               ion_rad=xx(1)
            endif
            if(index(line,'desy ').ne.0)then
               call readl(line,xx,nn)
               desy=xx(1)
            endif
            if(index(line,'desymaxat ').ne.0)then
               call readl(line,xx,nn)
               maxatdesy=idint(xx(1))
            endif
            if(index(line,'ex_open_HS ').ne.0)then
               call readl(line,xx,nn)
               ex_open_HS=xx(1)
            endif
            goto 110
         endif
      goto 100
 200  continue
      close(1)
      endif

      open(unit=1,file=filename)
 1000 read(1,'(a)',end=2000)line
         if(index(line,'$set').ne.0)then
            rdset=.true.
 1100       read(1,'(a)',end=2000)line
            if(index(line,'chrg ').ne.0.or.
     .         index(line,'charge ').ne.0)then
               call readl(line,xx,nn)
               ichrg=int(xx(1))
               chrg =ichrg
            endif
            if(index(line,'uhf ').ne.0)then
               call readl(line,xx,nn)
               nalphabeta=int(xx(1))
               spin =nalphabeta
            endif
            if(index(line,'etemp ').ne.0)then
               call readl(line,xx,nn)
               if(.not.comline1)etemp=xx(1)
               modflag(7)=1
            endif
            if(index(line,'gbsa ').ne.0.and.(.not.comline3))then
               i=index(line,'gbsa')+5
               if(index(line,'none').eq.0)then
                  modflag(17)=1
                  solvent(1:20)=line(i:i+20)
                  lgbsa=.true.
               endif
            endif
            if(index(line,'ion_st ').ne.0)then
               call readl(line,xx,nn)
               ionst=xx(1)
               if(ionst.gt.0.d0) lsalt=.true.
               modflag(21)=1
            endif
            if(index(line,'ion_rad ').ne.0)then
               call readl(line,xx,nn)
               ion_rad=xx(1)
            endif
            if(index(line,'desy ').ne.0)then
               call readl(line,xx,nn)
               desy=xx(1)
               modflag(29)=1
            endif
            if(index(line,'desymaxat ').ne.0)then
               call readl(line,xx,nn)
               maxatdesy=idint(xx(1))
               modflag(30)=1
            endif
            if(index(line,'ex_open_HS ').ne.0)then
               call readl(line,xx,nn)
               ex_open_HS=xx(1)
            endif
            if(index(line,'ex_open_LS ').ne.0)then
               call readl(line,xx,nn)
               ex_open_LS=xx(1)
            endif
            goto 1100
         endif
      goto 1000
 2000 continue
      close(1)

 9999 continue

      if(tstep_md.lt.0) then
         minmass=10000.
         do i=1,n
            if(ams(at(i)).lt.minmass) minmass=ams(at(i))
         enddo
         tstep_md=(minmass/1.0079d0)**(1./3.)
      endif

      if(rohf.eq.0)then
         ex_open=ex_open_LS
      else
         ex_open=ex_open_HS
      endif

      if(sphere.eq.1)write(*,'(''spherical cavity R :'',f9.2)') boxr
      if(sphere.eq.2)write(*,'(''ellipsoid axis     :'',3f9.2)')rabco
      if(sphere.eq.3)write(*,'(''ellipsoid axis o/i :'',6f9.2)')
     .               rabco,rabci

      if(rdset)then
               write(*,*)
               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(*,*) '! reading options from coord file and    !'
               write(*,*) '! overriding <.xtb4stdarc>               !'
               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif

      end
