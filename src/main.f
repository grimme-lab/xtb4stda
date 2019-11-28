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

      program XTBprog
      use gbobc, only: lgbsa,init_gbsa
      use ncoord, only: ncoord_erf
      implicit none

c symm stuff
c     include 'interface.f'

      integer n,nel,nbf,nao,nshell
      integer nopen,ncore,chrg,gsolvstate
      integer i,j,k,idum,nfrozh
      real*8 ,allocatable :: xyz (:,:)
      real*8 ,allocatable :: q   (:)
      real*8 ,allocatable :: qsh (:)
      real*8 ,allocatable :: chir(:)
      real*8 ,allocatable :: z   (:)
      real*8 ,allocatable :: cn  (:)
      real*8 ,allocatable :: sat (:)
      real*8 ,allocatable :: freq(:)
      real*8 ,allocatable :: wbo (:,:)
      real*8 ,allocatable :: g   (:,:)
      real*8 ,allocatable :: g1  (:,:)
      real*8 ,allocatable :: gdum(:,:)
      real*8 ,allocatable :: P   (:,:)
      integer,allocatable :: at(:)
      integer,allocatable :: frozh(:)
      real*8 xx(10),globpar(25)
      character*80 fname,fnv,arg(20),fnx
      character*128 atmp
      character*2 asym

c some method parameters
      real*8 gscal,kspd(6),tfac,kcn,zcnf,zqf,kcnsh(4)
      real*8 ken1,split,lshifta,lshift,wllscal,fpol,xbrad,xbdamp
      real*8 alphaj,d3a1,d3a2,d3s8,d3atm

      real*8 dum3,dum5,dum1,dum2,dum4,el,mowrcut,xtemp,egap,etot
      real*8 zero,t0,t1,w0,w1,empty2,acc,etot2,ipshift,eashift,g298
      real*8 qmdff_s6,one,two,d(3)
      real*8 :: etotl,etotr
      real*8,parameter :: step = 1e-5
      real*8,allocatable::dqdum(:,:,:)
      real x
      parameter (zero=0.0d0)
      parameter (one =1.0d0)
      parameter (two =2.0d0)

      logical ex,okbas,pop,dipol,wbocalc
      logical molden,rdchrg,moout,lindh
      logical epr,diff,test,vip,vea

c the basis setup stuff
      include 'ehtcommon.fh'
      include 'aoelementcommon.fh'
      include 'setcommon.fh'
      include 'spherecommon.fh'
      include 'scancommon.fh'
      include 'splitcommon.fh'
      include 'symcommon.fh'
      include 'stuff.fh'
      include 'd3common.fh' !GFN0
      include 'fixcommon.fh'

c OMP
      common /proc/ nproc
      integer nproc
      integer TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

                                                      call timing(t0,w0)
      write(*,*)
      write(*,'(7x,'' ______________________________________'')')
      write(*,'(7x,''|                                      |'')')
      write(*,'(7x,''|          ===================         |'')')
      write(*,'(7x,''|             xTB for sTDA             |'')')
      write(*,'(7x,''|          ===================         |'')')
      write(*,'(7x,''|               S. Grimme              |'')')
      write(*,'(7x,''|        Universitaet Bonn, MCTC       |'')')
      write(*,'(7x,''|         2015-19, Version 1.0         |'')')
      write(*,'(7x,''|     Wed Apr  3 09:14:55 CEST 2019    |'')')
      write(*,'(7x,''|______________________________________|'')')
      write(*,*)
      write(*,'(7x, '' This code comes without any warranty'')')
      write(*,'(7x, '' for non-commercial, academia use only.'')')
      write(*,'(7x, '' Preliminary test version.'')')
c     write(*,'(7x, '' Licensed to: Novartis Institutes'')')
c     write(*,'(7x, '' Test license to: Cresset, Cambridgeshire, UK'')')
c     write(*,'(7x, '' Test license to: R. Sure, BASF'')')
c     write(*,'(7x, '' Licensed to CreativeQuantum GmbH, Berlin'')')
      write(*,'(7x, '' Cite this work as:'')')
      write(*,'(7x,
     .'' S. Grimme  &  C. Bannwarth, JCP 145 (2016) 054103'')')
      write(*,'(7x, '' if GBSA is used additionally:'')')
      write(*,'(7x,
     .'' P. Shushkov, F. März, C. Bannwarth, S. Grimme, unpublished.''
     .)')
      write(*,*)
      write(*,'(7x, '' with help from'')')
      write(*,'(7x,
     .'' P. Shushkov, G. Brandenburg, S. Dohm, J. Pisarek,'')')
      write(*,'(7x,
     .'' F. März, M. Checinski, S. Ehrlich, S. Spicher, '')')
      write(*,'(7x,
     .'' P. Pracht, E. Caldeweyher, S. Ehlert, and C. Bannwarth.'')')
      write(*,*)
      write(*,'(7x, '' usage        :'')')
      write(*,'(7x, '' xtb4stda <coord_file> [options]'')')
      write(*,'(7x, '' where <coord_file> is a valid file of TM)'')')
      write(*,'(7x, '' (*coord, Bohr) or xmol (*xyz, Angstroem)'')')
      write(*,'(7x, '' format.'')')
      write(*,'(/7x,'' options:'')')
      write(*,'(7x,'' -chrg <int>   : molecular charge'')')
      write(*,'(7x,'' -uhf  <int>   : # of unpaired electrons(=2*S)'')')
      write(*,'(7x,
     .'' -nox          : skip second, extended part in sTDA-xTB'')')
      write(*,'(7x,'' -pop          : do population analysis'')')
      write(*,'(7x,''               : (default for -nox)'')')
      write(*,'(7x,'' -mowr <real>  : cut MO write above (def.3 Eh)'')')
      write(*,'(7x,'' -molden       : write formatted molden file'')')
      write(*,'(7x,'' -parx <file>  : read parameters for sTDA-xTB'')')
      write(*,'(7x,''                 calc (def:~/.param_stda2.xtb)'')')
      write(*,'(7x,
     .'' -parv <file>  : read parameters for vTB part in'')')
      write(*,'(7x,
     .''                 sTDA-xTB (def: ~/.param_stda1.xtb)'')')
      write(*,'(7x,'' -xtemp <real> : el. temp for xTB (def. 0 K)'')')
      write(*,'(7x,'' -etemp <real> : el. temp for GFN (def. 300 K)'')')
      write(*,'(7x,'' -fod          : calculate the FOD and write'')')
      write(*,'(7x,''                 molden.input file with'')')
      write(*,'(7x,''                 appropriate occupations'')')
      write(*,'(7x,''                 for plotting. Def. T=12500 K'')')
      write(*,*)
      write(*,'(7x,'' -acc <real>   : xTB accuracy (def. 1.0)'')')
      write(*,'(7x,'' -gbsa [string1] [string2]'')')
      write(*,'(7x,''                 use GBSA implicit solvent'')')
      write(*,'(7x,''                 for solvent [string1] and'')')
      write(*,'(7x,''                 solvation state [string2]='')')
      write(*,'(7x,''                 reference, bar1M (default=1M)'')')
      write(*,'(7x,'' additional GFN/opt/MD/siman options read from'')')
      write(*,'(7x,'' $XTB4STDAHOME/.xtb4stdarc or $set in input'')')
      write(*,*)
c     write(*,'(7x,''point charges (charge, xyz, atom symbol)'')')
c     write(*,'(7x,''are read from file <pcharge> (one line/charge)'')')
      write(*,'(7x,''spin and charge state information can be on:'')')
      write(*,'(7x,''<.CHRG> (charge) and <.UHF> (=nalpha-nbeta)'')')
      write(*,'(7x,''-uhf and -chrg override the file info.'')')
      write(*,'(7x,''useful machine settings:'')')
      write(*,'(7x,''setenv MKL_NUM_THREADS <NCORE_ON_YOUR_MACHINE>'')')
      write(*,'(7x,''setenv OMP_STACKSIZE 500m'')')
      write(*,'(7x,''limit stacksize unlimited'')')
      write(*,*)
      write(*,'(7x,''total energies in Eh, gaps/HL energies in eV'')')
      write(*,'(7x,''please read REVISION and HOWTO files carefully'')')

!$OMP PARALLEL PRIVATE(TID)
      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
         nproc = OMP_GET_NUM_THREADS()
         PRINT *, '============================='
         PRINT *, ' # OMP threads =', nproc
         PRINT *, '============================='
      END IF
!$OMP END PARALLEL

      xtb4stdahome=''
      call get_environment_variable('XTB4STDAHOME',xtb4stdahome)
      if(xtb4stdahome.eq.'') xtb4stdahome='~/'
      i=len(trim(xtb4stdahome))
      if(xtb4stdahome(i:i).ne.'/') xtb4stdahome(i+1:i+1)='/'
      write(*,*)
      write(*,*) 'xtb4stdahome directory:',trim(xtb4stdahome)
      write(*,*)

      arg=''
      do i=1, command_argument_count()
         call get_command_argument(i,arg(i))
      enddo
c coord filename
      fname = arg(1)
      inputname = fname

      pop     =.false.
      dipol   =.false.
      wbocalc =.false.
      rdchrg  =.false.
      epr     =.false.
      lmo     =.false.
      moout   =.false.
      molden  =.false.
      test    =.false.
      lgbsa   =.false.
      tsopt   =.false.
      moldfile=.false.
      vip     =.false.
      vea     =.false.
      lindh   =.false.
c which root in TS opt
      tsroot =0
c default parameter file names
      fnx=trim(xtb4stdahome) // trim('.param_stda2.xtb')
      fnv=''
c electronic temp
      etemp=-99.0d0 ! not set here
      xtemp=  0.0d0
c (virt)MO file write cut-off (Eh)
      mowrcut=3.
c opt level
      optlev =-99    ! not set here
      solvent='none'
      runtyp =-1
c frozen Hess
      nfrozh=0
c SCC acc. for SCC and grad
      acc   =1.0
c reference state for Gsolv (default 1 M gas/sol. = gsolv in cosmotherm)
      gsolvstate=0
c external code settings
      extcode=0
      extmode=0
c GFN version
      gfnver=-1 !GFN0
c C6 scaling of QMDFF dispersion to simulate solvent
      qmdff_s6=1.0d0
c Gtot for sdf output
      gtot=0
c default is no cube/ESP calc
      cube_cal=-1

      do i=2,20
         if(arg(i).ne.'')then
            write(*,*) 'argument ',i-1,':',trim(arg(i))
            if(index(arg(i),'-nodiff').ne.0)  runtyp=0
            if(index(arg(i),'-nox'   ).ne.0)  runtyp=0
            if(index(arg(i),'-pop'   ).ne.0)    pop = .true.
            if(index(arg(i),'-mold'  ).ne.0)  molden= .true.
            if(index(arg(i),'-dip'   ).ne.0)   dipol= .true.
            if(index(arg(i),'-wbo'   ).ne.0) wbocalc= .true.
            if(index(arg(i),'-rdchrg').ne.0)  rdchrg= .true.
            if(index(arg(i),'-esp'   ).ne.0)cube_cal= 2
            if(index(arg(i),'-stm'   ).ne.0)cube_cal= 3
            if(index(arg(i),'-lmo'   ).ne.0)then
                   lmo= .true.
                molden= .true.
                 moout= .true.
            endif
            if(index(arg(i),'-fod'   ).ne.0)then
               runtyp = 1
                  pop = .true.
                molden= .true.
                 xtemp=12500.
            endif
            if(index(arg(i),'-mowr'  ).ne.0)then
               call readl(arg(i+1),xx,j)
               mowrcut=xx(1)
            endif
            if(index(arg(i),'-xtemp' ).ne.0)then
               call readl(arg(i+1),xx,j)
               xtemp=xx(1)
            endif
            if(index(arg(i),'-etemp' ).ne.0)then
               call readl(arg(i+1),xx,j)
               etemp=xx(1)
            endif
            if(index(arg(i),'-acc' ).ne.0)then
               call readl(arg(i+1),xx,j)
               acc=xx(1)
               if(acc.lt.1.d-4.or.acc.gt.1000) acc=1.0d0
            endif
            if(index(arg(i),'-parv'  ).ne.0)then
               fnv=arg(i+1)
            endif
            if(index(arg(i),'-parx'  ).ne.0)then
               fnx=arg(i+1)
            endif
            if(index(arg(i),'-gbsa'  ).ne.0 .or.
     .         index(arg(i),'-g '  ).ne.0) then
               lgbsa=.true.
               atmp=arg(i+1)
               if(index(atmp,'-').eq.0.and.atmp(1:1).ne.' ')then
               solvent=arg(i+1)
               if(index(solvent,'none').ne.0) lgbsa=.false.
               endif
               atmp=arg(i+2)
               if(index(atmp,'reference').ne.0)gsolvstate=1
               if(index(atmp,'bar1M'    ).ne.0)gsolvstate=2
            endif
         endif
      enddo

      if(dipol.or.wbocalc) pop     =.true.
      if(molden          ) moout   =.true.   ! for vTB
      if(molden          ) moldfile=.true.   ! for GFN

      write(*,*)

c which elements are metals?
      call setmetal

c read molecule, # atoms first
      call rd0(fname,n)

      if(n.lt.1) stop 'no atoms!'
      if(n.gt.10000) stop 'too many atoms for this compilation!'

      allocate(at(n),xyz(3,n),wbo(n,n),q(n),z(n),
     .         cn(n),sat(n),g(3,n),frozh(n))

      call molbld(fname,arg,n,chrg,nopen,nel,at,z,xyz,lgbsa,
     .            nfrozh,frozh)

c     if(fit) acc=0.2 ! higher SCF accuracy during fit

      if(runtyp.lt.0) runtyp = 1 ! sTDA-xTB default if nothing is input
      if(runtyp.le.1.and.fnv.eq.'')
     .fnv=trim(xtb4stdahome) // trim('.param_stda1.xtb')

c----------------------------------------------------------------
c PARAMETERPARAMETERPARAMETERPARAMETERPARAMETERPARAMETERPARAMETER
c----------------------------------------------------------------
      inquire(file=fnv,exist=ex)
      if(.not.ex) then
      write(*,*) 'parameter file ',trim(fnv),' not found!'
      stop
      endif

      call rdelparam(1,fnv,globpar,gfnver)

      write(*,*)
      kspd(1:6)=globpar(1:6)
      wllscal  =globpar(7)
      gscal    =globpar(8)
      zcnf     =globpar(9)
      tfac     =globpar(10)
      kcn      =globpar(11)
      fpol     =globpar(12)
      ken1     =globpar(13)
      lshift   =globpar(14)
      lshifta  =globpar(15)
      split    =globpar(16)
      zqf      =globpar(17)
      write(*,*) '-----------------------------------------'
      write(*,*) '     charge density (VTB) calculation'
      write(*,*) '-----------------------------------------'
      write(*,*)
      write(*,'(''      method parameters     '')')
      write(*,'('' k(s)        :'',F8.4)') kspd(1)
      write(*,'('' k(p)        :'',F8.4)') kspd(2)
      write(*,'('' k(d)        :'',F8.4)') kspd(3)
      write(*,'('' k(f)        :'',F8.4)') kspd(4)
      write(*,'('' Tscal       :'',F8.4)') tfac
      write(*,'('' Gscal       :'',F8.4)') gscal
      write(*,'('' fpol        :'',F8.4)') fpol
      write(*,'('' Zcnf        :'',F8.4)') zcnf
      write(*,'('' Zqf         :'',F8.4)') zqf
      write(*,'('' kcn         :'',F8.4)') kcn
      write(*,'('' kEN1        :'',F8.4)') ken1
      write(*,'('' wllscal     :'',F8.4)') wllscal
      write(*,*)

c unrestricted atomic spin constants
c optimized scale factor for Mulliken spin densities
      call setwll_pbe (wllscal)

c init GBSA part
      if(lgbsa) then
         call init_gbsa(n,at,solvent,gsolvstate,temp_md)
      else
         write(*,*) '    -------------------------'
         write(*,*) '    ! NO SOLVENT MODEL USED !'
         write(*,*) '    -------------------------'
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                    calcs start here
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call xbasis0(n,at,nbf,nshell)

c EN charges and CN
      call iniqcn(n,nel,at,z,xyz,chrg,q,cn,ken1,runtyp,gfnver)

      write(*,*)
       write(*,
     . '(12x,''Z      q(EN)   CN      Cart. coordinates'')')
       do i=1,n
         write(*,'(i6,2x,a2,1x,f4.1,2(2x,f6.3),3f12.5)')
     .   i,asym(at(i)),z(i),q(i),cn(i),xyz(1:3,i)
       enddo
       write(*,*)
       call printbas(n,at)
       call xbasis (n,at,nbf,xyz,q,cn,
     .             zqf,zcnf,0.d0,okbas,diff,runtyp,gfnver) !GNF0
      if(.not.okbas) stop 'TB basis incomplete'

c dipole moment of PCs
      d = 0
      do i=1,n
         d(1)=d(1)+xyz(1,i)*q(i)
         d(2)=d(2)+xyz(2,i)*q(i)
         d(3)=d(3)+xyz(3,i)*q(i)
      enddo
      write(*,*)
      write(*,*)'dipole moment of classical point charges (au)'
      write(*,*)'    X       Y       Z   '
      write(*,'(3f9.4,''  total (Debye): '',f8.3,/)')
     . d(1:3),sqrt(d(1)**2+d(2)**2+d(3)**2)*2.5418

      if(extcode.gt.0.and.runtyp.eq.1.and.ntrans.gt.1)then
        call grdsym(xyz,n)
        call wrcoord( 6,n,xyz,at,0.0d0,'')
        call wrcoord(16,n,xyz,at,0.0d0,'coord')
        write(*,*)'symmetry adapted coordinates written to file <coord>'
        stop 'normal termination after symmetrysation'
      endif

         k=0
         do i=1,n
            do j=1,3
               k=k+1
               call random_number(dum1)
               rf(k)=xyz(j,i)+1.d-8*dum1
            enddo
         enddo

ccccccccccccccccccccccc
c make the VTB charges
ccccccccccccccccccccccc

      if(nopen.ne.0)then
      call uxtb(n,at,xyz,q,z,cn,sat,nel,nopen,nbf,el,wbo,
     .         kspd,gscal,tfac,kcn,lshift,lshifta,lshifta,
     .         mowrcut,.true.,dipol,.true.,.true.,moout,lmo,
     .         diff,1,etemp,epr)
      else
      call xtb(n,at,xyz,q,z,cn,sat,nel,nopen,nbf,el,wbo,
     .         kspd,gscal,tfac,kcn,lshift,
     .         mowrcut,.true.,dipol,.true.,.true.,moout,lmo,
     .         diff,1,etemp)
      endif
      if(runtyp.eq.0) goto 999

c----------------------------------------------------------------
c PARAMETERPARAMETERPARAMETERPARAMETERPARAMETERPARAMETERPARAMETER
c----------------------------------------------------------------

      inquire(file=fnx,exist=ex)
      if(.not.ex) then
      write(*,*) 'parameter file ',trim(fnx),' not found!'
      endif
      call rdelparam(0,fnx,globpar,gfnver)
      kspd(1:6)=globpar(1:6)
      wllscal  =globpar(7)
      gscal    =globpar(8)
      zcnf     =globpar(9)
      tfac     =globpar(10)
      kcn      =globpar(11)
      fpol     =globpar(12)
      ken1     =globpar(13)
      lshift   =globpar(14)
      lshifta  =globpar(15)
      split    =globpar(16)
      zqf      =globpar(17)

c unrestricted atomic spin constants
c optimized scale factor for OS exc.en.
      call setwll_pbe(wllscal)

c----------------------------------------------------------------
c PARAMETERPARAMETERPARAMETERPARAMETERPARAMETERPARAMETERPARAMETER
c----------------------------------------------------------------

      write(*,*)
      write(*,*) '--------------------------------------------'
      write(*,*) '         xTB calculation for sTDA...'
      write(*,*) '--------------------------------------------'
      write(*,*)
      write(*,*) 'reading parameter file ',trim(fnx)
      write(*,*)
      write(*,'(''      method parameters     '')')
      write(*,'('' k(s)        :'',F8.4)') kspd(1)
      write(*,'('' k(p)        :'',F8.4)') kspd(2)
      write(*,'('' k(d)        :'',F8.4)') kspd(3)
      write(*,'('' k(f)        :'',F8.4)') kspd(4)
      write(*,'('' k(R-V)      :'',F8.4)') kspd(5)
      write(*,'('' k(R-R)      :'',F8.4)') kspd(6)
      write(*,'('' Tscal       :'',F8.4)') tfac
      write(*,'('' Gscal       :'',F8.4)') gscal
      write(*,'('' fpol        :'',F8.4)') fpol
      write(*,'('' Zcnf        :'',F8.4)') zcnf
      write(*,'('' Zqf         :'',F8.4)') zqf
      write(*,'('' kcn         :'',F8.4)') kcn
      write(*,'('' Ryd Hsplit  :'',F8.4)') split
      write(*,'('' lshift(virt):'',F8.4)') lshift
      write(*,'('' lshift(OS)  :'',F8.4)') lshifta
      write(*,'('' wllscaling  :'',F8.4)') wllscal

      write(*,*)
      write(*,'('' mowrcut     :'',F6.3)') mowrcut
      write(*,*)

      call xbasis0(n,at,nbf,nshell)

c     read 'true' charges
      if(rdchrg)then
        allocate(chir(n))
        open(unit=11,file='charges.dat')
        do i=1,n
         read(11,*) q(i)
        enddo
c       call docm5(n,at,.false.,xyz,chir,q)
        deallocate(chir)
      endif

      call printbas(n,at)
      call xbasis (n,at,nbf,xyz,q,cn,
     .             zqf,zcnf,split,okbas,diff,runtyp)
      if(.not.okbas) stop 'TB basis incomplete'

      if(nopen.eq.0)then
      call xtb(n,at,xyz,q,z,cn,sat,nel,nopen,nbf,el,wbo,
     .         kspd,gscal,tfac,kcn,lshift,
     .         mowrcut,pop,dipol,wbocalc,molden,.true.,lmo,
     .         diff,0,xtemp)
      else
      call uxtb(n,at,xyz,q,z,cn,sat,nel,nopen,nbf,el,wbo,
     .         kspd,gscal,tfac,kcn,lshift,lshifta,lshifta,
     .         mowrcut,pop,dipol,wbocalc,molden,.true.,lmo,
     .         diff,0,xtemp,epr)
      endif


 999  write(*,*)
                                                      call timing(t1,w1)
                            write(*,'(''speedup'',f6.2)')(t1-t0)/(w1-w0)
                                        call prtime(6,t1-t0,w1-w0,'all')
      end

